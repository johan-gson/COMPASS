#include<cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include "Structures.h"
#include "Scores.h"
#include "input.h"
#include <cmath>

#ifdef _MSC_VER
// Different include paths in VS for gsl g++ the way I installed it
//#include <multimin/gsl_multimin.h>
//#include <fit/gsl_fit.h>
//#include <statistics/gsl_statistics.h>
//#include <statistics/gsl_statistics_double.h>
//#include <randist/gsl_randist.h>
#else
// Different include paths in linux for gsl g++
//#include <gsl/gsl_cdf.h>
//#include <gsl/gsl_sf_lngamma.h>
//#include <gsl_multimin.h>
#endif


//global variables
extern int n_cells;
extern int n_loci;
extern int n_regions;
extern std::vector<Cell> cells;
extern Data data;
extern Params parameters;


inline CellType parse_cell_type(const std::string& str) {
    if (str == "U") {
        return CT_U;
    }
    else if (str == "N") {
        return CT_N;
    }
    else if (str == "M") {
        return CT_M;
    }
    else {
        throw std::invalid_argument("Invalid cell type string: " + str);
    }
}

inline VariantType parse_variant_type(const std::string& str) {
    if (str == "UNKNOWN") {
        return VT_UNKNOWN;
    }
    else if (str == "GERMLINE") {
        return VT_GERMLINE;
    }
    else if (str == "SOMATIC") {
        return VT_SOMATIC;
    }
    else {
        throw std::invalid_argument("Invalid variant type string: " + str);
    }
}

inline CNType parse_cn_type(const std::string& str) {
    if (str == "UNKNOWN") {
        return CNT_UNKNOWN;
    }
    else if (str == "COPY_NEUTRAL") {
        return CNT_COPY_NEUTRAL;
    }
    else if (str == "CNV_DETECTED") {
        return CNT_CNV_DETECTED;
    }
    else {
        throw std::invalid_argument("Invalid region type string: " + str);
    }
}


void load_CSV(std::string base_name, std::string regionweights_file, bool use_CNA){
    std::ifstream file_variants(base_name+"_variants.csv");
    if(!file_variants.is_open()) throw std::runtime_error("Could not open variants file");
    // Read region counts (if use CNA)
    std::vector<std::vector<int>> region_counts{};
    data.region_to_name.clear();
    data.region_to_chromosome.clear();
    std::string line, val;
    if(use_CNA){
        std::ifstream file_region(base_name+"_regions.csv");
        if (!file_region.is_open()) throw std::runtime_error("Could not open region file");
        int region_index=0;
        while(std::getline(file_region,line,'\n')){
            std::stringstream ss(line);
            //index: chromosome and name of the region
            std::getline(ss, val, ',');
            int idx_split=0;
            std::vector<int> counts{};
            while (idx_split<val.size() && val[idx_split]!='_') idx_split++;
            if (idx_split<val.size()){
                // First column is chr_regionName
                data.region_to_name.push_back(val.substr(idx_split+1,val.size()-idx_split-1));
                data.region_to_chromosome.push_back(val.substr(0,idx_split));
            }
            else{
                // no chromosome and region_name were given: first column corresponds to the first cell
                data.region_to_name.push_back(std::to_string(region_index));
                data.region_to_chromosome.push_back(std::to_string(region_index));
                counts.push_back(stoi(val));
            }
            //second column is prior knowledge - is the region known to be copy neutral, or does it look like there is a CNV?
            std::getline(ss, val, ',');
            data.region_to_cn_type.push_back(parse_cn_type(val));
            //third column is theta
            std::getline(ss, val, ',');
            data.region_to_theta.push_back(atof(val.c_str()));

            while (std::getline(ss, val, ',')){
                counts.push_back(stoi(val));
            }
            region_counts.push_back(counts);
            region_index++;
        }
        n_regions = region_counts.size();
        file_region.close();
    }
    else{
        n_regions=0;
    }

    // Read variants
    data.locus_to_chromosome.clear();
    data.locus_to_position.clear();
    data.locus_to_reference.clear();
    data.locus_to_alternative.clear();
    data.variant_is_SNV.clear();
    data.locus_to_region.clear();
    data.region_to_loci.clear();
    data.locus_to_name.clear();
    data.locus_to_freq.clear();

    std::vector<std::vector<int>> ref_counts{};
    std::vector<std::vector<int>> alt_counts{};
    std::vector<std::vector<int>> genotypes{};
    std::vector<std::string> cell_names{};
    std::vector<int> ref_counts_variant{};
    std::vector<int> alt_counts_variant{};
    std::vector<int> genotypes_variant{};
    
    // First line: header
    std::getline(file_variants,line,'\n'); 
    std::stringstream header(line);
    std::vector<std::string> columns{};
    while (std::getline(header, val, ',')){
        columns.push_back(val);
        if (val != "CHR" && val!="POS" && val!="REF" && val!="ALT" && val!="REGION" && val!="NAME" && val!="FREQ" && val != "VARIANT_TYPE") {
            cell_names.push_back(val);
        }
    }


    // content
    while(std::getline(file_variants,line,'\n')){
        ref_counts_variant.clear();
        alt_counts_variant.clear();
        genotypes_variant.clear();
        std::stringstream ss(line);
        int column_count=0;
        while (std::getline(ss, val, ',')){
            
            if (columns[column_count]=="CHR"){
                data.locus_to_chromosome.push_back(val);
            }
            else if (columns[column_count]=="POS"){
                data.locus_to_position.push_back(stoi(val));
            }
            else if (columns[column_count]=="REF"){
                data.locus_to_reference.push_back(val);
            }
            else if (columns[column_count]=="ALT"){
                data.locus_to_alternative.push_back(val);
            }
            else if (columns[column_count]=="REGION"){
                int region_index = -1;
                for (int i=0;i<data.region_to_name.size();i++){
                    if (val==data.region_to_name[i]) region_index=i;
                }
                if (region_index==-1){
                    data.region_to_name.push_back(val);
                    data.region_to_chromosome.push_back(data.locus_to_chromosome[data.locus_to_chromosome.size()-1]);
                    region_index = n_regions;
                    n_regions++;
                }
                data.locus_to_region.push_back(region_index);
            }
            else if (columns[column_count]=="NAME"){
                data.locus_to_name.push_back(val);
            }
            else if (columns[column_count] == "FREQ") {
                data.locus_to_freq.push_back(stod(val));
            }
            else if (columns[column_count] == "VARIANT_TYPE") {
                data.locus_to_variant_type.push_back(parse_variant_type(val));
            }
            else{
                // for each cell, contains RO:AD:GT (:GT being optional)
                int pos = val.find(':');
                ref_counts_variant.push_back(stoi(val.substr(0,pos)));
                int pos2 = val.find(':',pos+1);
                if (pos2==std::string::npos){ // does not contain genotypes
                    alt_counts_variant.push_back(stoi(val.substr(pos+1,val.length()-pos-1)));
                    genotypes_variant.push_back(3);
                }
                else{
                    alt_counts_variant.push_back(stoi(val.substr(pos+1,pos2-pos-1)));
                    genotypes_variant.push_back(stoi(val.substr(pos2+1,1)));
                }
            }
            column_count++;
        }
        ref_counts.push_back(ref_counts_variant);
        alt_counts.push_back(alt_counts_variant);
        genotypes.push_back(genotypes_variant);  
    }
    file_variants.close();
    n_cells = ref_counts[0].size();
    n_loci = ref_counts.size();

    if (data.locus_to_reference.size()==0){
        data.variant_is_SNV = std::vector<bool>(n_loci,true);
    }
     else{
        for (int i=0;i<n_loci;i++){
            data.variant_is_SNV.push_back((data.locus_to_reference[i]=="A" || data.locus_to_reference[i]=="C" 
                                            || data.locus_to_reference[i]=="G" || data.locus_to_reference[i]=="T") 
                                &&(data.locus_to_alternative[i]=="A" || data.locus_to_alternative[i]=="C" 
                                            || data.locus_to_alternative[i]=="G" || data.locus_to_alternative[i]=="T"));
        }
    }

    if (data.locus_to_name.size()==0){
        if (data.locus_to_region.size()==n_loci){
            for (int i=0;i<n_loci;i++){
                data.locus_to_name.push_back(data.region_to_name[data.locus_to_region[i]]);
            }
        }
        else{
            for (int i=0;i<n_loci;i++) data.locus_to_name.push_back(std::to_string(i));
        }
    }

    if (use_CNA && data.region_to_name.size()==0){
        for (int k=0;k<n_regions;k++) data.region_to_name.push_back(std::to_string(k));
    }
    

    // In case no mapping from variants to regions were provided
    if (data.locus_to_region.size()==0){
        if (use_CNA) throw std::invalid_argument("Missing region information for variants. When using CNAs, the variants file must contain a column indicating to which region (amplicon or gene) each variant belongs.");
        for (int i=0;i<n_loci;i++){
            data.locus_to_region.push_back(i);
            data.region_to_name.push_back(data.locus_to_name[i]);
            data.region_is_reliable.push_back(false);
            data.region_to_chromosome.push_back(" ");
        }
        n_regions = n_loci;
    }
    // Map region index to loci index
    data.region_to_loci.resize(n_regions);
    for (int i=0;i<n_loci;i++){
        data.region_to_loci[data.locus_to_region[i]].push_back(i);
    }

    // In case no variant frequencies were provided
    if (data.locus_to_freq.size()==0) data.locus_to_freq = std::vector<double>(n_loci,0.0);

    //read cell metadata
    std::vector<CellType> cellTypes;
    std::ifstream file_cellmeta(base_name + "_cell_metadata.csv");
    if (file_cellmeta) {
        std::size_t i = 0;
        std::string line;
        std::getline(file_cellmeta, line);//skip header
        while (std::getline(file_cellmeta, line)) {
            // Skip empty lines
            if (line.empty()) continue;

            std::istringstream iss(line);
            std::string cellID, cellTypeStr;

            // Split the line by the comma delimiter
            if (std::getline(iss, cellID, ',') && std::getline(iss, cellTypeStr, ',')) {
                try {
                    // Convert the cell type string to enum
                    CellType ct = parse_cell_type(cellTypeStr);
                    // Store the results in the vectors
                    if (cellID != cell_names[i]) {
                        throw std::runtime_error(std::string("The cell ids in ") + base_name + "_cell_metadata.csv does not match that of the variant file. They must be in the same order!");
                    }
                    cellTypes.push_back(ct);
                    ++i;
                }
                catch (const std::invalid_argument& e) {
                    std::cerr << e.what() << " Skipping line: " << line << std::endl;
                    // Optionally, handle the error (e.g., log it) and continue
                }
            }
            else {
                std::cerr << "Error parsing line: " << line << std::endl;
                continue; // Skip malformed lines
            }
        }
        file_cellmeta.close();
        if (cellTypes.size() != cell_names.size()) {
            throw std::runtime_error(std::string("The cell ids in ") + base_name + "_cell_metadata.csv does not match that of the variant file. The number of cells differ.");
        }
    }
    else {
        std::cout << "No cell metadata file found - setting all cell types to unknown.";
        cellTypes = std::vector<CellType>(n_cells, CellType::CT_U);
    }
    // store by cell
    cells.clear();
    cells.reserve(n_cells);
    for (int j=0;j<n_cells;j++){
        cells.push_back(Cell{});
        cells[j].ref_counts.reserve(n_loci);
        cells[j].alt_counts.reserve(n_loci);
        cells[j].genotypes.reserve(n_loci);
        cells[j].GQ.reserve(n_loci);
        for (int i=0;i<n_loci;i++){
            cells[j].ref_counts.push_back(ref_counts[i][j]);
            cells[j].alt_counts.push_back(alt_counts[i][j]);
            cells[j].genotypes.push_back(genotypes[i][j]);
        }
        cells[j].name = cell_names[j];
        int total_count=0;
        if (use_CNA){
            for (int i=0;i<n_regions;i++){
                cells[j].region_counts.push_back(region_counts[i][j]);
                total_count+=region_counts[i][j];
            }
            cells[j].total_counts=total_count;
        }
        cells[j].cell_type = cellTypes[j];
    }
    data.predetermined_region_weights.clear();
    if (use_CNA){
        // Read region weights, if they are given as input
        if (regionweights_file!=""){
            data.predetermined_region_weights = std::vector<double>(n_regions,-1);
            std::string line, val;
            std::ifstream file_region(regionweights_file);
            if (!file_region.is_open()) throw std::runtime_error("Could not open file with region weights.");
            int region_index=0;
            while(std::getline(file_region,line,'\n')){
                int idx_comma=0;
                while (line[idx_comma]!=',') idx_comma++;
                std::string region_name = line.substr(0,idx_comma);
                double weight = stod(line.substr(idx_comma+1,line.size()-idx_comma-1));
                for (int k=0;k<n_regions;k++){
                    if (data.region_to_name[k] == region_name){
                        data.predetermined_region_weights[k] = weight;
                    }
                }
            }
        }
        // Filter regions with insufficient coverage
        if (parameters.filter_regions) filter_regions();
        else data.region_is_reliable = std::vector<bool>(n_regions,true);
    }
    else data.region_is_reliable = std::vector<bool>(n_regions,false);

    //fit_region_thetas(); //fit the thetas to the counts from each region - the dispersion varies a lot between regions
}





void filter_regions(){
    // Filter out regions for which many cells have 0 (or almost 0) reads
    double threshold = 1.0 / n_regions / 15.0;
    data.region_is_reliable.clear();
    bool regions_filtered=false;
    for (int k=0;k<n_regions;k++){
        int count_cells_below_threshold=0;
        double mean=0;
        for (int j=0;j<n_cells;j++){
            if (1.0*cells[j].region_counts[k] / cells[j].total_counts <= threshold) count_cells_below_threshold++;
            mean+= 1.0*cells[j].region_counts[k] / cells[j].total_counts / n_cells;
        }
        bool is_reliable = ((1.0*count_cells_below_threshold/n_cells <= 0.04) && (mean>=0.2/n_regions));
        if (data.predetermined_region_weights.size()>k && data.predetermined_region_weights[k]==-1) is_reliable=false;
        data.region_is_reliable.push_back(is_reliable);
        regions_filtered = regions_filtered || ((1.0*count_cells_below_threshold/n_cells > 0.04) || (mean<0.2/n_regions));
    }
    if (regions_filtered){
        std::cout<<"The following regions are excluded from the CNA inference because their coverage is too low: ";
        for (int k=0;k<n_regions;k++){
            if (!data.region_is_reliable[k]) std::cout<<data.region_to_name[k]<<",";
        }
        std::cout<<std::endl;
    }
    else{
        std::cout<<"No regions filtered "<<std::endl;
    }
    
}

void init_params(){
    parameters.sequencing_error_rate=0.02;
	parameters.omega_hom=50.0;
	parameters.omega_het=8.0;
	parameters.sequencing_error_rate_indel=0.06; // higher error rate and dispersion for indels because the allelic calls are less reliable
	parameters.omega_hom_indel = 15.0;
	parameters.omega_het_indel = 4.0;

    // The dropout rates are inferred for each SNV, using a beta distribution as prior
    parameters.prior_dropoutrate_mean=0.05; 
    parameters.prior_dropoutrate_omega=100; // concentration parameter for the beta distribution (higher values: force the dropout rates to be close to the prior mean)

	//parameters.theta=6.0;
	parameters.doublet_rate=0.08;

    parameters.use_doublets=true;
    parameters.filter_regions=true; // if true, only allow CNAs on regions with sufficient coverage
    parameters.filter_regions_CNLOH=true; // if filter_regions_CNLOH is false but filter_regions is true, CNLOH will still be possible in the regions with low coverage (but not gains and losses).
    parameters.verbose=true;

    // Tree prior
    parameters.node_cost=1.0; //higher value will result in fewer nodes
    parameters.CNA_cost=85.0; //Higher values will result in fewer CNAs
    parameters.LOH_cost=85.0; // CNLOH and losses resulting in a LOH have a higher penalty (the penalty for such an event is the sum of CNA_cost and LOH_cost)
    parameters.mut_notAtRoot_cost=10;
    parameters.mut_notAtRoot_freq_cost=100000; // Penalty for not placing SNVs present in the 1000G database at the root.
}

/*
#include <cmath>
#include <vector>

// Data structure to pass to the minimizer
struct data {
    const std::vector<int>* counts;
    const std::vector<double>* means;
};

// Log-likelihood function for the Negative Binomial distribution
double neg_log_likelihood(const gsl_vector* v, void* params) {
    double theta = gsl_vector_get(v, 0);
    struct data* d = (struct data*)params;
    const std::vector<int>& counts = *(d->counts);
    const std::vector<double>& means = *(d->means);

    double neg_log_likelihood = 0.0;
    for (std::size_t i = 0; i < counts.size(); ++i) {
        double mu = means[i];
        int k = counts[i];
        neg_log_likelihood -= std::lgamma(k + theta) - std::lgamma(theta) - std::lgamma(k + 1)
            + theta * std::log(theta / (theta + mu))
            + k * std::log(mu / (theta + mu));
    }
    return neg_log_likelihood;
}

// Function to fit theta using a minimizer
double fit_theta_gsl(const std::vector<int>& counts, const std::vector<double>& means) {
    // Initial guess for theta
    double theta_initial = 1.0;

    // Set up the GSL minimizer
    gsl_vector* x = gsl_vector_alloc(1);
    gsl_vector_set(x, 0, theta_initial);

    gsl_multimin_function minex_func;
    minex_func.n = 1;
    minex_func.f = neg_log_likelihood;
    struct data d = { &counts, &means };
    minex_func.params = &d;

    // Use the Nelder-Mead simplex algorithm for minimization
    const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc(T, 1);
    gsl_multimin_fminimizer_set(s, &minex_func, x, gsl_vector_alloc(1));

    // Set the precision and tolerance
    double tol = 1e-6;
    int iter = 0;
    int status;

    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status) break;

        double size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, tol);

    } while (status == GSL_CONTINUE && iter < 10000);

    double theta = gsl_vector_get(s->x, 0);

    gsl_vector_free(x);
    gsl_multimin_fminimizer_free(s);

    return theta;
}





// Fit theta parameter for negative binomial using maximum likelihood estimation
double fit_theta(std::size_t region, const std::vector<std::size_t>& normal_cell_indices, double region_proportion) {
    double theta = parameters.theta; // initial guess
    double epsilon = 1e-6; // convergence threshold
    double step_size = 0.01; // gradient descent step size
    double best_log_likelihood = -std::numeric_limits<double>::infinity();
    double last_log_likelihood = -std::numeric_limits<double>::infinity();
    double best_theta = theta;

    // Iterate to optimize theta
    for (int iter = 0; iter < 1000; ++iter) {
        double current_log_likelihood = 0.0;

        // Iterate through each cell
        for (std::size_t i = 0; i < normal_cell_indices.size(); ++i) {
            double mu = cells[normal_cell_indices[i]].total_counts * region_proportion;
            double log_theta = std::log(theta);
            double log_mu = std::log(mu);
            double log_theta_mu = std::log(theta + mu);
            double count = cells[normal_cell_indices[i]].region_counts[region];

            // Log-PMF for negative binomial
            double log_pmf = std::lgamma(count + theta) - std::lgamma(theta) - std::lgamma(count + 1)
                + theta * log_theta - (count + theta) * log_theta_mu
                + count * log_mu - std::lgamma(count + 1);

            current_log_likelihood += log_pmf;
        }

        if (current_log_likelihood > best_log_likelihood) {
            double improvement = current_log_likelihood - best_log_likelihood;
            best_log_likelihood = current_log_likelihood;
            best_theta = theta;
            if (improvement < epsilon) {
                break; //goal reached
            }
        }
        else {
            //break; // Convergence //we can't stop here
        }

        // Gradient descent step to update theta
        if (iter == 0) {
            theta += step_size * 10; //we move in a random direction since we don't have a gradient - we will get one in the next step
        }
        else {
            theta += step_size * (current_log_likelihood - last_log_likelihood);
        }
        last_log_likelihood = current_log_likelihood;
        
    }

    return best_theta;
}

void write_thetas(std::string filename) {
    std::ofstream of(filename);
    std::cout << "region\ttheta\n";
    for (std::size_t i = 0; i < data.region_to_theta.size(); ++i) {
        std::cout << data.region_to_name[i] << '\t' << data.region_to_theta[i] << '\n';
    }
}

void fit_region_thetas() {
    //init with the standard theta, in case of no/too few normal cells. It works pretty poorly though.
    data.region_to_theta = std::vector<double>(data.region_to_name.size(), parameters.theta);
    //get normal cells
    std::vector<std::size_t> normal_cell_indices;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        if (cells[i].cell_type == CellType::CT_N) {
            normal_cell_indices.push_back(i);
        }
    }
    if (normal_cell_indices.size() < 50) {
        std::cout << "Too few normal cells to estimate thetas for the counts distributions per region: " << normal_cell_indices.size() <<".\n" <<
            "At least 50 is required.\n" <<
            "Using the standard value - this leads to worse performance for inferral of CNA.\n";
        return;
    }

    //get region weights
    std::vector<double> region_proportions;//the fraction of counts that end up in each region, should sum to 1
    if (data.predetermined_region_weights.size() > 0) {
        region_proportions = data.predetermined_region_weights;
    } 
    else {
        //calculate it from the normal cells. First get the row sums for all regions
        std::vector<double> rs(0.0, data.region_to_name.size());
        double tot_sum = 0.0;
        for (std::size_t r = 0; r < data.region_to_name.size(); ++r) {
            for (std::size_t c = 0; c < normal_cell_indices.size(); ++r) {
                rs[r] += cells[normal_cell_indices[c]].region_counts[r];
            }
            tot_sum += rs[r];
        }
        region_proportions = rs;
        for (std::size_t r = 0; r < data.region_to_name.size(); ++r) {
            region_proportions[r] /= tot_sum;
        }
    }
    

    std::cout << "Fitting thetas...";

    for (std::size_t r = 0; r < data.region_to_theta.size(); ++r) {
        //What we get here is a vector of counts per cell and the means, i.e., mu,
        //i.e., the expected counts for that cell given its total counts and region
        //we calculate the mean as tot_counts*region_proportion
        std::vector<int> counts(normal_cell_indices.size(),0);
        std::vector<double> means(normal_cell_indices.size(), 0);
        for (std::size_t c = 0; c < normal_cell_indices.size(); ++c) {
            counts[c] = cells[normal_cell_indices[c]].region_counts[r];
            means[c] = cells[normal_cell_indices[c]].total_counts * region_proportions[r];
        }
        
        data.region_to_theta[r] = fit_theta_gsl(counts, means);
        

        //data.region_to_theta[r] = fit_theta(r, normal_cell_indices, region_proportions[r]);
    }
    std::cout << "done\n";
}
*/