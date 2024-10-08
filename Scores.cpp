#include <cmath>
#include <vector>
#include <cfloat>
#include <string>
#include <iostream>
#include <random>

#include "Scores.h"
#include "Node.h"
#include "Structures.h"


//global variables
extern int n_loci;
extern int n_cells;
extern int n_regions;
extern size_t n_amplicons;

extern std::vector<Cell> cells;
extern Data data;
extern Params parameters;



Scores::Scores(){
    std::map<long int,std::vector<std::vector<double>>> cache_likelihood_allelecounts_cells{};
    std::map<int,std::map<long int,std::vector<double>>> cache_dropoutrate_cellscores{};
    std::map<long int,std::vector<double>> cache_dropoutscores{};

    std::map<int,std::map<long int,std::vector<double>>> cache_dropoutrate_dropoutsref{};
    std::map<int,std::map<long int,std::vector<double>>> cache_dropoutrate_dropoutsalt{};
    std::map<int,double> cache_n_choose_k{};
    count_cache=0;
}

double Scores::log_sum_exp(const std::vector<double>& terms){
    double max = -DBL_MAX;
    for (int i=0;i<terms.size();i++){
        if (terms[i]>max) max=terms[i];
    }

    double sum=0;
    for (int i=0;i<terms.size();i++){
        sum+= std::exp(terms[i]-max);
    }
    return std::log(sum) + max;
}

std::vector<double> Scores::log_sum_exp_vector(const std::vector<std::vector<double>>& terms){
    std::vector<double> maxs(n_cells,-DBL_MAX);
    for (int i=0;i<terms.size();i++){
    	for (int j=0;j<n_cells;j++){
    		if (terms[i][j]>maxs[j]) maxs[j]=terms[i][j];
    	}
    }

    std::vector<double> sums(n_cells,0.0);
    for (int i=0;i<terms.size();i++){
    	for (int j=0;j<n_cells;j++){
    		sums[j]+= std::exp(terms[i][j]-maxs[j]);
    	}
    }
    for (int j=0;j<n_cells;j++) sums[j] = std::log(sums[j]) + maxs[j];
    return sums;
}

double Scores::log_n_choose_k(int n, int k){
    int hash = n+1000*k;
    if (cache_n_choose_k.count(hash)) return cache_n_choose_k[hash];
    double res=0;
    for (int i=0; i<k;i++){
        res+= std::log(n-i);
        res-= std::log(i+1);
    }
    cache_n_choose_k[hash] = res;
    return res;
}


std::vector<double> Scores::compute_SNV_loglikelihoods(int c_ref,int c_alt,int locus, double dropout_rate_ref, double dropout_rate_alt){
    // compute the SNV log-likelihood for one locus for all of the cells
    // c_ref and c_alt are the copy numbers of each allele in the genotype
    
    // If homozygous, the copy number of the only allele is irrelevant for the allelic proportion
    if ((c_ref == 0) && (c_alt == 0)) { //double loss, we expect a balanced noise, so set equal proportions the same way as wt
        c_ref = 1;
        c_alt = 1;
    }
    else if (c_ref==0) c_alt=1;
    else if (c_alt==0) c_ref=1;

    long int hash = c_ref+ 20*c_alt + 400*locus ;
    
    int discretized_dropout_rate_ref = std::round(dropout_rate_ref*1000);
    int discretized_dropout_rate_alt = std::round(dropout_rate_alt*1000);
    int hash_dropoutrates = discretized_dropout_rate_ref + 1000 * discretized_dropout_rate_alt;
    // if we already computed the scores for this dropout rate and this combination of locus, c_ref and c_alt.
    if (cache_dropoutrate_cellscores.count(hash_dropoutrates) && cache_dropoutrate_cellscores[hash_dropoutrates].count(hash)){
        return cache_dropoutrate_cellscores[hash_dropoutrates][hash];
    }
    dropout_rate_ref = 1.0*discretized_dropout_rate_ref/1000.0;
    dropout_rate_alt = 1.0*discretized_dropout_rate_alt/1000.0;
    std::vector<std::vector<double>> likelihood_alleles_cells{};
    likelihood_alleles_cells.resize((c_ref+1)*(c_alt+1)-1);

    std::vector<double> dropoutscores{};
    dropoutscores.resize(n_cells);
    int idx=0;
    for (int k=0;k<=c_ref;k++){
        for (int l=0;l<=c_alt;l++){
            if (k==0 && l==0) continue;
            long int hash_dropout = k+20*l + 400 * locus;
            if (cache_dropoutscores.count(hash_dropout)){
                dropoutscores = cache_dropoutscores[hash_dropout];                
            }
            else{ // does not depend on dropout rate so does not have to be computed often
                std::vector<double> seq_error_rates(n_cells,0.0);
                double f;
                double likelihood_dropout;
                double omega;
                double eps1; // sequencing errors from ref to alt (always rare)
                double eps2; // sequencing errors from alt to ref (seem quite common for indels)
                if (data.variant_is_SNV[locus]){
                    eps1 = parameters.sequencing_error_rate;
                    eps2 = parameters.sequencing_error_rate;
                    if (k==0 || l==0) omega = parameters.omega_hom;
                    else omega = parameters.omega_het;
                }
                else{
                    eps1 = parameters.sequencing_error_rate;
                    eps2 = parameters.sequencing_error_rate_indel;
                    if (k==0 || l==0) omega = parameters.omega_hom_indel;
                    else omega = parameters.omega_het_indel;
                }

                f = 1.0*l/(k+l) * (1-eps2) + 1.0*k/(k+l) * eps1; // frequency of the alt nucleotide
                for (int j=0;j<n_cells;j++){
                    int ref_count= cells[j].ref_counts[locus];
                    int alt_count=cells[j].alt_counts[locus];
                    likelihood_dropout = std::lgamma(alt_count + omega*f) + std::lgamma(ref_count + omega*(1-f)) - std::lgamma(alt_count+ref_count + omega)
                    -std::lgamma(omega*f) - std::lgamma(omega*(1-f)) + std::lgamma(omega);
                    //the term (a+r) choose a does not depend on the parameters, so does not need to be included here. + log_n_choose_k(ref_count+alt_count,alt_count)

                    dropoutscores[j] = likelihood_dropout;
                }
                cache_dropoutscores[hash_dropout] = dropoutscores;
            }
            // add dropout probability
            likelihood_alleles_cells[idx].resize(n_cells);
            double dropout_prob = log_n_choose_k(c_ref,k) + log_n_choose_k(c_alt,l) 
                                            +(c_ref-k) * std::log(dropout_rate_ref) + k * std::log(1-dropout_rate_ref)
                                            + (c_alt -l)*std::log(dropout_rate_alt) + l * std::log(1-dropout_rate_alt);
            // Cannot have a dropout of all the alleles
            double all_dropout_prob = std::pow(dropout_rate_ref,c_ref) * std::pow(dropout_rate_alt,c_alt);
            dropout_prob-= std::log(1-all_dropout_prob);
            for (int j=0;j<n_cells;j++){
                likelihood_alleles_cells[idx][j] = dropoutscores[j] + dropout_prob;
            }
            idx+=1;
        }
    }
    // Compute the score for each cell by summing over each dropout combination
    std::vector<double> scores_cells = log_sum_exp_vector(likelihood_alleles_cells);
    std::vector<double> dropoutsref(n_cells,0.0);
    std::vector<double> dropoutsalt(n_cells,0.0);
    double config_prob;

    // Compute the of average number of dropouts at this locus that occured in each cell
    // This is used for optimizing the dropout rate
    idx=0;
    for (int k=0;k<=c_ref;k++){
        for (int l=0;l<=c_alt;l++){
            if (k==0 && l==0) continue;
            for (int j=0;j<n_cells;j++){
                config_prob = std::exp(likelihood_alleles_cells[idx][j]-scores_cells[j]);
                dropoutsref[j]+=1.0*(c_ref-k) * config_prob;
                dropoutsalt[j]+=1.0*(c_alt-l) * config_prob;
            }
            idx+=1;
        }
    }

    cache_dropoutrate_cellscores[hash_dropoutrates][hash] = scores_cells;
    cache_dropoutrate_dropoutsref[hash_dropoutrates][hash] = dropoutsref;
    cache_dropoutrate_dropoutsalt[hash_dropoutrates][hash] = dropoutsalt;
    count_cache++;
    return scores_cells;
}

std::vector<double> Scores::compute_CNA_loglikelihoods(int region, const std::vector<double>& region_probabilities, double exp_cn, const std::vector<double>& normalization_factors){
    std::vector<double> region_proportions;
    for (std::size_t i = 0; i < region_probabilities.size(); ++i) {
        region_proportions.push_back(region_probabilities[i] * exp_cn / normalization_factors[i]);
    }
    
    // Compute the likelihood of the read count in the region, based on the negative binomial distribution (Gamma-Poisson)
    // theta is the scale parameter for the Gamma distribution.
    //note that there is a risk that 64 bits are not enough here if we have 6 or more samples, so a bit risky, then in theory 
    //two calculations could give the same hash, although unlikely
    //uint64_t hash = region;
    //const int max_discr_prop_val = 2000;
    //uint64_t multFactor = n_regions;
    //for (std::size_t i = 0; i < region_probabilities.size(); ++i) {
    //    uint64_t discretized_region_proportion = std::round(region_proportions[i] * max_discr_prop_val);
    //    hash += multFactor * discretized_region_proportion;
    //    multFactor *= max_discr_prop_val;
    //}
    //hmm, since we have region in the hash, we only need to use the first probability - the others will be the same as the first was when calculated last time
    const double max_discr_prop_val = 1000.0;//this number has a big impact on performance - if higher, it takes MUCH longer to run
    uint64_t discretized_region_proportion = uint64_t(std::round(region_proportions[0] * max_discr_prop_val));
    uint64_t hash = region + n_regions*discretized_region_proportion;
    if (cache_cnalikelihood_cells.count(hash)) return cache_cnalikelihood_cells[hash];

    //There are two options for how to calculate this: 
    // 1) We use the region counts
    // 2) We use the collective log likelihood of all amplicons (amplicon counts) within the region
    //Regardless of which, the caching on region proportion still works
    // The strategy is to use option 2 if we have amplicons and more than one amplicon for the region, otherwise option 1.

    std::vector<double> cnv_loglikelihoods(n_cells, 0.0);

    if (n_amplicons > 0 && data.region_to_amplicons[region].size() > 1) {
        //option 2 - amplicon counts
        std::vector<double> amplicon_proportions(data.sample_ids.size(), 0.0);
        for (auto amp : data.region_to_amplicons[region]) {
            const auto& thetas = data.amplicon_to_theta[amp]; //parameters.theta;
            //precalculate the possible amplicon_proportions for all cells
            for (std::size_t sample = 0; sample < data.sample_ids.size(); ++sample) {
                amplicon_proportions[sample] = data.amplicon_weights[amp][sample] * exp_cn / normalization_factors[sample];
            }
            //double amplicon_proportion = data.amplicon_weights[amp] * exp_cn / normalization_factor;
            for (int j = 0; j < n_cells; j++) {
                std::size_t sample = cells[j].sample;
                double theta = thetas[sample];
                double expected_read_count_amplicon = cells[j].total_counts * amplicon_proportions[sample];
                cnv_loglikelihoods[j] += std::lgamma(cells[j].amplicon_counts[amp] + theta - 1) + theta * std::log(theta / (theta + expected_read_count_amplicon))
                    + cells[j].amplicon_counts[amp] * std::log(expected_read_count_amplicon / (expected_read_count_amplicon + theta));
            }
        }
    }
    else {
        //option 1 - region counts
        const auto& thetas = data.region_to_theta[region]; //parameters.theta;
        for (int j = 0; j < n_cells; j++) {
            std::size_t sample = cells[j].sample;
            double theta = thetas[sample];
            double expected_read_count_region = cells[j].total_counts * region_proportions[sample];
            cnv_loglikelihoods[j] = std::lgamma(cells[j].region_counts[region] + theta - 1) + theta * std::log(theta / (theta + expected_read_count_region))
                + cells[j].region_counts[region] * std::log(expected_read_count_region / (expected_read_count_region + theta));
        }
    }
    cache_cnalikelihood_cells[hash] = cnv_loglikelihoods;
    count_cache++;
    return cnv_loglikelihoods;
}



std::vector<double> Scores::get_dropoutref_counts_genotype(int c_ref,int c_alt, int locus,double dropout_rate_ref,double dropout_rate_alt){
    if ((c_ref == 0) && (c_alt == 0)) { //double loss, we expect a balanced noise, so set equal proportions the same way as wt
        c_ref = 1;
        c_alt = 1;
    }
    else if (c_ref==0) c_alt=1;
    else if (c_alt==0) c_ref=1;
    if (c_ref==0 || c_alt==0){
        dropout_rate_ref=0.1;
        dropout_rate_alt=0.1;
    }
    long int hash = c_ref+ 20*c_alt + 400*locus;

    int discretized_dropout_rate_ref = std::round(dropout_rate_ref*1000);
    int discretized_dropout_rate_alt = std::round(dropout_rate_alt*1000);
    int hash_dropoutrate = discretized_dropout_rate_ref + discretized_dropout_rate_alt * 1000;
    // if we already computed the scores for this dropout rate and this combination of locus, c_ref and c_alt.
    if (cache_dropoutrate_dropoutsref.count(hash_dropoutrate) && cache_dropoutrate_dropoutsref[hash_dropoutrate].count(hash)){
        return cache_dropoutrate_dropoutsref[hash_dropoutrate][hash];
    }
    else{
        auto temp = compute_SNV_loglikelihoods(c_ref,c_alt,locus,dropout_rate_ref,dropout_rate_alt);
        return cache_dropoutrate_dropoutsref[hash_dropoutrate][hash];
    }
}

std::vector<double> Scores::get_dropoutalt_counts_genotype(int c_ref,int c_alt, int locus,double dropout_rate_ref,double dropout_rate_alt){
    if ((c_ref == 0) && (c_alt == 0)) { //double loss, we expect a balanced noise, so set equal proportions the same way as wt
        c_ref = 1;
        c_alt = 1;
    }
    else if (c_ref==0) c_alt=1;
    else if (c_alt==0) c_ref=1;
    if (c_ref==0 || c_alt==0){
        dropout_rate_ref=0.1;
        dropout_rate_alt=0.1;
    }
    long int hash = c_ref+ 20*c_alt + 400*locus ;

    int discretized_dropout_rate_ref = std::round(dropout_rate_ref*1000);
    int discretized_dropout_rate_alt = std::round(dropout_rate_alt*1000);
    int hash_dropoutrate = discretized_dropout_rate_ref + discretized_dropout_rate_alt * 1000;
    // if we already computed the scores for this dropout rate and this combination of locus, c_ref and c_alt.
    if (cache_dropoutrate_dropoutsalt.count(hash_dropoutrate) && cache_dropoutrate_dropoutsalt[hash_dropoutrate].count(hash)){
        return cache_dropoutrate_dropoutsalt[hash_dropoutrate][hash];
    }
    else{
        auto temp = compute_SNV_loglikelihoods(c_ref,c_alt,locus,dropout_rate_ref,dropout_rate_alt);
        return cache_dropoutrate_dropoutsalt[hash_dropoutrate][hash];
    }
}


void Scores::clear_cache(){
    cache_dropoutrate_cellscores.clear();
    cache_dropoutrate_dropoutsalt.clear();
    cache_dropoutrate_dropoutsref.clear();

    cache_cnalikelihood_cells.clear();

    count_cache=0;
}

void Scores::clear_cache_if_too_large(){
    // To avoid using too much memory
    if (count_cache>8000) clear_cache();
    //if (count_cache > 12000) clear_cache();
}





