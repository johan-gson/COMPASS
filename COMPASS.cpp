#include <iostream>
#include <cmath>
#include <cfloat>
#include <random>
#include <fstream>
#include <string.h>
#include <stdexcept>
#include <omp.h>

#include "Inference.h"
#include "Tree.h"
#include "Scores.h"
#include "input.h"
#include <chrono>

//For debugging
#include <csignal>
#include <iostream>
#include <cstdlib>
#ifndef _MSC_VER
#include <execinfo.h>
#endif

int n_cells = 0;
int n_loci = 0;
int n_regions = 0;
std::size_t n_amplicons = 0;
std::vector<Cell> cells;
Data data;
Params parameters;

#ifndef _MSC_VER

void printCallstack() {
    // Obtain and print the backtrace
    void* array[10];
    size_t size = backtrace(array, 10);
    char** symbols = backtrace_symbols(array, size);

    std::cerr << "Backtrace:" << std::endl;
    for (size_t i = 0; i < size; i++) {
        std::cerr << symbols[i] << std::endl;
    }

    free(symbols);  // Remember to free the allocated memory

}

//For debugging of errors such as floating point errors
void signalHandler(int signal) {
    if (signal == SIGFPE) {
        std::cerr << "Caught a floating point exception (SIGFPE)." << std::endl;
        printCallstack();
        exit(EXIT_FAILURE);
    }
}


#endif


int main(int argc, char* argv[]){
    auto startTime = std::chrono::high_resolution_clock::now();

    // Register signal handler for SIGFPE
    #ifndef _MSC_VER
    std::signal(SIGFPE, signalHandler); //for debugging
    #endif

    init_params();
    parameters.verbose=false;
    // Read command line arguments
    std::string input_file{};
    std::string regionweights_file{};
    int n_chains=4;
    int chain_length=5000;
    int burn_in = -1;
    double temperature=10;
    double betabin_overdisp = parameters.omega_het;
    bool use_CNA=true;
    bool output_simplified = true;
    std::string start_tree;
    std::string output{};
    data.sex = "female";
    //parameters.verbose=true;
    for (int i=1;i<argc-1;i++){
        std::string argument{argv[i]};
        if (strcmp(argv[i],"-i")==0){
            input_file = argv[i+1];
        }
        else if (strcmp(argv[i],"--regionweights")==0){
            regionweights_file = argv[i+1];
        }
        else if (strcmp(argv[i],"--nchains")==0){
            n_chains=atoi(argv[i+1]);
        }
        else if (strcmp(argv[i],"--chainlength")==0){
            chain_length=atoi(argv[i+1]);
        }
        else if (strcmp(argv[i],"--burnin")==0){
            burn_in=atoi(argv[i+1]);
        }
        else if (strcmp(argv[i],"--temperature")==0){
            temperature=atoi(argv[i+1]);
        }
        else if (strcmp(argv[i],"--overdisp")==0){
            betabin_overdisp=atof(argv[i+1]);
        }
        else if (strcmp(argv[i],"--doubletrate")==0){
            parameters.doublet_rate=atof(argv[i+1]);
        }
        else if (strcmp(argv[i],"--dropoutrate")==0){ // mean of the prior dropout rate
            parameters.prior_dropoutrate_mean=atof(argv[i+1]);
        }
        else if (strcmp(argv[i],"--dropoutrate_concentration")==0){ // concentration parameter for the beta binomial distribution for the dropout rates (higher values: dropout rates will be closer to the mean)
            parameters.prior_dropoutrate_omega=atof(argv[i+1]);
        }
        else if (strcmp(argv[i],"--seqerror")==0){ // sequencing error rate
            parameters.sequencing_error_rate=atof(argv[i+1]);
        }
        else if (strcmp(argv[i],"--nodecost")==0){ // Penalty for adding nodes in the tree
            parameters.node_cost=atof(argv[i+1]);
        }
        else if (strcmp(argv[i],"--cnacost")==0){ // Penalty for adding CNA events in the tree
            parameters.CNA_cost=atof(argv[i+1]);
        }
        else if (strcmp(argv[i],"--lohcost")==0){ // Penalty for adding loh events in the tree
            parameters.LOH_cost=atof(argv[i+1]);
        }
        else if (strcmp(argv[i], "--mutNotAtRootPen") == 0) {
            parameters.mut_notAtRoot_cost = atof(argv[i++ + 1]);
        }
        else if (strcmp(argv[i],"-o")==0){
            output=argv[i+1];
        }
        else if (strcmp(argv[i],"-d")==0){
            if (strcmp(argv[i+1],"0")==0) parameters.use_doublets=false;
        }
        else if (strcmp(argv[i],"--CNA")==0){
            if (strcmp(argv[i+1],"0")==0) use_CNA=false;
        }
        else if (strcmp(argv[i],"--CNV")==0){
            if (strcmp(argv[i+1],"0")==0) use_CNA=false;
        }
        else if (strcmp(argv[i],"--filterregions")==0){
            if (strcmp(argv[i+1],"0")==0){
                parameters.filter_regions=false;
                parameters.filter_regions_CNLOH=false;
            }
        }
        else if (strcmp(argv[i],"--filterregionsCNLOH")==0){
            if (strcmp(argv[i+1],"0")==0){
                parameters.filter_regions_CNLOH=false;
            }
        }
        else if (strcmp(argv[i],"--verbose")==0){
            if (strcmp(argv[i+1],"1")==0) parameters.verbose=true;
        }
        else if (strcmp(argv[i],"--sex")==0){
           data.sex= std::string(argv[i+1]);
        }
        else if (strcmp(argv[i], "--prettyplot") == 0) {
            if (strcmp(argv[i + 1], "0") == 0) output_simplified = false;
        }
        else if (strcmp(argv[i], "--start_tree") == 0) {
            start_tree = argv[i + 1];
        }
        else if (argument.substr(0,1)=="-"){
            std::cout<<" Unrecognized argument: " <<argv[i]<<std::endl;
            throw std::invalid_argument("Invalid argument: "+ argument);
        }
    }
    if (input_file.size()==0){
        if (argc==2){
            input_file = argv[1];
            std::cout << "Will use "<<argv[1]<<" as input."<<std::endl;
        }
        else{
            throw std::invalid_argument("No input was provided. Please provide one with the -i option.");
        }
    }

    if (output.size()==0){
        std::cout << "No output name was provided. COMPASS will use the same basename as the input for the output." <<std::endl;
    }
    if (burn_in==-1){
        burn_in=chain_length/2;
    }
    std::cout << "Loading data...\n";
    load_CSV(input_file,regionweights_file,use_CNA); 
    std::cout << "Done...\n";

    parameters.omega_het = std::min(parameters.omega_het,betabin_overdisp);
    parameters.omega_het_indel = std::min(parameters.omega_het_indel,betabin_overdisp);

    // Get the name of the file, without directory
    std::string input_name = input_file;
    int name_start=0;
    for (int i=0;i<input_file.size();i++){
        if (input_file[i]=='/') name_start=i+1;
    }

    input_name=input_file.substr(name_start,input_file.size()-name_start);

	std::vector<double> results{};
    results.resize(n_chains);
    std::vector<Tree> best_trees{};
    best_trees.resize(n_chains);
    if (n_chains<omp_get_num_procs()) omp_set_num_threads(n_chains);
    else omp_set_num_threads(omp_get_num_procs());

    std::cout<<"Starting "<<std::to_string(n_chains)<< " MCMC chains in parallel"<<std::endl;
    #pragma omp parallel for
	for (int i=0;i<n_chains;i++){
		std::srand(i);
        Inference infer{ "",temperature,i,start_tree};
        best_trees[i] = infer.find_best_tree(use_CNA, chain_length, burn_in);
        results[i] = best_trees[i].log_score;
	}
    double best_score = -DBL_MAX;
    int best_score_index = -1;
    for (int i = 0; i < n_chains; i++) {
        if (best_score < results[i]) {
            best_score = results[i];
            best_score_index = i;
        }
    }
    auto& best_tree = best_trees[best_score_index];
    //It is a bit impractical, but the scores are deleted in the loop above - we create new ones here - will cost a bit of time
    Scores scores;
    best_tree.set_cache_scores(&scores);
    best_tree.resort_nodes_in_tree();
    if (output_simplified) best_tree.to_dot(output, true);
    else best_tree.to_dot(output, false);

    std::string gv_filename(output);
    if (output.size() <= 3 || (output.size() > 3 && output.substr(output.size() - 3) != ".gv")) {
        gv_filename = output + "tree.gv";
    }
    auto finishTime = std::chrono::high_resolution_clock::now();

    std::cout << "Completed in " << std::chrono::duration_cast<std::chrono::seconds> (finishTime - startTime).count() << 
        " seconds! The output was written to " << output << ". You can visualize the tree by running: dot -Tpng " << 
        gv_filename << " -o " << output << "output.png" << std::endl;

    return 0;
}
