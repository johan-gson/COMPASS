// TestStuffVS.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <cmath>
#include <cfloat>
#include <random>
#include <fstream>
#include <string.h>
#include <stdexcept>
#include <omp.h>

#include "../Inference.h"
#include "../Tree.h"
#include "../Scores.h"
#include "../input.h"
//For debugging
//#include <csignal>
//#include <iostream>
//#include <execinfo.h>
//#include <cstdlib>

int n_cells = 0;
int n_loci = 0;
int n_regions = 0;
std::size_t n_amplicons = 0;
std::vector<Cell> cells;
Data data;
Params parameters;


Tree build_debug_tree_X1() {
    //this code follows the same strategy as when reading a tree from file, but instead builds a hardcoded one
    //that is pefectly suited for testing something
    auto* cache_scores = new Scores(); //This will cause a memory leak, but it is just for test so it doesn't matter
    Tree tree(cache_scores, true, true);

    tree.add_node_no_randomization(0);

    //this is hardcoded to specific sample 
    tree.get_node(1)->add_mutation(0);
    tree.get_node(1)->add_mutation(1);

    tree.init_debug_tree(true);

    return tree;
}

Tree build_debug_tree_X2() {
    //this code follows the same strategy as when reading a tree from file, but instead builds a hardcoded one
    //that is pefectly suited for testing something
    auto* cache_scores = new Scores(); //This will cause a memory leak, but it is just for test so it doesn't matter
    Tree tree(cache_scores, true, true);

    tree.add_node_no_randomization(0);
    tree.add_node_no_randomization(0);

    //this is hardcoded to specific sample 
    tree.get_node(1)->add_mutation(0);
    tree.get_node(2)->add_mutation(1);

    tree.init_debug_tree(true);

    return tree;
}


Tree build_debug_tree_1() {
    //this code follows the same strategy as when reading a tree from file, but instead builds a hardcoded one
    //that is pefectly suited for testing something
    auto* cache_scores = new Scores(); //This will cause a memory leak, but it is just for test so it doesn't matter
    Tree tree(cache_scores, true, true);

    //We create the following tree
    //0 parent - all germline variants
    //1 child 1 connected to parent - has "NOTCH1:chr9:136496196:CAG/C" and EGR2:chr10:62813488:G/T
    //2 grandchild 1 - has CHD2:chr15:93009179:T/C
    //3 grandchild 2 - has LOH 36 (chr17:178000:10981999)

    //tree.add_node_no_randomization(0);
    //skip  these for now
    //tree.add_node_no_randomization(1);
    //tree.add_node_no_randomization(1);

    //this is hardcoded to specific sample 
    tree.get_node(0)->add_mutation(0); //NOTCH1:chr9:136496196:CAG/C
    tree.get_node(0)->add_mutation(1); //EGR2:chr10:62813488:G/T
    //tree.get_node(2)->add_mutation(3); //CHD2:chr15:93009179:T/C
    //tree.get_node(3)->add_mutation(0); //NOTCH1:chr9:136496196:CAG/C
    //place these in the root to not have them affect anything
    tree.get_node(0)->add_mutation(3); //CHD2:chr15:93009179:T/C
    //all germline
    for (std::size_t i = 3; i < data.locus_to_name.size(); ++i) {
        tree.get_node(0)->add_mutation(i);
    }
    //tree.get_node(1)->add_CNA(std::make_tuple(36, -1, std::vector<int>{1, 0}), false);

    tree.init_debug_tree(true);

    return tree;
}

Tree build_debug_tree_with_bial() {
    Tree tree = build_debug_tree_1();
    tree.add_node_no_randomization(0);
    tree.get_node(1)->add_CNA(std::make_tuple(18, -2, std::vector<int>{1}), false);

    tree.init_debug_tree(true);

    return tree;
}

Tree build_debug_tree_with_prior_init_CNV() {
    Tree tree = build_debug_tree_1();
    tree.add_node_no_randomization(0);

    int region = 24;
    int type = -1;

    std::vector<int> alleles{};
    //at a 33% chance, initiate the alleles according to the prior supplied as input by the user.
    //otherwise just randomize
    bool bInit = true;
    double alleleHastings = bInit ? 1.0 / 3.0 : 2.0 / 3.0;

    for (int i = 0; i < data.region_to_loci[region].size(); i++) {
        int allele = std::rand() % 2; //0: CNA affects ref allele; 1: CNA affects alt allele
        alleleHastings = alleleHastings * 2.0;
        if (bInit) {
            switch (data.locus_to_cna_allele_prior[data.region_to_loci[region][i]]) {
            case CNAP_DECREASED:
                if (type <= 0) {
                    //we are adding a loss of some kind. We then assume that a decreased CNV means that the mutation is lost
                    allele = 1;
                }
                else {
                    //this is a gain, we then assume that the ref allele is increased (since AF is decreased)
                    allele = 0;
                }
                alleleHastings = alleleHastings / 2;//reduce hastings, no randomization
                break;
            case CNAP_INCREASED:
                if (type <= 0) {
                    //we are adding a loss of some kind. We then assume that an increased CNV means that the reference is lost
                    allele = 0;
                }
                else {
                    //this is a gain, we then assume that the alt allele is increased (since AF is increased)
                    allele = 1;
                }
                alleleHastings = alleleHastings / 2;//reduce hastings, no randomization
                break;
            case CNAP_UNKNOWN:
                //Do nothing, we use the random value as is, and don't change the hastings ratio
                break;
            }
        }

        // Can only gain/lose lose one allele if we had at least one copy of it !
        //if (allele == 0 && nodes[parent]->get_n_ref_allele(data.region_to_loci[region][i]) == 0) allele = 1;
        //if (allele == 1 && nodes[parent]->get_n_alt_allele(data.region_to_loci[region][i]) == 0) allele = 0;
        alleles.push_back(allele);
    }


    tree.get_node(1)->add_CNA(std::make_tuple(region, type, alleles), false);

    tree.init_debug_tree(true);

    return tree;
}


Tree build_debug_tree_1_with_CNV_node() {
    //this code follows the same strategy as when reading a tree from file, but instead builds a hardcoded one
    //that is pefectly suited for testing something
    auto* cache_scores = new Scores(); //This will cause a memory leak, but it is just for test so it doesn't matter
    Tree tree(cache_scores, true, true);

    //We create the following tree
    //0 parent - all germline variants
    //1 child 1 connected to parent - has "NOTCH1:chr9:136496197:AGG/A" and EGR2:chr10:62813488:G/T
    //2 grandchild 1 - has CHD2:chr15:93009179:T/C
    //3 grandchild 2 - has NOTCH1:chr9:136496196:CAG/C

    tree.add_node_no_randomization(0);

    //this is hardcoded to specific sample 
    tree.get_node(1)->add_mutation(1); //NOTCH1:chr9:136496197:AGG/A
    tree.get_node(1)->add_mutation(2); //EGR2:chr10:62813488:G/T
    //place these in the root to not have them affect anything
    tree.get_node(0)->add_mutation(3); //CHD2:chr15:93009179:T/C
    tree.get_node(0)->add_mutation(0); //NOTCH1:chr9:136496196:CAG/C

    //all germline
    for (std::size_t i = 4; i < data.locus_to_name.size(); ++i) {
        tree.get_node(0)->add_mutation(i);
    }

    //now insert an extra node with CNV
    tree.add_node_no_randomization(1);
    tree.get_node(2)->add_CNA(std::tuple<int, int, std::vector<int>>(1, -1, std::vector<int>{0, 1, 1, 1, 0}), false);//chr1 loss (region 1), SNPs 4,5,6,7,8, with ref/alt 0,1,1,1,0 

    tree.init_debug_tree(true);

    return tree;
}

Tree build_debug_tree_1_with_CNLOH_after_LOH() {
    auto tree = build_debug_tree_1_with_CNV_node();
    tree.add_node_no_randomization(2);
    tree.get_node(3)->add_CNA(std::tuple<int, int, std::vector<int>>(1, 0, std::vector<int>{0, 1, 1, 1, 0}), false);//chr1 loss (region 1), SNPs 4,5,6,7,8, with ref/alt 0,1,1,1,0 
    tree.init_debug_tree(true);
    return tree;
}

Tree build_debug_tree_1_chr9_LOH() {
    auto tree = build_debug_tree_1();
    tree.add_node_no_randomization(1);
    tree.get_node(2)->add_CNA(std::tuple<int, int, std::vector<int>>(18, -1, std::vector<int>{0}), false);//This is the place with the double loss
    tree.init_debug_tree(true);
    return tree;
}

Tree build_debug_tree_1_chr9_double_loss() {
    auto tree = build_debug_tree_1();
    tree.add_node_no_randomization(1);
    tree.get_node(2)->add_CNA(std::tuple<int, int, std::vector<int>>(18, -2, std::vector<int>{0}), false);//This is the place with the double loss
    tree.init_debug_tree(true);
    return tree;
}

//This tree should give a penalty in the prior since the CNLOH removes a previous mutation
Tree build_debug_tree_1_with_CNLOH_node_that_removes_mut() {
    //this code follows the same strategy as when reading a tree from file, but instead builds a hardcoded one
    //that is pefectly suited for testing something
    auto* cache_scores = new Scores(); //This will cause a memory leak, but it is just for test so it doesn't matter
    Tree tree(cache_scores, true, true);

    //We create the following tree
    //0 parent - all germline variants
    //1 child 1 connected to parent - has "NOTCH1:chr9:136496197:AGG/A" and EGR2:chr10:62813488:G/T
    //2 grandchild 1 - has CHD2:chr15:93009179:T/C
    //3 grandchild 2 - has NOTCH1:chr9:136496196:CAG/C

    tree.add_node_no_randomization(0);

    //this is hardcoded to specific sample 
    tree.get_node(1)->add_mutation(1); //NOTCH1:chr9:136496197:AGG/A
    tree.get_node(1)->add_mutation(2); //EGR2:chr10:62813488:G/T
    //place these in the root to not have them affect anything
    tree.get_node(0)->add_mutation(3); //CHD2:chr15:93009179:T/C
    tree.get_node(0)->add_mutation(0); //NOTCH1:chr9:136496196:CAG/C

    //all germline
    for (std::size_t i = 4; i < data.locus_to_name.size(); ++i) {
        tree.get_node(0)->add_mutation(i);
    }

    //now insert an extra node with CNV
    tree.add_node_no_randomization(1);
    tree.get_node(2)->add_CNA(std::tuple<int, int, std::vector<int>>(22, 0, std::vector<int>{1, 1}), false);//chr1 CHLOH (region 22), SNPs 0,1 with ref/alt 1,1 

    tree.init_debug_tree(true);

    return tree;
}



Tree build_debug_tree_CN_1() {
    auto* cache_scores = new Scores(); //This will cause a memory leak, but it is just for test so it doesn't matter
    Tree tree(cache_scores, true, true);

    //We create the following tree
    //0 parent - all germline variants
    //1 child 1 connected to parent - no variants

    tree.add_node_no_randomization(0);
    tree.get_node(0)->add_mutation(0);
    tree.get_node(0)->add_mutation(1);
    tree.get_node(1)->add_mutation(2);

    tree.init_debug_tree(true);

    return tree;
}

//std::vector<std::vector<double>> rootScores;
//std::vector<std::vector<double>> nodeScores;
std::vector<std::vector<double>> savedScoresPerCellAndLoci;//temp debug

int main(int argc, char* argv[])
{
    init_params();
    parameters.verbose = false;
    // Read command line arguments
    std::string input_file{};
    std::string regionweights_file{};
    int n_chains = 4;
    int chain_length = 5000;
    int burn_in = -1;
    double temperature = 10;
    double betabin_overdisp = parameters.omega_het;
    bool use_CNA = true;
    bool output_simplified = true;
    std::string output{};
    data.sex = "female";
    //parameters.verbose=true;
    for (int i = 1; i < argc - 1; i++) {
        std::string argument{ argv[i] };
        if (strcmp(argv[i], "-i") == 0) {
            input_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "--regionweights") == 0) {
            regionweights_file = argv[i + 1];
        }
        else if (strcmp(argv[i], "--nchains") == 0) {
            n_chains = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "--chainlength") == 0) {
            chain_length = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "--burnin") == 0) {
            burn_in = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "--temperature") == 0) {
            temperature = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i], "--overdisp") == 0) {
            betabin_overdisp = atof(argv[i + 1]);
        }
        else if (strcmp(argv[i], "--doubletrate") == 0) {
            parameters.doublet_rate = atof(argv[i + 1]);
        }
        else if (strcmp(argv[i], "--dropoutrate") == 0) { // mean of the prior dropout rate
            parameters.prior_dropoutrate_mean = atof(argv[i + 1]);
        }
        else if (strcmp(argv[i], "--dropoutrate_concentration") == 0) { // concentration parameter for the beta binomial distribution for the dropout rates (higher values: dropout rates will be closer to the mean)
            parameters.prior_dropoutrate_omega = atof(argv[i + 1]);
        }
        else if (strcmp(argv[i], "--seqerror") == 0) { // sequencing error rate
            parameters.sequencing_error_rate = atof(argv[i + 1]);
        }
        else if (strcmp(argv[i], "--nodecost") == 0) { // Penalty for adding nodes in the tree
            parameters.node_cost = atof(argv[i + 1]);
        }
        else if (strcmp(argv[i], "--cnacost") == 0) { // Penalty for adding CNA events in the tree
            parameters.CNA_cost = atof(argv[i + 1]);
        }
        else if (strcmp(argv[i], "--lohcost") == 0) { // Penalty for adding loh events in the tree
            parameters.LOH_cost = atof(argv[i + 1]);
        }
        else if (strcmp(argv[i], "--mutNotAtRootPen") == 0) {
            parameters.mut_notAtRoot_cost = atof(argv[i++ + 1]);
        }
        else if (strcmp(argv[i], "-o") == 0) {
            output = argv[i + 1];
        }
        else if (strcmp(argv[i], "-d") == 0) {
            if (strcmp(argv[i + 1], "0") == 0) parameters.use_doublets = false;
        }
        else if (strcmp(argv[i], "--CNA") == 0) {
            if (strcmp(argv[i + 1], "0") == 0) use_CNA = false;
        }
        else if (strcmp(argv[i], "--CNV") == 0) {
            if (strcmp(argv[i + 1], "0") == 0) use_CNA = false;
        }
        else if (strcmp(argv[i], "--filterregions") == 0) {
            if (strcmp(argv[i + 1], "0") == 0) {
                parameters.filter_regions = false;
                parameters.filter_regions_CNLOH = false;
            }
        }
        else if (strcmp(argv[i], "--filterregionsCNLOH") == 0) {
            if (strcmp(argv[i + 1], "0") == 0) {
                parameters.filter_regions_CNLOH = false;
            }
        }
        else if (strcmp(argv[i], "--verbose") == 0) {
            if (strcmp(argv[i + 1], "1") == 0) parameters.verbose = true;
        }
        else if (strcmp(argv[i], "--sex") == 0) {
            data.sex = std::string(argv[i + 1]);
        }
        else if (strcmp(argv[i], "--prettyplot") == 0) {
            if (strcmp(argv[i + 1], "0") == 0) output_simplified = false;
        }
        else if (argument.substr(0, 1) == "-") {
            std::cout << " Unrecognized argument: " << argv[i] << std::endl;
            throw std::invalid_argument("Invalid argument: " + argument);
        }
    }
    if (input_file.size() == 0) {
        if (argc == 2) {
            input_file = argv[1];
            std::cout << "Will use " << argv[1] << " as input." << std::endl;
        }
        else {
            throw std::invalid_argument("No input was provided. Please provide one with the -i option.");
        }
    }

    if (output.size() == 0) {
        std::cout << "No output name was provided. COMPASS will use the same basename as the input for the output." << std::endl;
    }
    if (burn_in == -1) {
        burn_in = chain_length / 2;
    }

    load_CSV(input_file, regionweights_file, use_CNA);

    parameters.omega_het = std::min(parameters.omega_het, betabin_overdisp);
    parameters.omega_het_indel = std::min(parameters.omega_het_indel, betabin_overdisp);

    // Get the name of the file, without directory
    std::string input_name = input_file;
    int name_start = 0;
    for (int i = 0; i < input_file.size(); i++) {
        if (input_file[i] == '/') name_start = i + 1;
    }

    input_name = input_file.substr(name_start, input_file.size() - name_start);

    std::vector<double> results{};
    results.resize(n_chains);
    std::vector<Tree> best_trees{};
    best_trees.resize(n_chains);
    if (n_chains < omp_get_num_procs()) omp_set_num_threads(n_chains);
    else omp_set_num_threads(omp_get_num_procs());

    std::cout << "Starting to test\n";

    std::srand(0);
    savedScoresPerCellAndLoci.resize(n_loci); //temp debug

    //the purpose of this test is to read a tree and write it again, and see if the files differ
    //Tree tree_test("C:/Code/Projects/RichterMissionBio/RichterMissionBio/data/WGS_data/DFCI-5573/DFCI-5573-Merged/samtools_mpileup/compass/output_before_extra_run_1/tree.gv", true);
    //tree_test.to_dot("C:/Code/Projects/RichterMissionBio/RichterMissionBio/data/WGS_data/DFCI-5573/DFCI-5573-Merged/samtools_mpileup/compass/output_before_extra_run_1/tree_rewritten.gv", false);

    //Tree tree("C:/Code/Projects/RichterMissionBio/RichterMissionBio/data/WGS_data/DFCI-5573/DFCI-5573-RT-01/compass/DFCI-5573-RT-01_tree.gv", true);
    //Tree tree = build_debug_tree_1();
    //Tree tree = build_debug_tree_X1();
    Tree tree("./output_before_extra_run_8/tree.gv", true);

    //Tree tree = build_debug_tree_1_chr9_LOH();
    //Tree tree = build_debug_tree_CN_1();

    tree.compute_attachment_scores(true, true);

    tree.compute_likelihood(true); //no idea about the param
    tree.compute_prior_score();
    tree.update_full_score();
    auto old_score = tree.log_score;
    std::cout << "bef score: " << tree.log_score << "\n";
    //now add a CNV event that we know are there
    //Tree tree2 = build_debug_tree_1_with_CNV_node();
    //Tree tree2 = build_debug_tree_1_with_CNLOH_node_that_removes_mut();
    //Tree tree2 = build_debug_tree_1_with_CNLOH_after_LOH();
    //bool res = tree2.rec_check_max_one_event_per_region_per_lineage();
    //Tree tree2 = build_debug_tree_with_bial();//build_debug_tree_1_chr9_double_loss();
    //Tree tree2 = build_debug_tree_with_prior_init_CNV();
    //Tree tree2 = build_debug_tree_X2();
    Tree tree2("./output_before_extra_run_8/tree.gv", true);
    //tree2.nodes[1]->remove_CNA(std::make_tuple(5, -1, std::vector<int>()));
    //auto CNA_to_move = std::make_tuple(28, 1, std::vector<int>{0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1});
    //tree2.nodes[6]->remove_CNA(CNA_to_move);
    //tree2.nodes[1]->add_CNA(CNA_to_move, false);
    ////Now move the nodes around - node 4 should be under 2
    //tree2.parents[4] = 2;
    ////node 5 should be under 4 (branched)
    //tree2.parents[5] = 4;
    ////node 3 should be under 4 and have a loss of the event above
    //tree2.parents[3] = 4;
    //auto CNA_loss = std::make_tuple(28, -1, std::vector<int>{0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1});
    ////also move all of the old events to node 4

    //tree2.nodes[3]->add_CNA(CNA_loss, false);
    //tree2.compute_children();
    //auto xxs = std::vector<int>{ 1 };
    
    //move CNA 13:
    auto CNA = tree2.nodes[3]->get_CNA_from_region(13);
    tree2.nodes[3]->remove_CNA(CNA);
    tree2.nodes[4]->add_CNA(CNA, false);
    //move all X:
    CNA = tree2.nodes[2]->get_CNA_from_region(44);
    tree2.nodes[2]->remove_CNA(CNA);
    tree2.nodes[4]->add_CNA(CNA, false);
    //CNA = tree2.nodes[2]->get_CNA_from_region(45);
    //tree2.nodes[2]->remove_CNA(CNA);
    //tree2.nodes[4]->add_CNA(CNA, false);
    //CNA = tree2.nodes[2]->get_CNA_from_region(46);
    //tree2.nodes[2]->remove_CNA(CNA);
    //tree2.nodes[4]->add_CNA(CNA, false);
    //CNA = tree2.nodes[2]->get_CNA_from_region(47);
    //tree2.nodes[2]->remove_CNA(CNA);
    //tree2.nodes[4]->add_CNA(CNA, false);
    //also remove from node 5:
    //CNA = tree2.nodes[5]->get_CNA_from_region(47);
    //tree2.nodes[5]->remove_CNA(CNA);

    //Add the third parallel node as a child to node 3
    //There should be a loss28 with the opposite allele compared to node 4 and a loss34 with opposite allele compared to Node 5
    tree2.add_node_no_randomization(3);
    std::size_t ind = tree2.nodes.size() - 1; //index of the new node
    tree2.nodes[ind]->add_CNA(flip_CNA_alleles(tree2.nodes[4]->get_CNA_from_region(28)), false);
    tree2.nodes[ind]->add_CNA(flip_CNA_alleles(tree2.nodes[5]->get_CNA_from_region(34)), false);
    
    //test to remove the small node 8
    tree2.delete_node(8);

    //test to Add LOH or CNLOH 12 (whatever fits best) to node 3
    //tree2.nodes[3]->add_CNA(std::make_tuple(9, -1, std::vector<int>{1}), false);
    tree2.nodes[3]->add_CNA(std::make_tuple(9, 0, std::vector<int>{1}), false);
    

    //tree2.nodes[3]->add_CNA(std::make_tuple(9, -1, std::vector<int>{1}), false);

    
    //Tree tree2(tree);
    //auto pNode = tree2.get_node(1);
    //Node nodeOld(*pNode);
    //pNode->add_CNA(std::tuple<int, int, std::vector<int>>(0, -1, std::vector<int>{0}), false);//chr1 loss (region 1), SNPs 4,5,6,7,8, with ref/alt 0,1,1,1,0 
    //pNode->add_CNA(std::tuple<int, int, std::vector<int>>(1, -1, std::vector<int>{0, 1, 1, 1, 0}), false);//chr1 loss (region 1), SNPs 4,5,6,7,8, with ref/alt 0,1,1,1,0 
    tree2.compute_likelihood(true);
    tree2.compute_prior_score();
    tree2.update_full_score();
    auto new_score = tree2.log_score;
    auto x = new_score - old_score;
    std::cout << "score diff: " << x << "bef score : " << old_score << " after score : " << new_score << "\n";
    //cell 13 is good

    //extract likelihood for all cells, no cnv counts involved
    double tree1CellsLl = 0, tree2CellsLl = 0;
    for (std::size_t i = 0; i < tree.cells_loglik.size(); ++i) {
        tree1CellsLl += tree.cells_loglik[i];
        tree2CellsLl += tree2.cells_loglik[i];
    }
//    auto diffAttScoresSNV = tree2.nodes[2]->attachment_scores_SNV;
//    for (std::size_t i = 0; i < tree2.nodes[2]->attachment_scores_SNV.size(); ++i) {
//        diffAttScoresSNV[i] -= tree2.nodes[0]->attachment_scores_SNV[i];
//    }
    
    //look at the cell likelihoods for the two cases in node 2
//    std::vector<double> snvDiff(tree2.nodes[2]->attachment_scores_SNV.size(), 0);
//    for (std::size_t i = 0; i < tree2.nodes[2]->attachment_scores_SNV.size(); ++i) {
//        snvDiff[i] = tree2.nodes[2]->attachment_scores_SNV[i] - tree.nodes[1]->attachment_scores_SNV[i];
//    }

    //tree2.to_dot("./output/", false);
    //CNV from the counts is close to nothing - the diff between node and base are like two orders of magnitude lower than that from the variants

    //prior:
    auto priorDiff = tree2.log_prior_score - tree.log_prior_score;

    tree2.to_dot("./output/", false);

    return 0;
}
