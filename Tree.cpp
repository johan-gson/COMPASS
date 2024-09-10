#include <random>
#include <cmath>
#include <cfloat>
#include <vector>
#include <stack>
#include <iostream>
#include <fstream>

#include "Tree.h"
#include "Node.h"
#include "Scores.h"
#include "input.h"

#include <assert.h>


//global variables
extern int n_cells;
extern int n_loci;
extern int n_regions;
extern std::vector<Cell> cells;
extern Data data;
extern Params parameters;

Tree::Tree(Scores* cache, bool use_CNA, bool noRandomization):
    hastings_ratio(-1.0),
    use_CNA(use_CNA)
{
    cache_scores = cache;
    if (noRandomization) {
        n_nodes = 1; //just parent
    }
    else {
        n_nodes = 3 + (std::rand() % 8);
    }
    
    dropout_rates = std::vector<double>(n_loci,0.05);
    dropout_rates_ref = std::vector<double>(n_loci,0.05);
    dropout_rates_alt = std::vector<double>(n_loci,0.05);

    //if we have normal cells and germline variants, calculate the dropout rates to a better start value.


    //Generate a random tree
    parents.resize(n_nodes);
    parents[0]=-1; //root
    if (!noRandomization) {
        for (int i = 1; i < n_nodes; i++) {
            parents[i] = std::rand() % i; // new node attached randomly to one of the existing nodes.
        }
    }
    children.resize(n_nodes);
    compute_children();

    //Initialize the nodes.
    nodes.resize(n_nodes);
    node_probabilities.resize(n_nodes);
    for (int i=0;i<n_nodes;i++){
        nodes[i] = new Node(cache_scores);
        node_probabilities[i] = 1.0/n_nodes;
    }

    //randomly assign each somatic mutation to a node and compute the genotypes
    //if unknown, place it randomly
    //if germline, place it in root
    //if somatic, place it randomly except in root
    if (!noRandomization) {
        for (int i = 0; i < n_loci; i++) {
            switch (data.locus_to_variant_type[i]) {
            case VariantType::VT_UNKNOWN:
                nodes[std::rand() % n_nodes]->add_mutation(i);
                break;
            case VariantType::VT_GERMLINE:
                nodes[0]->add_mutation(i);
                break;
            case VariantType::VT_SOMATIC:
                nodes[1 + (std::rand() % (n_nodes - 1))]->add_mutation(i); //randomize but skip the first
                break;
            }
        }
    }
    compute_nodes_genotypes();

    // Initialize utils
    cells_attach_loglik.resize(n_cells);
    cells_loglik.resize(n_cells);
    cells_attach_prob.resize(n_cells);
    best_attachments.resize(n_cells);

    compute_likelihood();
    compute_prior_score();
    update_full_score();
}

Tree::Tree(){
    //nodes.clear();
    //doublets.clear();
    cells_attach_loglik.resize(n_cells);
    cells_loglik.resize(n_cells);
    cells_attach_prob.resize(n_cells);
    best_attachments.resize(n_cells);
}


Tree::Tree(const Tree& source):  
    // copy constructor
    use_CNA(source.use_CNA),
    n_nodes(source.n_nodes),
    parents(source.parents),
    children(source.children),
    DFT_order(source.DFT_order),
    node_probabilities(source.node_probabilities),
    dropout_rates(source.dropout_rates),
    dropout_rates_ref(source.dropout_rates_ref),
    dropout_rates_alt(source.dropout_rates_alt),
    region_probabilities(source.region_probabilities),
    candidate_regions(source.candidate_regions),
    regions_successor(source.regions_successor),
    log_prior_score(source.log_prior_score),
    log_likelihood(source.log_likelihood),
    log_score(source.log_score),
    cache_scores(source.cache_scores)
{
    // Deep copy of the nodes
    for (Node* node: source.nodes){
        nodes.push_back(new Node(*node));
    }

    // Initialize utils
    //cells_attach_loglik.resize(n_cells);
    //cells_loglik.resize(n_cells);
    //cells_attach_prob.resize(n_cells);
    cells_attach_loglik = source.cells_attach_loglik;
    cells_loglik = source.cells_loglik;
    cells_attach_prob = source.cells_attach_prob;
    best_attachments = source.best_attachments;
}



Tree::~Tree(){
    for (Node* node: nodes){
        delete node;
    }
    for (Node* doublet: doublets){
        delete doublet;
    }
}

Tree& Tree::operator=(const Tree& source){
    // delete the old nodes
    for (Node* n: nodes){
        delete n;
    }
    nodes.clear();
    if (parameters.use_doublets){
        for (Node* doublet: doublets){
            delete doublet;
        }
        doublets.clear();
    }
    use_CNA = source.use_CNA;
    cache_scores = source.cache_scores;
    n_nodes = source.n_nodes;
    parents = source.parents;
    children = source.children;
    DFT_order=source.DFT_order;
    node_probabilities = source.node_probabilities;
    
    dropout_rates = source.dropout_rates;
    dropout_rates_ref = source.dropout_rates_ref;
    dropout_rates_alt = source.dropout_rates_alt;
   
    region_probabilities = source.region_probabilities;
    candidate_regions = source.candidate_regions;
    regions_successor = source.regions_successor;
    
    log_prior_score = source.log_prior_score;
    log_likelihood = source.log_likelihood;
    log_score = source.log_score;
    best_attachments = source.best_attachments;
    cells_attach_loglik = source.cells_attach_loglik;
    cells_loglik = source.cells_loglik;
    cells_attach_prob = source.cells_attach_prob;
    
    for (int i=0;i<n_nodes;i++){
        nodes.push_back(new Node(*source.nodes[i]));
    }
    return *this;
}

void Tree::compute_children(){
    // Compute the list of children of each node, from the parent vector.
    children.resize(n_nodes);
    for (int i=0; i <n_nodes; i++){
        children[i].clear();
    }
    for (int i=1; i <n_nodes; i++){ //start from 1 because the root has no parent  
        children[parents[i]].push_back(i);
    }

    // Also compute the DFT order
    std::stack<int> stk;
    DFT_order.clear();
    stk.push(0);
    while (!stk.empty()) {
        int top = stk.top();
        DFT_order.push_back(top);
        stk.pop();
        for (int child: children[top]) {
            stk.push(child);
        };
    }
}

bool Tree::is_ancestor(int potential_ancestor, int potential_descendant){
    if (potential_ancestor==0) return true;
    int ancestor = potential_descendant;
    while (ancestor!=0){
        if (ancestor==potential_ancestor) return true;
        ancestor = parents[ancestor];
    }
    return false;
}


void Tree::compute_nodes_genotypes(){
    // Perform a depth-first traversal and compute the genotype of each node, based on the genotype of its parent.
    std::stack<int> stk;
    stk.push(0); 
    while (!stk.empty()) {
        int top = stk.top();
        stk.pop();
        for (int child: children[top]) {
            stk.push(child);
        };
        if (top!=0) nodes[top]->update_genotype(nodes[parents[top]], top == 0);
        else nodes[top]->update_genotype(nullptr, top == 0);
    }
    if (parameters.use_doublets){
        for (Node* doublet: doublets){
            delete doublet;
        }
        doublets.clear();
        for (int k=0;k<n_nodes;k++){ //k: first node in the doublet
            for (int l=k;l<n_nodes;l++){ //l: second node in the doublet
                // the first node in the doublet should be the one which is earlier is the DFT order
                // when we compute the scores of the doublet from its parent, we only update the updated loci of the second node
                int i=0;
                while (DFT_order[i]!=k && DFT_order[i]!=l) i++;
                if (DFT_order[i]==k) doublets.push_back(new Node(*nodes[k],*nodes[l]));
                else doublets.push_back(new Node(*nodes[l],*nodes[k]));
                
            }
        }
    }
}




void Tree::compute_attachment_scores(bool use_doublets_local, bool recompute_CNA_scores){
    // Compute the attachment scores of the nodes below node_index (including node_index) by performing a DFT
    std::stack<int> stk;
    stk.push(0);
    while (!stk.empty()) {
        int top = stk.top();
        stk.pop();
        for (int child: children[top]) {
            stk.push(child);
        };
        //std::cout << "Node: " << top;
        if (top == 0) // root: compute score from scratch
            nodes[top]->compute_attachment_scores(use_CNA,dropout_rates_ref,dropout_rates_alt,region_probabilities);
        else // start from the parent score, and only compute the difference on loci/regions affected by events
            nodes[top]->compute_attachment_scores_parent(use_CNA,nodes[parents[top]], dropout_rates_ref,dropout_rates_alt,region_probabilities,recompute_CNA_scores);
    }

    if (use_doublets_local){
        for (int k=0;k<n_nodes;k++){
            for (int l=k;l<n_nodes;l++){ 
                // traverse in DFT order, so that the scores of the parent are always already computed
                int n1 = DFT_order[k]; // first node of the doublet
                int n2 = DFT_order[l]; // second node of the doublet
                int idx; // index in the doublets vector. The vector only contains pairs (k,l) with k<=l
                if (n2>n1) idx = n_nodes * n1 - (n1*(n1-1))/2 + n2-n1;
                else idx = n_nodes * n2 - (n2*(n2-1))/2 + n1-n2;
                // compute doublet (0,0) from scratch
                if (n1==0 && n2==0) doublets[idx]->compute_attachment_scores(use_CNA,dropout_rates_ref,dropout_rates_alt,region_probabilities);
                // compute doublet (n1,n1) from (parent(n1),n1) (or (n1,parent(n1) )
                else if (n1==n2){
                    int idx_parent;
                    if (parents[n1]>n1) idx_parent = n_nodes * n1 - (n1*(n1-1))/2 + parents[n1]-n1;
                    else idx_parent = n_nodes * parents[n1] - (parents[n1]*(parents[n1]-1))/2 + n1- parents[n1];
                    doublets[idx]->compute_attachment_scores_parent(use_CNA,doublets[idx_parent],dropout_rates_ref,dropout_rates_alt,region_probabilities,recompute_CNA_scores);
                } 
                // compute doublet (n1,n2) from (n1,parent(n2)) (or (parent(n2),n1) )
                else {
                    int idx_parent;
                    if (parents[n2]>n1) idx_parent = n_nodes * n1 - (n1*(n1-1))/2 + parents[n2]-n1;
                    else idx_parent = n_nodes * parents[n2] - (parents[n2]*(parents[n2]-1))/2 + n1- parents[n2];
                    doublets[idx]->compute_attachment_scores_parent(use_CNA,doublets[idx_parent],dropout_rates_ref,dropout_rates_alt,region_probabilities,recompute_CNA_scores);
                }
            }
        }
    }

}



void Tree::compute_likelihood(bool allow_diff_dropoutrates){
    // Compute the likelihood by marginalizing over the attachment points

    // Remove empty nodes
    int n=1;
    while (n<n_nodes){  
        if (nodes[n]->is_empty()){
            delete_node(n);
            n=1;
        }
        else{
            n+=1;
        }
    }

    compute_nodes_genotypes();
    cache_scores->clear_cache_if_too_large();

    // EM algorithm to find best node probabilities and dropout rates
    // start with uniform node probabilities and the previous dropout rates
    for (int k=0;k<n_nodes;k++) node_probabilities[k]=1.0/n_nodes;

    avg_diff_nodeprob=10.0; // how much the node probabilities changed between 2 EM steps
    avg_diff_dropoutrates=10.0; // how much the dropout rates changed between 2 EM steps

    bool use_doublets_EM = false;
    bool recompute_CNA_scores=true; // Compute CNA scores only once, because they do not depend on the dropout rates.
    int n_loops=0;
    while((avg_diff_nodeprob>0.0005|| avg_diff_dropoutrates>0.0001) && n_loops<100){
        if (avg_diff_dropoutrates>0.0005) compute_attachment_scores(use_doublets_EM,recompute_CNA_scores); // attachment scores of cells to nodes do not depend on node probabilities
        compute_cells_likelihoods(use_doublets_EM);
        EM_step(use_doublets_EM,false);
        recompute_CNA_scores=false;
        n_loops++;
    }
    // See if likelihood can be improved by allowing, for some loci, the 2 alleles to have different dropout rates
    if (allow_diff_dropoutrates){
        EM_step(use_doublets_EM,true);
        n_loops=0;
        while((avg_diff_nodeprob>0.0005|| avg_diff_dropoutrates>0.0001) && n_loops<40){
            if (avg_diff_dropoutrates>0.0005) compute_attachment_scores(use_doublets_EM,recompute_CNA_scores); // attachment scores of cells to nodes do not depend on node probabilities
            compute_cells_likelihoods(use_doublets_EM);
            EM_step(use_doublets_EM,true);
            n_loops++;
        }
    }
    if (parameters.use_doublets && !use_doublets_EM){
        // if we did not use doublets for the EM algorithm, we need to recompute the scores with the doublets
        compute_attachment_scores(parameters.use_doublets,true); 
        compute_cells_likelihoods(parameters.use_doublets);
    }
    log_likelihood=0;
    for (int j=0;j<n_cells;j++){
        log_likelihood+=cells_loglik[j];
    }
}


void Tree::compute_cells_likelihoods(bool use_doublets_local){
    int n_attachment_points = n_nodes;
    if (use_doublets_local) n_attachment_points = (n_nodes*(n_nodes+3)) / 2; //can attach to a node or to a doublet

    for (int j=0; j<n_cells; j++){ 
        double max_loglik=-DBL_MAX;
        int best_attach=-1;
        cells_attach_loglik[j].resize(n_attachment_points);
        for (int k=0; k<n_nodes;k++){
            cells_attach_loglik[j][k] = nodes[k]->attachment_scores[j] + std::log(node_probabilities[k]);
            if (use_doublets_local) cells_attach_loglik[j][k]+=std::log(1-parameters.doublet_rate);
            if (cells_attach_loglik[j][k] > max_loglik){
                max_loglik=cells_attach_loglik[j][k];
                best_attach=k;
            }
        }
        if (use_doublets_local){
            int idx=0;
            for (int k=0;k<n_nodes;k++){
                for (int l=k;l<n_nodes;l++){
                    cells_attach_loglik[j][n_nodes+idx] = doublets[idx]->attachment_scores[j] + std::log(node_probabilities[k]) 
                    + std::log(node_probabilities[l]) +std::log(parameters.doublet_rate);
                    if (k!=l)  cells_attach_loglik[j][n_nodes+idx] += std::log(2); // the doublet (l,k) has the same probability as (k,l)
                    if ( cells_attach_loglik[j][n_nodes+idx] > max_loglik){
                        max_loglik =  cells_attach_loglik[j][n_nodes+idx];
                        best_attach = n_nodes+idx;
                    }
                    idx++;
                }
            }
        }
        cells_loglik[j] = Scores::log_sum_exp(cells_attach_loglik[j]);
        best_attachments[j] = best_attach;
    }
}

void Tree::EM_step(bool use_doublets_local, bool allow_diff_dropoutrates){
    int n_attachment_points = n_nodes;
    if (use_doublets_local) n_attachment_points = (n_nodes*(n_nodes+3)) / 2; //can attach to a node or to a doublet
    nodes_attachment_counts.resize(n_nodes);
    for (int k=0;k<n_nodes;k++) nodes_attachment_counts[k]=0.0;
    //E-step: compute probabilities that cell j is attached to node k, number of cells attached to each node and number of dropouts
    for (int j=0; j<n_cells; j++){ 
        cells_attach_prob[j].resize(n_attachment_points);

        for (int k=0;k<n_nodes;k++){
            cells_attach_prob[j][k] = std::exp(cells_attach_loglik[j][k]-cells_loglik[j]);
            nodes_attachment_counts[k]+=cells_attach_prob[j][k]; // add posterior probability that cell j is attached to node k
        }

        if (use_doublets_local){
            int idx=0;
            for (int k=0;k<n_nodes;k++){
                for (int l=k;l<n_nodes;l++){
                    cells_attach_prob[j][n_nodes+idx] = std::exp(cells_attach_loglik[j][n_nodes+idx]-cells_loglik[j]);
                    nodes_attachment_counts[k] += cells_attach_prob[j][n_nodes+idx]/2.0;
                    nodes_attachment_counts[l] += cells_attach_prob[j][n_nodes+idx]/2.0;
                    idx++;
                }
            }
        }
    }
    // count dropouts 
    double doublet_rate_local = parameters.doublet_rate;
    if (!use_doublets_local) doublet_rate_local=0.0;
    std::vector<double> dropoutref_counts(n_loci,0.0);
    std::vector<double> dropoutalt_counts(n_loci,0.0);
    std::vector<double> nrefalleles(n_loci,0.0);
    std::vector<double> naltalleles(n_loci,0.0);

    for (int i = 0; i < n_loci; i++) {
        std::map<std::pair<int, int>, std::vector<double>> genotypes_prob; // for each cell, probability that it is attached to a node with a given genotype for locus i
        // First, group nodes depending on their genotype at locus i
        for (int k = 0; k < n_nodes; k++) {
            int n_ref = nodes[k]->get_n_ref_allele(i);
            int n_alt = nodes[k]->get_n_alt_allele(i);
            std::pair<int, int> p = std::make_pair(n_ref, n_alt); //p: genotype at locus i for node k
            if (!genotypes_prob.count(p)) {
                genotypes_prob[p] = std::vector<double>(n_cells, 0.0);
            }
            for (int j = 0; j < n_cells; j++) {
                genotypes_prob[p][j] += cells_attach_prob[j][k] / (1.0 - doublet_rate_local); // do not use doublets for inferring dropout rates (faster)
            }
        }

        // Count dropouts for each genotype at locus i
        for (const auto& m : genotypes_prob) {
            std::pair<int, int> p = m.first;
            int n_ref = p.first;
            int n_alt = p.second;
            if (n_ref > 0 && n_alt > 0) { // only consider dropouts for heterozygous genotypes. 
                const std::vector<double>& dropoutsref = cache_scores->get_dropoutref_counts_genotype(n_ref, n_alt, i, dropout_rates_ref[i],
                    dropout_rates_alt[i]);
                const std::vector<double>& dropoutsalt = cache_scores->get_dropoutalt_counts_genotype(n_ref, n_alt, i, dropout_rates_ref[i],
                    dropout_rates_alt[i]);
                for (int j = 0; j < n_cells; j++) {
                    // if we have no reads at this position in this cell, we can't tell if a dropout occurred.
                    if (cells[j].ref_counts[i] + cells[j].alt_counts[i] > 4) {
                        dropoutref_counts[i] += genotypes_prob[p][j] * dropoutsref[j];
                        dropoutalt_counts[i] += genotypes_prob[p][j] * dropoutsalt[j];
                        nrefalleles[i] += genotypes_prob[p][j] * n_ref;
                        naltalleles[i] += genotypes_prob[p][j] * n_alt;
                    }
                }
            }
        }
    }


    //M-step: optimize the node probabilities and dropout rates
    avg_diff_nodeprob = 0.0;
    for (int k = 0; k < n_nodes; k++) { // update node probabilities
        double prev = node_probabilities[k];
        node_probabilities[k] = nodes_attachment_counts[k] / n_cells;
        avg_diff_nodeprob += std::abs(node_probabilities[k] - prev) / n_nodes;
    }

    avg_diff_dropoutrates = 0.0;
    for (int i = 0; i < n_loci; i++) { // update dropout rates
        double prev_ref = dropout_rates_ref[i];
        double prev_alt = dropout_rates_alt[i];
        dropout_rates[i] = ((parameters.prior_dropoutrate_mean * parameters.prior_dropoutrate_omega - 1) * 2 + dropoutref_counts[i] + dropoutalt_counts[i])
            / ((parameters.prior_dropoutrate_omega - 2) * 2 + nrefalleles[i] + naltalleles[i]);
        if (allow_diff_dropoutrates) {
            // See if likelihood can be improved by allowing 2 different dropout rates
            dropout_rates_ref[i] = (parameters.prior_dropoutrate_mean * parameters.prior_dropoutrate_omega - 1 + dropoutref_counts[i])
                / (parameters.prior_dropoutrate_omega - 2 + nrefalleles[i]);
            dropout_rates_alt[i] = (parameters.prior_dropoutrate_mean * parameters.prior_dropoutrate_omega - 1 + dropoutalt_counts[i])
                / (parameters.prior_dropoutrate_omega - 2 + naltalleles[i]);
            double diff_dropoutrates_lik = dropoutref_counts[i] * std::log(dropout_rates_ref[i]) + (nrefalleles[i] - dropoutref_counts[i]) * std::log(1 - dropout_rates_ref[i])
                + dropoutalt_counts[i] * std::log(dropout_rates_alt[i]) + (naltalleles[i] - dropoutalt_counts[i]) * std::log(1 - dropout_rates_alt[i])
                + (parameters.prior_dropoutrate_mean * parameters.prior_dropoutrate_omega - 1) * (std::log(dropout_rates_ref[i]) + std::log(dropout_rates_alt[i]))
                + ((1.0 - parameters.prior_dropoutrate_mean) * parameters.prior_dropoutrate_omega - 1) * (std::log(1.0 - dropout_rates_ref[i]) + std::log(1.0 - dropout_rates_alt[i]));
            double same_dropoutrate_lik = (dropoutref_counts[i] + dropoutalt_counts[i]) * std::log(dropout_rates[i])
                + (nrefalleles[i] + naltalleles[i] - dropoutref_counts[i] - dropoutalt_counts[i]) * std::log(1 - dropout_rates[i])
                + (parameters.prior_dropoutrate_mean * parameters.prior_dropoutrate_omega - 1) * 2 * std::log(dropout_rates[i])
                + ((1.0 - parameters.prior_dropoutrate_mean) * parameters.prior_dropoutrate_omega - 1) * 2 * std::log(1.0 - dropout_rates[i]);
            if (same_dropoutrate_lik > diff_dropoutrates_lik - 60) {
                dropout_rates_ref[i] = dropout_rates[i];
                dropout_rates_alt[i] = dropout_rates[i];
            }
        }
        else {
            dropout_rates_ref[i] = dropout_rates[i];
            dropout_rates_alt[i] = dropout_rates[i];
        }
        // Set minimum and maximum dropout rate
        if (dropout_rates[i] < 0.01) dropout_rates[i] = 0.01;
        if (dropout_rates_ref[i] < 0.01) dropout_rates_ref[i] = 0.01;
        if (dropout_rates_alt[i] < 0.01) dropout_rates_alt[i] = 0.01;

        if (dropout_rates[i] > 0.50) dropout_rates[i] = 0.50;
        if (dropout_rates_ref[i] > 0.50) dropout_rates_ref[i] = 0.50;
        if (dropout_rates_alt[i] > 0.50) dropout_rates_alt[i] = 0.50;

        avg_diff_dropoutrates += (std::abs(dropout_rates_ref[i] - prev_ref) + std::abs(dropout_rates_alt[i] - prev_alt)) / n_loci;

    }

}

double Tree::rec_check_max_one_event_per_region_per_lineage(int node, std::vector<std::vector<int>> previousEvents) {
    // check that each region is affected, in one lineage, by at most one CNA.
    // we do however allow for two cases: 
    // 1) Loss above in the tree followed by a gain, leading to a CNLOH. This is a common scenario in biology.
    // 2) LOH above followed by a second LOH, resulting in a double loss. This is only allowed 
    //    if parameters.allow_double_allele_loss is true
    double penalty = 0;
    if (node == 0) {
        previousEvents.resize(data.region_to_name.size());
    }

    for (auto CNA : nodes[node]->get_CNA_events()) {
        previousEvents[std::get<0>(CNA)].push_back(std::get<1>(CNA));
        if (previousEvents[std::get<0>(CNA)].size() > 2) {
            return 100000;
        }
        else if (previousEvents[std::get<0>(CNA)].size() == 2) {
            if (previousEvents[std::get<0>(CNA)][0] != -1) { //must start with a loss
                return 100000;
            }
            if ((previousEvents[std::get<0>(CNA)][1] == 0) || //CHLOH after loss, not ok
                ((previousEvents[std::get<0>(CNA)][1] == -1) && !parameters.allow_double_allele_loss)) { //Loss + loss only ok with param set, loss + gain always ok, will automatically be CNLOH
                return 10000;
            }
            //we still penalize these cases a bit
            penalty += parameters.two_CNA_in_lineage_penalty;
        }
    }

    for (int child : children[node]) {
        penalty += rec_check_max_one_event_per_region_per_lineage(child, previousEvents);
    }
    return penalty;
}

int Tree::rec_get_number_of_cnloh_removing_somatic_mut(int node, std::vector<int> muts) const {
    int CNA_events_removing_muts = 0;
    
    //add the new non_germline mutations to muts unless we are in root
    if (node != 0) {
        //since all germline variants should be in base, we don't check the type - there will not be any in these vectors
        muts.insert(muts.end(), nodes[node]->get_mutations().begin(), nodes[node]->get_mutations().end());
    }
    
    // check if we have any CNLOH that removes the muts in this node
    for (auto CNA : nodes[node]->get_CNA_events()) {
        if (std::get<1>(CNA) == 0) { // check if copy neutral event (CNLOH)
            const auto& loci = data.region_to_loci[std::get<0>(CNA)];
            std::vector<int> lociAltLost;
            for (std::size_t i = 0; i < loci.size(); ++i) {
                if (data.locus_to_variant_type[loci[i]] != VariantType::VT_GERMLINE) { //this check could be unneccesary, since they exist in root only
                    if (std::get<2>(CNA)[i] == 1) { //check if alt allele lost
                        lociAltLost.push_back(loci[i]);
                    }
                }
            }
            //check for intersection between muts and loci - if any exist, we have a hit
            if (std::find_first_of(lociAltLost.begin(), lociAltLost.end(),
                muts.begin(), muts.end()) != lociAltLost.end()) {
                ++CNA_events_removing_muts;
            }
        }
    }

    //sum up the bad CNLOH from the children
    for (int child : children[node]) {
        CNA_events_removing_muts += rec_get_number_of_cnloh_removing_somatic_mut(child, muts);
    }
    return CNA_events_removing_muts;
}

std::set<int> Tree::rec_get_nodes_in_lineages_with_2_CNA(int region, int node, std::set<int> parents, std::vector<int> previousEvents) const {
    // find all nodes participating in a lineage with at least 2 CNA for a certain region

    if (node == 0) {
        previousEvents.resize(data.region_to_name.size());
    }
    std::set<int> results;

    for (auto CNA : nodes[node]->get_CNA_events()) {
        if (std::get<0>(CNA) == region) {
            previousEvents.push_back(std::get<1>(CNA));
            if (previousEvents.size() >= 2) {
                //we know all children + this node + the parents are in a lineage with 2, so just sum them up
                results.insert(parents.begin(), parents.end());
                results.insert(node);
                auto desc = rec_get_descendents(node);
                results.insert(desc.begin(), desc.end());
                return results;
            }
        }
    }

    std::set<int> new_parents(parents);
    new_parents.insert(node);
    for (int child : children[node]) {
        auto desc = rec_get_nodes_in_lineages_with_2_CNA(region, child, new_parents, previousEvents);
        results.insert(desc.begin(), desc.end());
    }
    //if we found anything down the tree, we need to add this node and the parents since they are also of that lineage
    if (!results.empty()) {
        results.insert(parents.begin(), parents.end());
        results.insert(node);
    }
    return results;
}

std::set<int> Tree::rec_get_descendents(int node) const {
    std::set<int> results;
    for (int child : children[node]) {
        results.insert(child);
        auto desc = rec_get_descendents(child);
        results.insert(desc.begin(), desc.end());
    }
    return results;
}

std::set<std::tuple<int, int>> Tree::rec_get_all_cna_events_in_tree(int node) const {
    std::set<std::tuple<int, int>> results;
    //first add CNAs from this node
    for (const auto& cna : nodes[node]->get_CNA_events()) {
        results.insert(std::make_tuple(std::get<0>(cna), std::get<1>(cna)));
    }
    //then recursively that from all children 
    for (int child : children[node]) {
        auto child_cnv = rec_get_all_cna_events_in_tree(child);
        results.insert(child_cnv.begin(), child_cnv.end());
    }
    return results;
}

double Tree::rec_get_cna_in_multiple_branches_penalty(int node) const {
    double penalty = 0.0;
    //we only look for the same cna in different children, so the CNAs of this node are not involved
    //we only compare region and type, not the SNPs for the CNV
    std::set<std::tuple<int, int>> CNAs;
    for (int child : children[node]) {
        auto child_CNAs = rec_get_all_cna_events_in_tree(child);
        for (auto CNA : child_CNAs) {
            if (CNAs.count(CNA) != 0) {
                penalty += parameters.cna_in_multiple_branches_penalty;
            }
        }
        CNAs.insert(child_CNAs.begin(), child_CNAs.end());
    }

    return penalty;
}



void Tree::compute_prior_score(){
    // Penalize number of nodes
    log_prior_score=-n_nodes*(4+n_loci)*parameters.node_cost; 
    // Forbid empty nodes and penalize nodes with only CNAs (since the order of CNAs is generally less reliable)
    for (int i=1;i<n_nodes;i++){
        if (nodes[i]->get_number_mutations()==0 && nodes[i]->get_number_CNA()==0) log_prior_score-=100000;
        if (nodes[i]->get_number_mutations()==0) log_prior_score-=(8+n_loci)*parameters.node_cost;
        if (nodes[i]->get_number_mutations()==0 && nodes[i]->get_number_effective_LOH(nodes[parents[i]])==0) log_prior_score-=2*(8+n_loci)*parameters.node_cost;
    }

    //give bonus to germline variants in the root
    for (int mut : nodes[0]->get_mutations()){
        switch (data.locus_to_variant_type[mut]) {
        case VariantType::VT_UNKNOWN:
            //do nothing, handled below
            break;
        case VariantType::VT_GERMLINE:
            //give bonus to having germline in root
            log_prior_score += parameters.germline_root_bonus;
            break;
        case VariantType::VT_SOMATIC:
            //do nothing, handled below
            break;
        }
    }

    //give bonus to mut variants in the root, and penalize unknowns not in root
    for (int k=1;k<n_nodes;k++){
        for (int mut : nodes[k]->get_mutations()){
            switch (data.locus_to_variant_type[mut]) {
            case VariantType::VT_UNKNOWN:
                log_prior_score += parameters.mut_notAtRoot_cost; //most unknown variants are at the root, so penalize diverging from that
                //penalize to put variants with high germline frequency outside the root
                if (data.locus_to_freq[mut] > 0.0001) log_prior_score -= data.locus_to_freq[mut] * parameters.mut_notAtRoot_freq_cost;
                break;
                break;
            case VariantType::VT_GERMLINE:
                //do nothing, handled above
                break;
            case VariantType::VT_SOMATIC:
                //give bonus since they are not in root
                log_prior_score += parameters.somatic_non_root_bonus;
                break;
            }
        }
    }


    // Penalize cells depending on cell type:
    // If they are normal, give a high penalty for being anywhere else except in the root - i.e., a bonus for being in the root
    // If they are malignant - give a penalty for being in the root. These are usually more uncertain in classification (some could be normal), so penalize a bit less
    
    for (int j = 0; j < n_cells; j++) {
        auto x = std::exp(cells_attach_loglik[j][0]);
        switch (cells[j].cell_type) {
        case CellType::CT_N:
            log_prior_score += cells_attach_prob[j][0] * parameters.normal_not_at_root_penalty; //unclear if this should be logged or tranformed in a different way
            break;
        case CellType::CT_M:
            log_prior_score -= cells_attach_prob[j][0] * parameters.malignant_at_root_penalty; //unclear if this should be logged or tranformed in a different way
            break;
        case CellType::CT_U:
            //do nothing if the cell type is unknown
            break;
        }
    }
    


    // Higher penalty when there are more cells
    double ncells_coef = 0.2 + 1.0*n_cells/8000.0;

    // Penalize CNA events
    log_prior_score-= ncells_coef* (parameters.LOH_cost+parameters.CNA_cost) /10.0 *  nodes[0]->get_number_CNA();  // Smaller penalty for CNLOH events at the root - this never happens anymore
    for (int n=1;n<n_nodes;n++){ 
        log_prior_score-= ncells_coef*parameters.CNA_cost * nodes[n]->get_number_disjoint_CNA(regions_successor); 
        // Higher penalty for CNAs resulting in LOH, because they have a bigger impact on the likelihood
        log_prior_score-= ncells_coef*parameters.LOH_cost * nodes[n]->get_number_disjoint_LOH(regions_successor);  
    }

    // Cannot have a CNA event at the root.
    if (nodes[0]->get_number_CNA_noncopyneutral()>0) log_prior_score-= 100000; 
    // One lineage cannot have more than one CNA affecting each region (but it is still possible to have events affecting the same region in parallel branches)
    log_prior_score -= rec_check_max_one_event_per_region_per_lineage();
    //If you have a somatic mutation, it is an unlikely event to have it removed by a CNLOH - the algorthm tends to do this. Penalize this
    log_prior_score -= rec_get_number_of_cnloh_removing_somatic_mut()* parameters.cnloh_removing_somatic_mut_penalty;
    //An undesired behavior is that CNAs tend to end up in two children instead of in the root, presumably since it is possible to adapt better to noise
    //We therefore penalize to have the same event in two child branches of the same node. It is still possible, but requires more evidence now.
    log_prior_score -= rec_get_cna_in_multiple_branches_penalty();


    // Dropout rates
    // with probability lambda, both alleles have the same dropout rate. With probability 1-lambda, they have the same dropout rate
    for (int i=0;i<n_loci;i++){
        if (std::abs(dropout_rates_ref[i]-dropout_rates_alt[i])>0.001){
            log_prior_score-=70; 
        }
        // Dropout rates are sampled from a beta distribution.
        log_prior_score+= (parameters.prior_dropoutrate_mean * parameters.prior_dropoutrate_omega-1) 
                                * (std::log(dropout_rates_ref[i]) + std::log(dropout_rates_alt[i])) 
                        + ((1.0-parameters.prior_dropoutrate_mean) * parameters.prior_dropoutrate_omega-1) 
                                * (std::log(1.0-dropout_rates_ref[i]) + std::log(1.0-dropout_rates_alt[i]));
    }
}

void Tree::update_full_score(){
    log_score = log_likelihood + log_prior_score;
}


void Tree::to_dot(std::string filename, bool simplified){
    // Save the tree structure in dot format (for visualization)

    std::vector<std::string> colors{"lightcoral","skyblue3","sandybrown","paleturquoise3","thistle","darkolivegreen3","lightpink","mediumpurple",
                    "darkseagreen3","navajowhite","gold"};

    // If filename ends with .gv: only output the tree in graphviz format. Otherwise output tree and cell assignments to nodes.
    bool full_output = true;
    if (filename.size()>3 &&  filename.substr(filename.size()-3)==".gv"){
        full_output = false;
    }
    std::string basename(filename);
    if (full_output){
        filename = filename + "_tree.gv";
    }
    else{
        basename = filename.substr(0,filename.size()-3);
    }
    std::ofstream out_file(filename);

    out_file <<"digraph G{"<<std::endl;
    out_file <<"node [color=dimgray fontsize=24 fontcolor=black fontname=Helvetica penwidth=5];"<<std::endl;
    for (int i=1;i<n_nodes;i++){
        if (parameters.verbose) std::cout<<i<< " is a child of "<<parents[i]<<std::endl;
        out_file<<parents[i]<<" -> "<<i<<" [color=dimgray penwidth=4 weight=2];"<<std::endl;
    }

    // Identify mutations at the root which are not affected by a CNA
    std::set<int> excluded_mutations{};
    if (simplified){
        if (n_nodes>0){
            for (int m: nodes[0]->get_mutations()){
                bool affected_by_event=false;
                for (int n=0;n<n_nodes;n++){
                    for (auto CNA: nodes[n]->get_CNA_events()){
                        if (data.locus_to_region[m]==std::get<0>(CNA)) affected_by_event = true;
                    }
                }
                if (!affected_by_event) excluded_mutations.insert(m);
            }
        }
    }

    for (int i=0;i<n_nodes;i++){
        if (simplified) out_file<<i<<"[label=<"<<nodes[i]->get_label_simple(excluded_mutations)<<">];"<<std::endl;
        else out_file<<i<<"[label=<"<<nodes[i]->get_label()<<">];"<<std::endl;
    }

    for (int k=0;k<n_nodes;k++){
        out_file<<k<<" -> "<<k+n_nodes<<" [dir=none style=dashed weight=1 penwidth=5 color="<<colors[k%colors.size()]<<"];"<<std::endl;
    }
    std::vector<int> count_nodes(n_nodes,0);
    int total=0;
    for (int j=0;j<n_cells;j++){
        if (best_attachments[j]>=0 && best_attachments[j]<n_nodes){
             count_nodes[best_attachments[j]]++;
             total++;
        }
    }
    for (int k=0;k<n_nodes;k++){
        double size = std::sqrt(100.0*count_nodes[k]/total) /3.0;
        out_file<<k+n_nodes<<"[label=\""<<count_nodes[k]<<" cells\\n"<<std::round(100.0*count_nodes[k]/total)<<"\\%\""<<" style = filled width="<<size
        <<" height="<<size<<" color="<<colors[k%colors.size()]<<"];"<<std::endl;
    }

    out_file <<"}"<<std::endl;
    out_file.close();
    if (parameters.verbose) std::cout<<"Node probabilities"<<std::endl;
    for (int n=0;n<n_nodes;n++){
        if (parameters.verbose) std::cout<<n<<": "<<node_probabilities[n]<<std::endl;
    }
    if (parameters.verbose) std::cout<<"Dropout rates"<<std::endl;
    for (int i=0;i<n_loci;i++){
        if (parameters.verbose) std::cout<<i<<" ("<<data.locus_to_name[i]<<"): "<<dropout_rates[i]<<" (ref:" <<dropout_rates_ref[i]<<", alt:"<<dropout_rates_alt[i]<<")"<<std::endl;
    }

    if (full_output){
        // Recompute assignment probabilities, only taking singlets into account
        std::vector<std::vector<double>> cells_attach_loglik_singlet; 
        std::vector<double> cells_loglik_singlet;
        std::vector<int> best_attachments_singlet;
        cells_attach_loglik_singlet.resize(n_cells);
        cells_loglik_singlet.resize(n_cells);
        best_attachments_singlet.resize(n_cells);
        for (int j=0;j<n_cells;j++){
            cells_attach_loglik_singlet[j].resize(n_nodes);
            double best_attach_score=-DBL_MAX;
            for (int k=0;k<n_nodes;k++){
                cells_attach_loglik_singlet[j][k] = cells_attach_loglik[j][k];
                if (cells_attach_loglik_singlet[j][k]>best_attach_score){
                    best_attach_score = cells_attach_loglik_singlet[j][k];
                    best_attachments_singlet[j] = k;
                }
            }
            cells_loglik_singlet[j] = Scores::log_sum_exp(cells_attach_loglik_singlet[j]);
        }

        // Assignments of cells to nodes 
        std::ofstream out_file_cell_assignments(basename+"_cellAssignments.tsv");
        std::ofstream out_file_cell_assignment_probs(basename+"_cellAssignmentProbs.tsv");

        //Header
        out_file_cell_assignments <<"cell\tnode\tdoublet"<<std::endl;
        out_file_cell_assignment_probs <<"cell";
        for (int k=0; k < n_nodes;k++){
            out_file_cell_assignment_probs<<"\tNode "<<k;
        }
        /*if (parameters.use_doublets){
            std::cout<<cells_attach_loglik[0].size()<<std::endl;
                for (int k=0;k<n_nodes;k++){
                    for (int l=k;l<n_nodes;l++){
                        out_file_cell_assignment_probs << "\tDoublet "<<k<<","<<l;
                    }
                }
        }*/
        out_file_cell_assignment_probs<<std::endl;
        
        // Content
        for (int j=0;j<n_cells;j++){
            out_file_cell_assignments << cells[j].name<<"\t"<<best_attachments_singlet[j];
            if (best_attachments[j]>=n_nodes) out_file_cell_assignments<<"\tyes"<<std::endl;
            else out_file_cell_assignments<<"\tno"<<std::endl;
            
            out_file_cell_assignment_probs<<cells[j].name;
            for (int k=0; k < n_nodes;k++){
                out_file_cell_assignment_probs<<"\t"<<std::exp(cells_attach_loglik_singlet[j][k]-cells_loglik_singlet[j]);
            }
            /*if (parameters.use_doublets){
                int idx=0;
                for (int k=0;k<n_nodes;k++){
                    for (int l=k;l<n_nodes;l++){
                        out_file_cell_assignment_probs << "\t"<<std::exp(cells_attach_loglik[j][n_nodes+idx]-cells_loglik[j]);
                        idx++;
                    }
                }
            }*/
            out_file_cell_assignment_probs<<std::endl;
        }
         out_file_cell_assignments.close();
         out_file_cell_assignment_probs.close();

        // ----------------------------------
        // Tree in json format
        std::ofstream out_file_json(basename+"_tree.json");
        out_file_json<<"{"<<std::endl;
        out_file_json<<"\"nodes\":["<<std::endl;
        for (int k=0;k<n_nodes;k++){
            out_file_json<<"\t{"<<std::endl;
            out_file_json<<"\t\t\"name\": \"Node "<<k<<"\","<<std::endl;
            if (parents[k]>=0){
                out_file_json<<"\t\t\"parent\": \"Node "<<parents[k]<<"\","<<std::endl;
            }
            else{
                out_file_json<<"\t\t\"parent\": \"-\","<<std::endl;
            }
            
            //SNV
            out_file_json<<"\t\t\"SNV\": [";
            std::vector<int> SNVs = nodes[k]->get_mutations();
            if (SNVs.size()>0){
                out_file_json<<"\""<<data.locus_to_name[SNVs[0]]<<"\"";
                for (int i=1;i<SNVs.size();i++){
                    out_file_json<<",\""<<data.locus_to_name[SNVs[i]]<<"\"";
                }
            }
            out_file_json<<"],"<<std::endl;

            //CNA
            out_file_json<<"\t\t\"CNA\": [";
            std::set<std::tuple<int,int,std::vector<int>>> CNAs = nodes[k]->get_CNA_events();
            if (CNAs.size()>0){
                bool first=true;
                for (auto CNA: CNAs){
                    if (!first) {
                         out_file_json<<",";
                    }
                    out_file_json<<"\"";
                    if (std::get<1>(CNA)>0){
                         out_file_json<<"Gain ";
                    }
                    else if (std::get<1>(CNA) == -1) {
                        out_file_json << "Loss ";
                    }
                    else if (std::get<1>(CNA) == -2) {
                        out_file_json << "Bial. loss ";
                    }
                    else{
                        out_file_json<<"CNLOH ";
                    }
                     out_file_json<<data.region_to_name[std::get<0>(CNA)];
                     out_file_json<<"\"";
                     first=false;
                }
            }
            out_file_json<<"]"<<std::endl;

            out_file_json<<"\t}";
            if (k<n_nodes-1){
                out_file_json<<",";
            }
            out_file_json<<std::endl;
        }

        out_file_json<<"]"<<std::endl;
        out_file_json<<"}"<<std::endl;
        out_file_json.close();

        // ----------------------------------
        // Node genotypes 

        std::ofstream out_file_genotypes(basename+"_nodes_genotypes.tsv");

        // Header
        out_file_genotypes<<"node";
        for (std::string name : data.locus_to_name){
            out_file_genotypes << "\t"<<name;
        }
        out_file_genotypes<<std::endl;

        // Content
        for (int k=0;k<n_nodes;k++){
            out_file_genotypes<<"Node "<<k;
            for (int i=0;i<n_loci;i++){
                if (nodes[k]->get_n_alt_allele(i)==0){
                     out_file_genotypes<<"\t0";
                }
                else if (nodes[k]->get_n_ref_allele(i)>0 ){
                     out_file_genotypes<<"\t1";
                }
                else{
                    out_file_genotypes<<"\t2";
                }
            }
            out_file_genotypes<<std::endl;
        }

        out_file_genotypes.close();

        // ----------------------------------
        // Node copy numbers
        if (use_CNA){

            std::ofstream out_file_copynumbers(basename+"_nodes_copynumbers.tsv");

            // Header
            out_file_copynumbers<<"node";
            for (std::string name : data.region_to_name){
                out_file_copynumbers << "\t"<<name;
            }
            out_file_copynumbers<<std::endl;

            // Content
            for (int k=0;k<n_nodes;k++){
                out_file_copynumbers<<"Node "<<k;
                for (int i=0;i<n_regions;i++){
                    out_file_copynumbers<<"\t"<<nodes[k]->get_cn_region(i);
                }
                out_file_copynumbers<<std::endl;
            }

            out_file_copynumbers.close();
        }
    }
    
}


Tree::Tree(std::string gv_file, bool use_CNA_arg): //Create tree from a graphviz file
    hastings_ratio(-1.0)     
{
    for (int i=0;i<n_loci;i++){
        dropout_rates.push_back(0.05);
        dropout_rates_ref.push_back(0.05);
        dropout_rates_alt.push_back(0.05);
    }
    cache_scores = new Scores();
    std::ifstream file(gv_file);
    std::string line;
    //skip first 2 lines
    getline (file, line);
    getline (file, line);
    n_nodes=1;
    parents.resize(1);
    parents[0]=-1;
    // Read parents
    bool finished_reading_parents=false;
    while (!finished_reading_parents) {
        getline (file, line);
        if (line.find("->")==std::string::npos) finished_reading_parents=true;
        else{
            int idx=1;
            while (line[idx]!=' ') idx++;
            int parent = stoi(line.substr(0,idx));
            idx++;
            while (line[idx]!=' ') idx++;
            idx++;
            int idx2=idx+1;
            while (line[idx2]!=' '&& line[idx]!=';') idx2++;
            int child = stoi(line.substr(idx,idx2-idx));
            if (child+1 >n_nodes) n_nodes = child+1;
            if (parent+1 > n_nodes) n_nodes = parent+1;
            parents.resize(n_nodes);
            parents[child] = parent;
        }
    }
    // Create nodes 
    for (int i=0;i<n_nodes;i++){
        nodes.push_back(new Node(cache_scores));
    }
    node_probabilities.resize(n_nodes);
    // Read labels (events) and fill the nodes with the events
    bool finished_reading_labels=false;
    while (!finished_reading_labels){
        int idx=1;
        while (idx < line.size() && line[idx]!='[') idx++;
        if (idx>=line.size() || line[idx+1]!='l') finished_reading_labels=true; //empty line
        else{
            int node =  stoi(line.substr(0,idx));
            while (line[idx]!='<') idx++;
            idx++;
            while (line[idx]==' ') idx++;
            bool finished_reading_line=false;
            if (line[idx]=='>') finished_reading_line=true;
            while (!finished_reading_line){
                if (line[idx]=='<') idx+=3; // remove <B>
                if ((line[idx]=='C' && line[idx+2]=='L')){ // CNLOH event
                    idx+=6;
                    int idx2 = idx+1;
                    while (line[idx2]!=':') idx2++;
                    int region = stoi(line.substr(idx,idx2-idx));
                    idx2++;
                    while (line[idx2]!=':') idx2++; // skip region name
                    idx = idx2+1;
                    std::vector<int> lost_alleles{};
                    while (line[idx]!='<' && line[idx]!='b' && line[idx]!='/'){
                        if (line[idx]=='0') lost_alleles.push_back(0);
                        else lost_alleles.push_back(1);
                        idx+=2;
                    }
                    idx--;
                    nodes[node]->add_CNA(std::make_tuple(region,0,lost_alleles), node == 0);
                }
                else if (line[idx]=='L' && line[idx+1]=='o'){ // Loss
                    idx+=5;
                    int idx2 = idx+1;
                    while (line[idx2]!=':') idx2++;
                    int region =  stoi(line.substr(idx,idx2-idx));
                    idx = idx2+1;
                    while (line[idx]!=':') idx++;
                    std::vector<int>alleles{};
                    idx++;
                    while (line[idx]!='<' && line[idx]!='b' && line[idx]!='/'){
                        if (line[idx]=='0') alleles.push_back(0);
                        else alleles.push_back(1);
                        idx+=2;
                    }
                    idx--;
                    nodes[node]->add_CNA(std::make_tuple(region,-1,alleles), node == 0);
                }
                else if (line[idx]=='G' && line[idx+1]=='a'){ // Gain
                    idx+=5;
                    int idx2 = idx+1;
                    while (line[idx2]!=':') idx2++;
                    int region =  stoi(line.substr(idx,idx2-idx));
                    idx = idx2+1;
                    while (line[idx]!=':') idx++;
                    std::vector<int>alleles{};
                    idx++;
                    while (line[idx]!='<' && line[idx]!='b' && line[idx]!='/'){
                        if (line[idx]=='0') alleles.push_back(0);
                        else alleles.push_back(1);
                        idx+=2;
                    }
                    idx--;
                    nodes[node]->add_CNA(std::make_tuple(region,1,alleles), node == 0);
                }
                else{ // somatic mutation
                    int idx2= idx+1;
                    while (line[idx2] != ':') idx2++;
                    int locus = atoi(line.substr(idx, idx2 - idx).c_str());
                    nodes[node]->add_mutation(locus);
                    //std::replace(strLocName.begin(), strLocName.end(), '_', ':'); //a bit risky, if variants contain _
                    //for some reason, g++ cannot compile the line below, so I replaced it with a loop
                    //std::size_t locus = std::find(data.locus_to_name.begin(), data.locus_to_name.end(), (std::string)strLocName) - data.locus_to_name.begin();
                    //if (locus <= data.locus_to_name.size()) {
                    //    nodes[node]->add_mutation(locus);
                    //}
                    //else {
                    //    throw std::runtime_error("Couldn't read tree");
                    //}
                    //for (std::size_t i = 0; i < data.locus_to_name.size(); ++i) {
                    //    if (data.locus_to_name[i] == strLocName) {
                    //        nodes[node]->add_mutation(i);
                    //    }
                    //}
                }
                while (line[idx]!='<') idx++;

                //HTML tags: events separated by <br/>
                if (line[idx+1]=='/') idx+=4; // </B>
                idx+=5;
                if (line[idx]=='>') finished_reading_line=true; 
            }
         } 
        if (!std::getline (file, line)) finished_reading_labels=true;
     }
    // Initialize utils
    cells_attach_loglik.resize(n_cells);
    cells_loglik.resize(n_cells);
    cells_attach_prob.resize(n_cells);
    best_attachments.resize(n_cells);
    use_CNA=false;
    compute_children();
    compute_likelihood(true);
    if(use_CNA_arg && select_regions()){
        use_CNA=true;
        compute_likelihood(true);
    }
    compute_prior_score();
    update_full_score();

    // Close the file
    file.close();

}

bool Tree::select_regions(int index) {
    // Find regions which might contain a non-copy-neutral CNA event. Return true if it was possible to estimate node regions (if the node contains enough cells), false otherwise.
    candidate_regions.clear();
    region_probabilities.resize(n_regions);
    int root = 0;

    // Compute number of cells attached to each node and make sure that there is a sufficient number of cells attached to the root
    std::vector<int> nodes_nbcells(n_nodes, 0);
    for (int j = 0; j < n_cells; j++) {
        if (best_attachments[j] < n_nodes && best_attachments[j] >= 0) {
            nodes_nbcells[best_attachments[j]]++;
        }
    }

    // Either used predetermined region weights or estimate them using cells attached to the root.
    if (data.predetermined_region_weights.size() > 0) {
        for (int k = 0; k < n_regions; k++) region_probabilities[k] = data.predetermined_region_weights[k];
    }
    else {
        if (nodes_nbcells[0] < std::max(40.0, 0.015 * n_cells)) { // not enough cells attached to the root
            // If some nodes have CNLOH, try to use one node without CNLOH in its ancestors as the root.
            std::vector<bool> CNLOH_in_ancestors(n_nodes, false);
            bool CNLOH_in_tree = false;
            int node_newroot = -1;
            for (int n : DFT_order) {
                if (n > 0) {
                    if (nodes[n]->get_number_CNA() > 0 || CNLOH_in_ancestors[parents[n]]) {
                        CNLOH_in_ancestors[n] = true;
                        CNLOH_in_tree = true;
                    }
                    else {
                        if (nodes_nbcells[n] >= std::max(40.0, 0.015 * n_cells) && node_newroot == -1) {
                            node_newroot = n;
                        }
                    }
                }
            }
            if (node_newroot == -1 || (!CNLOH_in_tree)) {
                if (index >= 0) std::cout << "Chain " << std::to_string(index) << ": ";
                std::cout << "In the tree inferred without CNAs, there were not enough cells attached to the root (" << nodes_nbcells[0] << ") to estimate the region weights, so COMPASS could not attempt to find CNAs." << std::endl;
                return false; // not enough cells attached to the root: cannot find CNAs
            }
            else {
                root = node_newroot;
            }
        }
    }
    //  Compute average region probability in each node
    std::vector<std::vector<std::vector<double>>> nodes_regionprobs{};
    nodes_regionprobs.resize(n_regions);
    for (int k = 0; k < n_regions; k++) {
        if (!data.region_is_reliable[k]) continue;
        std::vector<double> nodes_averages(n_nodes, 0.0);
        nodes_regionprobs[k].resize(n_nodes);
        for (int n = 0; n < n_nodes; n++) nodes_regionprobs[k][n].clear();
        for (int j = 0; j < n_cells; j++) {
            if (best_attachments[j] < n_nodes && best_attachments[j] >= 0) {
                nodes_regionprobs[k][best_attachments[j]].push_back(1.0 * cells[j].region_counts[k] / cells[j].total_counts);
            }
        }

        if (data.predetermined_region_weights.size() == 0) {
            if (nodes_regionprobs[k][root].size() < std::max(40.0, 0.015 * n_cells)) {
                if (index >= 0) std::cout << "Chain " << std::to_string(index) << ": ";
                std::cout << "In the tree inferred without CNAs, there were not enough cells attached to the root (" << nodes_nbcells[0] << ") to estimate the region weights, so COMPASS could not attempt to find CNAs.." << std::endl;
                return false; // not enough cells attached to the root: cannot find CNAs
            }
            double rootprob = 0;
            for (double prob : nodes_regionprobs[k][root]) {
                rootprob += prob / nodes_regionprobs[k][root].size();
            }
            region_probabilities[k] = rootprob;
        }
    }
    // Select regions whose probability is different between the root and another node
    for (int k = 0; k < n_regions; k++) {
        if (!data.region_is_reliable[k]) continue;
        //If there is prior information that a region is of interest even though we couldn't see any copy number change in the counts,
        //still make it a candidate. The germline variants may be more powerful than the region counts in determining CN
        if (data.region_to_cn_type[k] == CNType::CNT_CNV_DETECTED) {
            candidate_regions.push_back(k);
        }
        else if (data.region_to_cn_type[k] == CNType::CNT_UNKNOWN) { //we skip segments defined as CN NEUTRAL
            for (int n = 1; n < n_nodes; n++) {
                if (n != root) {
                    if (nodes_regionprobs[k][n].size() >= std::max(40.0, 0.03 * n_cells)) {
                        double prob = 0;
                        for (double d : nodes_regionprobs[k][n]) prob += d / nodes_regionprobs[k][n].size();
                        if ((prob > region_probabilities[k] * 1.275 || region_probabilities[k] > prob * 1.35) && (region_probabilities[k] > 0.05 / n_regions || (!parameters.filter_regions))) {
                            candidate_regions.push_back(k);
                            break;
                        }
                    }
                }
            }
        }
    }
    if (candidate_regions.size() == 0) {
        if (index >= 0) std::cout << "Chain " << std::to_string(index) << ": ";
        std::cout << "In the tree inferred without CNAs, there were no regions whose coverage in one node were different from the root, so COMPASS could not identify any CNAs." << std::endl;
        return false;
    }



    regions_successor = std::vector<int>(n_regions,-1);
    for (int k=0;k<candidate_regions.size()-1;k++){
        int k1 = candidate_regions[k];
        int k2 = candidate_regions[k+1];
        if (data.region_to_chromosome.size()>0 && data.region_to_chromosome[k1] == data.region_to_chromosome[k2]){
            regions_successor[k1]=k2;
        }
    }
    if (true || parameters.verbose){
        if (index >=0) std::cout<<"Chain "<<std::to_string(index)<<": ";
        std::cout<<"After the first phase (inferring the best tree without CNAs), COMPASS identified the following candidate regions which might contain CNAs: ";
        for (int i: candidate_regions){
            std::cout<<data.region_to_name[i]<<",";
        }
        std::cout<<std::endl;
    }
    
    return true;
}


void Tree::add_node(int parent){
    // Add a new node below an existing one, and randomly reassign the children of the parent node.
    n_nodes++;
    nodes.push_back(new Node(cache_scores));
    parents.push_back(parent);
    node_probabilities.resize(n_nodes);
    for (int child: children[parent]){
        if (std::rand()%2==0) parents[child] = n_nodes-1;
    }
    compute_children();
}

void Tree::delete_node(int node){
    for (int child: children[node]){
        parents[child] = parents[node];
    }
    delete nodes[node];
    // Change the index of the last node.
    if (node!=n_nodes-1){
        nodes[node] = nodes[n_nodes-1];
        for (int n=0;n<n_nodes;n++){
            if (parents[n]==n_nodes-1) parents[n] = node;
        }
        parents[node] = parents[n_nodes-1];
    }
    n_nodes--;
    nodes.resize(n_nodes);
    parents.resize(n_nodes);
    node_probabilities.resize(n_nodes);
    compute_children();
}


void Tree::prune_reattach(){
    // Prune a subtree, and reattach it somewhere else in the tree
    check_root_cnv();
    if (n_nodes<=1){ // Need at least 2 nodes to be able to prune and reattach.
        hastings_ratio=0.0;
        return;
    }

    // Randomly select one node (apart from the root) to prune and reattach
    int pruned_node_id = ( std::rand() % (n_nodes-1) ) + 1;

    // Find all nodes that are not a descendant of the pruned node
    // Perform a DFT from the pruned node to identify all its descendants
    std::vector<bool> is_descendant(n_nodes,false);
    std::stack<int> stk;
    stk.push(pruned_node_id);
    while (!stk.empty()) {
        int top = stk.top();
        stk.pop();
        for (int child:children[top]) {
            stk.push(child);
        }
        is_descendant[top] = true;
    }
    std::vector<int> candidate_attachment_points;
    for (int i=0;i<n_nodes;i++){
        if (!is_descendant[i])
            candidate_attachment_points.push_back(i);
    }

    // Sample one attachment point among all nodes that are not descendant of the pruned node
    int attachment_point_id = candidate_attachment_points[std::rand()%candidate_attachment_points.size()];

    // Update the tree accordingly, and recompute the genotypes
    parents[pruned_node_id] = attachment_point_id;
    compute_children(); //update the children list according to the new parent vector

    hastings_ratio=1.0;
    check_root_cnv();
}

void Tree::swap_node_labels(){
    check_root_cnv();
    // Select the 2 nodes and copy the nodes below them
    if (n_nodes<=1){
        hastings_ratio=0.0;
        return;
    }
    int node1 = std::rand()%n_nodes;
    int node2 = std::rand()%(n_nodes-1);
    if (node2==node1) node2 = n_nodes-1;

    //if node 1 or 2 was the root, move all germline events and CNVs
    int smallest = std::min(node1, node2);
    if (smallest == 0) {
        Node* pSmallest = nodes[smallest];
        Node* pLargest = nodes[std::max(node1, node2)];
        pSmallest->move_germline_mutations_and_cnvs(pLargest);
    }

    // Exchange the nodes
    Node* temp = nodes[node1];
    nodes[node1] = nodes[node2];
    nodes[node2] = temp;

    
    

    hastings_ratio=1.0;
    check_root_cnv();
}

void Tree::move_SNV(){
    check_root_cnv();

    // sample the source node (from which the event will be moved) among the nodes which have at least one mutation.
    std::vector<int> nodes_with_event{};
    for (int i=0;i<n_nodes;i++){
        //look for non-germline SNP
        for (std::size_t j = 0; j < nodes[i]->get_number_mutations();++j) {
            if (data.locus_to_variant_type[nodes[i]->get_mutations()[j]] != VariantType::VT_GERMLINE) {
                nodes_with_event.push_back(i);
            }
        }
    }
    int initial_nb_nodes_with_events = nodes_with_event.size();
    int new_nb_nodes_with_events = initial_nb_nodes_with_events;
    int source_node = nodes_with_event[std::rand() % nodes_with_event.size()];
    hastings_ratio=1.0 * initial_nb_nodes_with_events * nodes[source_node]->get_number_non_germline_mutations();

    // Select the destination node (either an existing or a new node)
    int destination_node;
    double new_node_prob=0.2;
    if (n_nodes<=2) new_node_prob = 1.0; //we don't want to move them to root, so we need at least 3 nodes
    if (1.0*std::rand()/RAND_MAX<new_node_prob){
        // Move the event to a new node.
        int parent = std::rand()%n_nodes; // Select the parent of the new node
        hastings_ratio*= 1/new_node_prob * n_nodes * std::pow(2,children[parent].size());
        add_node(parent);
        destination_node = n_nodes-1;
    }
    else{
        // Move the event to an existing node
        destination_node = (std::rand()%(n_nodes-2)) + 1; //+1 for skipping root
        if (destination_node==source_node) destination_node = n_nodes-1;
        hastings_ratio*= 1/(1.0-new_node_prob) * (n_nodes-2);
        
    }

    // Remove the mutation from the source node and add it to the destination node
    int mutation = nodes[source_node]->remove_random_mutation();
    nodes[destination_node]->add_mutation(mutation);
    if (nodes[destination_node]->get_number_mutations()==1) new_nb_nodes_with_events++;

    // If some CNAs affect this locus, possibly change the alleles affected by the CNA.
    std::vector<bool> below_mutation(n_nodes,false);
    for (int n: DFT_order){
        if (n==destination_node || (n>0 &&below_mutation[parents[n]]) ) below_mutation[n]=true;
        //may cause a corruption if change_alleles_CNA_locus changes anything!
        std::vector<int> mutsToChange;
        for (auto CNA: nodes[n]->get_CNA_events()){
            if (std::get<0>(CNA)==data.locus_to_region[mutation]){
                mutsToChange.push_back(mutation); //don't change in here, will do insert/erase, may mess up iterators in the loop
            }
        }
        for (int mut : mutsToChange) {
            nodes[n]->change_alleles_CNA_locus(mut, below_mutation[n]);
        }
    }

    // If the source node is now empty, remove it.
    if (source_node !=0 && nodes[source_node]->is_empty()){
        // Delete the source node
        // To reverse the move, we will need to select add new node, select the right parent for the new destination node, and reassign the children correctly.
        hastings_ratio*= new_node_prob / n_nodes/ std::pow(2,children[source_node].size());
        delete_node(source_node);
        //need to decrease the destination node index if the source node was before it
        if (source_node < destination_node) --destination_node;
        new_nb_nodes_with_events--;
    }
    else{
        // The source node still exists. To reverse the move, we need to select not adding a new node, and select the right destination node.
        hastings_ratio*= (1.0-new_node_prob) / (n_nodes-1);
        if (nodes[source_node]->get_number_mutations()==0) new_nb_nodes_with_events--;
    }

    // To reverse the move, we must select the right source node, the right event in this node
    hastings_ratio*= 1.0/new_nb_nodes_with_events / nodes[destination_node]->get_number_mutations();
    check_root_cnv();
}

void Tree::split_merge_node(){
    check_root_cnv();
    double default_merge_probability = 0.20;
    double merge_probability = default_merge_probability;
    if (n_nodes<=1) merge_probability =0.0; //cannot merge nodes if there is only one node.

    if ( (1.0*std::rand())/RAND_MAX <= merge_probability){ // merge nodes
        // Select one node which is not the root
        int node1 = std::rand()%(n_nodes-1) + 1;
        int node2 = parents[node1];
        if ((node1 == 2) && (n_nodes == 3)) {
            int iii = 10;
        }
        // Move the events of node 1 into node 2 
        while (nodes[node1]->get_number_non_germline_mutations()>0){
            nodes[node2]->add_mutation(nodes[node1]->remove_random_mutation());
        }
        while (nodes[node1]->get_number_CNA()>0){
            std::tuple<int,int,std::vector<int>> CNA = nodes[node1]->remove_random_CNA();
            if (node2 != 0) { //don't move copy number events to root. Instead just let them be deleted.
                nodes[node2]->add_CNA(CNA, node2 == 0);
            }
        }

        delete_node(node1);
        if (node1 < node2) --node2;
//        if (node1 != n_nodes - 1 && node2 == n_nodes - 1) node2 = node1; // the index of node 2 might have changed following the deletion.

        // Hastings ratio
        // For merge, there is one possibility to merge the events and set the new parents
        // For split, there are 2**n_events possibilities to split the events between the 2 nodes
        // and 2**n_children possibilities to set the parent of the children of the split node
        if (n_nodes==1) hastings_ratio = 1.0/merge_probability;
        else hastings_ratio = (1.0-merge_probability)/merge_probability;
        
        hastings_ratio *= 1.0/std::pow(2.0,nodes[node2]->get_number_mutations() + nodes[node2]->get_number_CNA());
        hastings_ratio *= 1.0/ std::pow(2.0,children[node2].size()); 
    }
    else{ // Split node
        int node = std::rand()%n_nodes;
        add_node(node); // The new node is a child of the original node, and the children of the original node are randomly distributed between the 2 nodes.

        //Randomly split the events between the two nodes (choosing a random number of events to move)
        if (nodes[node]->get_number_non_germline_mutations()>0){
            int n_events_moved = std::rand()%(nodes[node]->get_number_non_germline_mutations()+1) ; // we can move from 0 to all of the mutations
            for (int i=0;i<n_events_moved;i++){
                nodes[n_nodes-1]->add_mutation(nodes[node]->remove_random_mutation());
            }
        }
        if (nodes[node]->get_number_CNA()>0){
            int n_events_moved = std::rand()%(nodes[node]->get_number_CNA()+1) ;
            for (int i=0;i<n_events_moved;i++){
                std::tuple<int,int,std::vector<int>> CNA = nodes[node]->remove_random_CNA();
                nodes[n_nodes-1]->add_CNA(CNA, n_nodes - 1 == 0);
            }
        }

        // Hastings ratio
        int n_events = nodes[node]->get_number_mutations()+nodes[n_nodes-1]->get_number_mutations()
                        + nodes[node]->get_number_CNA()+nodes[n_nodes-1]->get_number_CNA();
        if (n_nodes==2) hastings_ratio = default_merge_probability;
        else hastings_ratio = merge_probability / (1.0-merge_probability);
        
        hastings_ratio *= std::pow(2.0,n_events);
        hastings_ratio *= std::pow(2.0,children[node].size()+children[n_nodes-1].size()); 
    }
    check_root_cnv();
}


void Tree::add_remove_CNA(bool use_CNA){
    check_root_cnv();
    double default_add_probability=0.70;
    double add_probability = default_add_probability;
    
    std::vector<int> nodes_with_events{};
    for (int i=0;i<n_nodes;i++){
        if (nodes[i]->get_number_CNA()>0) nodes_with_events.push_back(i);
    }
    if (nodes_with_events.size() == 0) {
        add_probability = 1.0; //cannot remove a CNA event if none exists
        //if it is not possible to add or remove CMV without messing up the root, we do nothing.
        if (n_nodes <= 1) {
            hastings_ratio = 0.0;
            return;
        }
    }
    else if (n_nodes <= 1) {
        add_probability = 0.0; //we don't want to add CNVs to the root
    }

    std::vector<int> candidate_regions_CNLOH{};
    for (int k=0;k<n_regions;k++){
        if (data.region_to_loci[k].size()>0 && (data.region_is_reliable[k] || (!parameters.filter_regions_CNLOH) )) candidate_regions_CNLOH.push_back(k);
    }
    if (candidate_regions.size() + candidate_regions_CNLOH.size()==0){
        hastings_ratio=0.0;
        return;
    }
    
    
    if ( (1.0*std::rand())/RAND_MAX <= add_probability){ // add CNA event
        // Select node, type (gain, loss, CNLOH), region and alleles
        int node = (std::rand()%(2*n_nodes-1))+1; // can add the CNA to an existing node, or to a new node below an existing node. And, don't add to the root
        int type=0;
        int n_possible_types = 1;
        if (use_CNA && candidate_regions.size() >0){
            n_possible_types=3;
            type = 1 - (std::rand()%3);
        }
        int region=-1;
        int n_candidate_regions=-1;
        if (type==0){
            region = candidate_regions_CNLOH[std::rand()%candidate_regions_CNLOH.size()];
            n_candidate_regions = candidate_regions_CNLOH.size();
        }
        else{
            region = candidate_regions[std::rand()%candidate_regions.size()];
            n_candidate_regions = candidate_regions.size();
        }
        
        int parent;
        if (node<n_nodes){ // Add event to an existing node
            parent = node;
            hastings_ratio=1.0;
        }
        else{ // Create a new node
            parent = node - n_nodes;
            node = n_nodes;
            add_node(parent);
            hastings_ratio = std::pow(2,children[parent].size()); 
        }
        std::vector<int> alleles{};
        for (int i=0;i<data.region_to_loci[region].size();i++){
            int allele =std::rand()%2; //0: CNA affects ref allele; 1: CNA affects alt allele
            // Can only gain/lose lose one allele if we had at least one copy of it !
            if (allele==0 && nodes[parent]->get_n_ref_allele(data.region_to_loci[region][i])==0) allele=1;
            if (allele==1 && nodes[parent]->get_n_alt_allele(data.region_to_loci[region][i])==0) allele=0;
            alleles.push_back(allele);
        }
        nodes[node]->add_CNA(std::make_tuple(region,type,alleles), node == 0);
        // There are 2*n_nodes-1 possibilities for where to place the CNA, 1 or 3 possibilities for the type, n_candidate_regions possibilities for the regions,
        // 2**n_alleles possibilities for the alleles and in case a new node was created, 2**(n_children of the parent) possibilities to assign the children of the parent
        // To reverse the move, we need to select remove, select the same node, and select the right CNA event
        int n_nodes_with_events = nodes_with_events.size();
        if (nodes[node]->get_number_CNA()==1) n_nodes_with_events++;
        hastings_ratio *= (1.0-default_add_probability) /n_nodes_with_events  / nodes[node]->get_number_CNA() 
                            / add_probability * (2*n_nodes-1) * n_candidate_regions * n_possible_types * std::pow(2,alleles.size());
    }
    else{ // remove CNA event
        // Select a node which has a CNA event
        int node = nodes_with_events[std::rand()%nodes_with_events.size()];
        auto CNA = nodes[node]->remove_random_CNA(); 
        int n_alleles = std::get<2>(CNA).size();
        int n_possible_types=1;
        if (use_CNA && candidate_regions.size()>0) n_possible_types=3;
        int n_candidate_regions= candidate_regions.size();
        if (std::get<1>(CNA)==0) n_candidate_regions = candidate_regions_CNLOH.size();
        hastings_ratio=1.0;
        //not sure how to deal with this, this line was below before, which clearly doesn't work in the case the node gets deleted
        //so, at least better this way
        hastings_ratio *= add_probability / 2.0 / n_nodes / n_candidate_regions / n_possible_types / std::pow(2, n_alleles)
            / (1.0 - add_probability) * nodes_with_events.size() * (nodes[node]->get_number_CNA() + 1);
        if (node!=0 && nodes[node]->is_empty()){
            // The node is now empty--> delete it.
            // Set the parent of the children of the deleted node.
            int parent = parents[node];
            if (node!=n_nodes-1 && parent==n_nodes-1) parent = node;
            delete_node(node);
            hastings_ratio = 1.0/std::pow(2,children[parent].size());
        }
    }
    check_root_cnv();
}


void Tree::move_CNA(){
    check_root_cnv();
    // sample the source node (from which the event will be moved) among the nodes which have at least one CNA event.
    std::vector<int> nodes_with_event{};
    for (int i=0;i<n_nodes;i++){
        if (nodes[i]->get_number_CNA()>0) nodes_with_event.push_back(i);
    }
    int initial_nb_nodes_with_events = nodes_with_event.size();
    int new_nb_nodes_with_events = initial_nb_nodes_with_events;
    if (initial_nb_nodes_with_events==0){ 
        hastings_ratio=0.0;
        return;
    }
    
    double new_node_prob=0.20;
    int source_node = nodes_with_event[std::rand() % nodes_with_event.size()];
    hastings_ratio=1.0 * initial_nb_nodes_with_events * nodes[source_node]->get_number_CNA();
    std::tuple<int,int,std::vector<int>> CNA = nodes[source_node]->remove_random_CNA();
    if (source_node !=0 && nodes[source_node]->is_empty()){
        // Delete the source node
        // To reverse the move, we will need to select add new node, select the right parent for the new destination node, and reassign the children correctly.
        hastings_ratio*= new_node_prob / n_nodes/ std::pow(2,children[source_node].size());
        delete_node(source_node);
        new_nb_nodes_with_events--;
    }
    else{
        // The source node still exists. To reverse the move, we need to select not adding a new node, and select the right destination node.
        hastings_ratio*= (1-new_node_prob) / (n_nodes-1);
        if (nodes[source_node]->get_number_CNA()==0) new_nb_nodes_with_events--;
    }
    
    int destination_node;
    if (n_nodes<=2) new_node_prob = 1.0; //we never move CNVs to base, so if we only have one other node, we must create a new one
    if (1.0*std::rand()/RAND_MAX<new_node_prob){
        // Move the event to a new node.
        int parent = std::rand()%n_nodes; // Select the parent of the new node
        hastings_ratio*= 1.0/new_node_prob * n_nodes * std::pow(2,children[parent].size());
        add_node(parent);
        destination_node = n_nodes-1;
    }
    else{
        // Move the event to an existing node
        destination_node = std::rand()%(n_nodes-2)+1;//skip base, also leave a gap in the end for the case of sampling the source node
        if (destination_node==source_node) destination_node = n_nodes-1;
        hastings_ratio*= 1.0/(1.0-new_node_prob) * (n_nodes-1);
    }
    // Add the event
    nodes[destination_node]->add_CNA(CNA, destination_node == 0);
    if (nodes[destination_node]->get_number_CNA()==1) new_nb_nodes_with_events++;

    // To reverse the move, we must select the right source node, the right event in this node
    hastings_ratio*= 1.0/new_nb_nodes_with_events / nodes[destination_node]->get_number_CNA();
    check_root_cnv();
}

void Tree::merge_or_duplicate_CNA(){
    check_root_cnv();
    // If there are two indentical CNAs in the tree, can merge them into one CNA at their most recent common ancestor
    // If one node containing a CNA has multiple children, can duplicate this CNA and place both copies in parallel branches.

    // Find duplicate CNAs
    std::multiset<std::tuple<int,int,std::vector<int>>> CNAs_in_tree;
    for (int n=0;n<n_nodes;n++){
        for (auto CNA: nodes[n]->get_CNA_events()) CNAs_in_tree.insert(CNA);
    }
    std::vector<std::tuple<int,int,std::vector<int>>> duplicate_CNAs{};
    for (auto CNA: CNAs_in_tree){
        if (CNAs_in_tree.count(CNA)>1){
            // Only add the CNA if it was not already in the vector
            bool new_CNA=true;
            for (auto CNA_already_in_vector: duplicate_CNAs){
                new_CNA = new_CNA && (CNA!=CNA_already_in_vector);
            }
            if (new_CNA) duplicate_CNAs.push_back(CNA);
        }
    }

    // Find nodes containing a CNA event and having multiple children
    std::vector<int> nodes_with_CNA_and_multiple_children{};
    for (int n=0;n<n_nodes;n++){
        if (nodes[n]->get_number_CNA()>0 && children[n].size()>1) nodes_with_CNA_and_multiple_children.push_back(n);
    }

    if (duplicate_CNAs.size()==0 && nodes_with_CNA_and_multiple_children.size()==0){
        // No CNA can be merged or duplicated
        hastings_ratio=0.0;
        return;
    } 


    int event_index = std::rand()%(duplicate_CNAs.size() + nodes_with_CNA_and_multiple_children.size());
    if (event_index < duplicate_CNAs.size()){
        // Merge 2 CNAs
        auto CNA = duplicate_CNAs[event_index];
        std::vector<int> nodes_containing_CNA{};
        for (int n=0;n<n_nodes;n++){
            for (auto CNA_node: nodes[n]->get_CNA_events()){
                if (CNA==CNA_node) nodes_containing_CNA.push_back(n);
            }
        }
        int index1 = std::rand()%nodes_containing_CNA.size();
        int index2 = std::rand()%(nodes_containing_CNA.size()-1);
        if (index1==index2) index2=nodes_containing_CNA.size()-1;
        int node1 = nodes_containing_CNA[index1];
        int node2 = nodes_containing_CNA[index2];
        // Find their most recent common ancestor
        int ancestor = node1;
        while (!is_ancestor(ancestor,node2)) ancestor = parents[ancestor];
        nodes[node1]->remove_CNA(CNA);
        nodes[node2]->remove_CNA(CNA);
        int new_n_nodes_with_CNA_and_multiple_children = nodes_with_CNA_and_multiple_children.size();
        if ((ancestor == 0) || (std::rand()%2==0)){
            add_node(ancestor);
            nodes[n_nodes-1]->add_CNA(CNA, n_nodes-1 == 0);
            new_n_nodes_with_CNA_and_multiple_children+=1;
        }
        else{
            nodes[ancestor]->add_CNA(CNA, ancestor == 0);
            if (nodes[ancestor]->get_number_CNA()==1) new_n_nodes_with_CNA_and_multiple_children+=1;
        }

        // Hastings ratio
        // In order to merge, we had to:
        //  * Select merge (probability: duplicate_CNAs.size() / (duplicate_CNAs.size() + nodes_with_CNA_and_multiple_children.size()))
        //  * Select the CNA among duplicate_CNAs.size() possibilities
        //  * Select the 2 nodes containing this CNA: nodes_containing_CNA.size() choose 2 possibilities
        //  * Select whether to move the CNA to their MRCA or create a new node below it (2 possibilities)
        // In order to reverse this move, we would have to:
        //  * Select duplicate (new_n_nodes_with_CNA_and_multiple_children / (duplicate_CNAs.size()-1 + new_n_nodes_with_CNA_and_multiple_children))
        //       (new_n_nodes_with_CNA_and_multiple_children because after having performed the move, there might be one more node with CNA and multiple children)
        //  * Select the right node among nodes_with_CNA_and_multiple_children_new.size() possibilities
        //  * Select the right CNA among nodes[ancestor]->get_number_CNA() possibilities
        //  * Select the right 2 children: children[ancestor].size() possibilities
        //  * For each child, select the right descendant to put the event.

        int n_descendants1=0; // Number of descendants in the subtree below "ancestor" containing node1
        int node=node1;
        while (node!=ancestor){
            node = parents[node];
            n_descendants1+=1;
        }
        std::stack<int> stk;
        stk.push(node1);
        while (!stk.empty()) {
            int top = stk.top();
            stk.pop();
            for (int child: children[top]) stk.push(child);
            if (top!=node1) n_descendants1++;
        }
        int n_descendants2=0;
        node=node2;
        while (node!=ancestor){
            node = parents[node];
            n_descendants2+=1;
        }
        stk.push(node2);
        while (!stk.empty()) {
            int top = stk.top();
            stk.pop();
            for (int child: children[top]) stk.push(child);
            if (top!=node2) n_descendants2++;
        }
        // Select merge
        hastings_ratio= 1.0 / duplicate_CNAs.size() * (duplicate_CNAs.size() + nodes_with_CNA_and_multiple_children.size());
        // Select which CNAs to merge, the nodes containing that CNA and whether to put the CNA to the MRCA or a new node below it.
        hastings_ratio*= duplicate_CNAs.size() * std::exp(cache_scores->log_n_choose_k(nodes_containing_CNA.size(),2) *2) *2.0;
        // Reverse: Select duplicate
        hastings_ratio*= 1.0* new_n_nodes_with_CNA_and_multiple_children / (duplicate_CNAs.size()-1 + new_n_nodes_with_CNA_and_multiple_children);
        // Select the node, the CNA and the descendants that will get a copy of the CNA
        hastings_ratio*= 1.0 / new_n_nodes_with_CNA_and_multiple_children / nodes[ancestor]->get_number_CNA() / children[ancestor].size() / n_descendants1 / n_descendants2;
    }
    else{
        // Duplicate 1 CNA

        //Select the node and the CNA event to duplicate
        int node_ancestor = nodes_with_CNA_and_multiple_children[std::rand()%nodes_with_CNA_and_multiple_children.size()];
        std::vector<std::tuple<int,int,std::vector<int>>> CNAs_node{};
        for (auto CNA: nodes[node_ancestor]->get_CNA_events()) CNAs_node.push_back(CNA);
        auto CNA = CNAs_node[std::rand()%CNAs_node.size()];

        // Select the 2 nodes where to add the CNA
        int child1=children[node_ancestor][std::rand()%children[node_ancestor].size()];
        int child2=children[node_ancestor][std::rand()%children[node_ancestor].size()];
        while (child1==child2) child2=children[node_ancestor][std::rand()%children[node_ancestor].size()];

        std::vector<int> descendants1{};
        std::stack<int> stk;
        stk.push(child1);
        while (!stk.empty()) {
            int top = stk.top();
            stk.pop();
            for (int child: children[top]) stk.push(child);
            descendants1.push_back(top);
        }
        std::vector<int> descendants2{};
        stk.push(child2);
        while (!stk.empty()) {
            int top = stk.top();
            stk.pop();
            for (int child: children[top]) stk.push(child);
            descendants2.push_back(top);
        }
        int node1 = descendants1[std::rand()%descendants1.size()];
        int node2 = descendants2[std::rand()%descendants2.size()];

        nodes[node_ancestor]->remove_CNA(CNA);
        nodes[node1]->add_CNA(CNA, node1 == 0);
        nodes[node2]->add_CNA(CNA, node2 == 0);

        int new_n_nodes_with_CNA_and_multiple_children = nodes_with_CNA_and_multiple_children.size();
        if (nodes[node_ancestor]->get_number_CNA()==0) new_n_nodes_with_CNA_and_multiple_children--;

        // Hastings Ratio
        std::vector<int> nodes_containing_CNA{};
        for (int n=0;n<n_nodes;n++){
            for (auto CNA_node: nodes[n]->get_CNA_events()){
                if (CNA==CNA_node) nodes_containing_CNA.push_back(n);
            }
        }

        hastings_ratio= 1.0 / nodes_with_CNA_and_multiple_children.size() * (duplicate_CNAs.size() + nodes_with_CNA_and_multiple_children.size());
        hastings_ratio*= 
        hastings_ratio*= 1.0 * nodes_with_CNA_and_multiple_children.size() * (nodes[node_ancestor]->get_number_CNA()+1.0)
                             * children[node_ancestor].size() *descendants1.size() * descendants2.size();
        //Reverse
        hastings_ratio*= 1.0 *(duplicate_CNAs.size()+1.0) / (duplicate_CNAs.size()+1 + new_n_nodes_with_CNA_and_multiple_children);
        hastings_ratio*=  1.0 / (duplicate_CNAs.size()+1.0) / std::exp(cache_scores->log_n_choose_k(nodes_containing_CNA.size(),2)) /2.0;

        //wait with deleting node to last - will mess up indices otherwise
        if (node_ancestor > 0 && nodes[node_ancestor]->is_empty()) delete_node(node_ancestor);
    }
    check_root_cnv();
}

void Tree::exchange_Loss_CNLOH(){
    check_root_cnv();
    std::vector<int> nodes_with_event{};
    for (int n=0;n<n_nodes;n++){
        if (nodes[n]->get_number_LOH()>0) nodes_with_event.push_back(n);
    }
    if (nodes_with_event.size()==0) hastings_ratio=0.0;
    else{
        int node = nodes_with_event[std::rand() % nodes_with_event.size()];
        hastings_ratio = nodes[node]->exchange_Loss_CNLOH(candidate_regions);
    }
    check_root_cnv();
}

void Tree::change_alleles_CNA(){
    check_root_cnv();
    hastings_ratio=1.0;
    std::vector<int> nodes_with_event{};
    for (int n=0;n<n_nodes;n++){
        if (nodes[n]->get_number_CNA_mut()>0) nodes_with_event.push_back(n);
    }
    if (nodes_with_event.size()==0) hastings_ratio=0.0;
    else{
        int node = nodes_with_event[std::rand() % nodes_with_event.size()];
        nodes[node]->change_alleles_CNA();
    }
    check_root_cnv();
}

void Tree::exchange_Loss_Double_loss() {
    check_root_cnv();
    std::vector<int> nodes_with_event{};
    std::vector<std::vector<int>> node_regions{};


    //we only allow this if there is an LOH/double and there is not two LOH in any lineage where the node participates
    //this holds for both flips, since a LOH + double LOH should never happen, heavily penalized in prior
    for (int n = 0; n < n_nodes; n++) {
        std::vector<int> ok_CNA_segments;
        for (auto& CNA: nodes[n]->get_CNA_events()) {
            if (std::get<1>(CNA) < 0) {
                std::set<int> forbidden = rec_get_nodes_in_lineages_with_2_CNA(std::get<0>(CNA));
                if (forbidden.count(n) == 0) {
                    ok_CNA_segments.push_back(std::get<0>(CNA));
                    break;
                }
            }
        }
        if (!ok_CNA_segments.empty()) {
            nodes_with_event.push_back(n);
            node_regions.push_back(ok_CNA_segments);
        }
    }
    if (nodes_with_event.size() == 0) hastings_ratio = 0.0;
    else {
        int ind = std::rand() % nodes_with_event.size();
        hastings_ratio = nodes[nodes_with_event[ind]]->exchange_Loss_Double_loss(node_regions[ind]);
    }
    check_root_cnv();
}

double Tree::get_regionprobs_variance(){
    double mean = 0;
    int n_reliable_regions=0;
    for (int k=0;k<n_regions;k++){
        if (data.region_is_reliable[k]){
            n_reliable_regions++;
            mean+=region_probabilities[k] * n_regions;
        }
    }
    mean = mean / n_reliable_regions;
    std::cout<<"mean "<<mean<<std::endl;
    double variance=0;
    for (int k=0;k<n_regions;k++){
        if (data.region_is_reliable[k]){
            std::cout<<k<<": "<<region_probabilities[k] * n_regions<<std::endl;
            variance+= std::pow(region_probabilities[k]*n_regions - mean ,2);
        }
    }
    variance = variance/n_reliable_regions;
    return variance;
}


void Tree::add_node_no_randomization(int parent) {
    // Add a new node below an existing one, and randomly reassign the children of the parent node.
    n_nodes++;
    nodes.push_back(new Node(cache_scores));
    parents.push_back(parent);
    node_probabilities.resize(n_nodes);
    compute_children();
}

void Tree::init_debug_tree(bool use_CNA_arg) {
    cells_attach_loglik.resize(n_cells);
    cells_loglik.resize(n_cells);
    cells_attach_prob.resize(n_cells);
    best_attachments.resize(n_cells);
    use_CNA = false; //not sure why these tricks are needed
    compute_children();
    compute_likelihood(true);
    if (use_CNA_arg && select_regions()) {
        use_CNA = true;
        compute_likelihood(true);
    }
    compute_prior_score();
    update_full_score();
}
