#ifndef DEF_TREE
#define DEF_TREE

#include <vector>
#include <set>
#include <string>
#include "Scores.h"
#include "Node.h"
#include "Structures.h"
#include <algorithm>

class Tree{
    //private:
    public: //temp debug
        // Tree structure
        bool use_CNA;
        int n_nodes;
        std::vector<Node*> nodes;
        std::vector<Node*> doublets; 
        std::vector<int> parents; //parents[n]: parent of node n. The root is always at index 0.
        std::vector<std::vector<int>> children; // children[n]: vector of the children of node n
        std::vector<int> DFT_order;

        // Parameters
        std::vector<double> node_probabilities; //prior attachment probabilities to each node.
        std::vector<double> dropout_rates; // dropout rate for each locus.
        std::vector<double> dropout_rates_ref; // dropout rate of the reference allele
        std::vector<double> dropout_rates_alt; // dropout rate of the alternative allele
        std::vector<double> region_probabilities; // probability for a read to fall on each region, when there are 2 copies of each region

        std::vector<int> candidate_regions; // list of amplicons which might be affected by a CNV.
        std::vector<int> regions_successor; // amplicons_successor[k] is l if l is the next candidate amplicon and in the same chromosome as k,
                                          // or -1 if the next amplicon is not on the same chromosome

        // Utils
        Scores* cache_scores;
        std::vector<std::vector<double>> cells_attach_loglik; //index j,k gives the complete (hidden) log likelihood of the data and cell j being attached to node k
        std::vector<double> cells_loglik; // index j gives the observed log likelihood for cell j (marginalizing over attachment points)
        std::vector<std::vector<double>> cells_attach_prob; //index j,k gives the posterior probability that cell j is attached to node k
        std::vector<int> best_attachments;
        std::vector<double> nodes_attachment_counts; //for each node, the expected number of cells attached to it
        double avg_diff_nodeprob; // how much the node probabilities changed between 2 iterations of the EM algorithm
        double avg_diff_dropoutrates; // how much the dropout rates changed between 2 iterations of the EM algorithm


    public:
        double hastings_ratio; // probability to move from previous tree to this tree / probability to move from this tree to the previous tree
        double log_prior_score;
        double log_likelihood;
        double log_score; //complete score of the tree, including the prior and the likelihood

        Tree(Scores* cache, bool use_CNA, bool noRandomization = false); // constructor
        Tree(); // empty constructor
        Tree(const Tree& source); // copy constructor (deep copy)
        Tree(std::string gv_file, bool use_CNA=true); // create tree from graphviz file
        ~Tree();
        Tree& operator=(const Tree&); // Assignment operator: deep copy
        
        void compute_children(); // compute the list of children of each node, from the parent vector. Also compute the DFT order
        void compute_nodes_genotypes(); // update the genotypes of the nodes below the node given in argument. (maybe just give node index as argument...)
        
        

        bool is_ancestor(int potential_ancestor, int potential_descendant); 
        double rec_check_max_one_event_per_region_per_lineage(int node = 0, std::vector<std::vector<int>> previousEvents = std::vector<std::vector<int>>());
        int rec_get_number_of_cnloh_removing_somatic_mut(int node = 0, std::vector<int> muts = std::vector<int>()) const;
        std::set<int> rec_get_nodes_in_lineages_with_2_CNA(int region, int node = 0, std::set<int> parents = std::set<int>(), 
            std::vector<int> previousEvents = std::vector<int>()) const;
        std::set<int> rec_get_descendents(int node) const;


        void compute_attachment_scores(bool use_doublets_local,bool recompute_CNA_scores);
        void compute_likelihood(bool allow_diff_dropoutrates=true);
        void compute_prior_score();
        void update_full_score();


        //utils for computing likelihood
        void compute_cells_likelihoods(bool use_doublets_local); //update the vectors cells_attach_loglik and cells_loglik
        void compute_cells_likelihoods_singlets();
        void compute_cells_likelihoods_doublets();
        std::vector<std::vector<double>> get_cells_likelihoods_doublets();
        void EM_step(bool use_doublets_local, bool allow_diff_dropoutrates); // compute the attachment probabilities of each cell (E step) and update the node probabilities and dropout rates (M step)


        void to_dot(std::string filename, bool simplified); //save the tree structure to the dot format (for visualization)

        void find_CNA();
        void allow_CNA(){use_CNA=true;}
        bool select_regions(int index=-1);
        bool contains_candidate_regions(){return candidate_regions.size()>0;}

        //MCMC moves
        void add_node (int parent); // Add a node below an existing one, and randomly assign the children of the parent to the new node or their initial parent
        void delete_node(int node); // delete a node, setting the parent of its children to the node's parent, and adapting the node indices
        void prune_reattach();
        void swap_node_labels();
        //void add_remove_mutation();
        void move_SNV();
        void split_merge_node();
        void add_remove_CNA(bool use_CNA);
        void move_CNA();
        void merge_or_duplicate_CNA();
        void exchange_Loss_CNLOH();
        void change_alleles_CNA();
        void exchange_Loss_Double_loss();

        double get_regionprobs_variance();

        int get_num_nodes() const { return n_nodes; }
        void check_root_cnv() {
            if (nodes[0]->get_CNA_events().size() > 0) {
                std::cout << "Warning: CNA in the root!!!\n";
            }
            //check that all germline variants are in root
            const auto& muts = nodes[0]->get_mutations();
            for (std::size_t i = 0; i < data.locus_to_variant_type.size(); ++i) {
                if (data.locus_to_variant_type[i] == VariantType::VT_GERMLINE) {
                    if (std::find(muts.begin(), muts.end(), i) == muts.end()) {
                        std::cout << "Warning: Germline event not in root!!!\n";
                    }
                }
            }
        }

        //for debugging
        Node* get_node(int i) { return nodes[i]; }
        void add_node_no_randomization(int parent);
        void init_debug_tree(bool use_CNA_arg);
};

#endif