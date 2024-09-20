#ifndef DEF_NODE
#define DEF_NODE

#include <vector>
#include <set>
#include <string>

#include "Structures.h"
#include "Scores.h"

extern Data data;

class Node {
    private: 
        std::vector<int> mutations; //somatic mutations always transform from ref to alt. 


        std::set<std::tuple<int,int,std::vector<int>>> CNA_events; // triplet (region_index,type,affected alleles)
        // the second element is +1 in case of copy number gain, -1 in case of copy number loss, and 0 in case of copy neutral event
        // the third element is a vector of length the number of loci in this region (potentially 0), 
        // and it contains 0 if the ref allele is deleted (Loss) / duplicated (Gain) / lost (copy neutral), or 1 for the alt allele
        // I use an (ordered) set, with the default (lexicographic) order, which here sorts the events by position.
        
        std::set<int> affected_loci; // loci for which the copy number of the ref or alt alleles is changed in the node.
        std::set<int> affected_regions; // regions affected by a CNA event (except copy neutral) in the node
        
        std::vector<int> n_ref_allele; //array of number of copies of the ref allele, for each variable locus
        std::vector<int> n_alt_allele; //array of number of copies of the alt allele, for each variable locus
        std::vector<int> cn_regions; // copy number of each region

        Scores* cache_scores;


    public:
        std::vector<double> attachment_scores_SNV;
        std::vector<double> attachment_scores_CNA;
        std::vector<double> attachment_scores;

        Node(Scores* cache); //constructor
        Node(Node& source); //copy constructor
        Node(Node& source1,Node& source2); // create doublet from 2 nodes
        void init_structures();
        ~Node();

        // Compute genotype of the node, depending on the genotype of its parent
        void update_genotype(Node* parent, bool isroot);
        // compute attachment scores from scratch
        void compute_attachment_scores(bool use_CNA,const std::vector<double>& dropout_rates_ref,
                                        const std::vector<double>& dropout_rates_alt, const std::vector <std::vector<double>>& region_probabilities);
        // just compute difference with parent
        void compute_attachment_scores_parent(bool use_CNA,Node* parent,const std::vector<double>& dropout_rates_ref,
                                        const std::vector<double>& dropout_rates_alt, const std::vector <std::vector<double>>& region_probabilities,bool recompute_CNA_scores);

        // MCMC moves for the nodes
        void add_mutation(int locus){mutations.push_back(locus);}
        int remove_random_mutation(); // removes a random mutation, and return the index of the mutation
        void add_CNA(std::tuple<int,int,std::vector<int>> CNA, bool isRoot, bool silent = false){
            //if adding a loss when there already is a loss of the same segment, change to a double loss to keep the tree more readable
            bool inserted = false;
            if (std::get<1>(CNA) == -1) {
                for (auto CNA2 : CNA_events) {
                    if ((std::get<1>(CNA2) == -1) && (std::get<0>(CNA2) == std::get<0>(CNA2))) {
                        remove_CNA(CNA2);
                        add_CNA(std::make_tuple(std::get<0>(CNA2), -2, std::get<2>(CNA2)), isRoot);
                        inserted = true;
                        break; //important to break here, the loop is messed up when inserting/removing
                    }
                }
            }
            if (!inserted) {
                CNA_events.insert(CNA);
            }
            if (isRoot && !silent) {
                std::cout << "adding CNV to base\n";
            }
        }
        std::tuple<int,int,std::vector<int>> remove_random_CNA();
        void remove_CNA(std::tuple<int,int,std::vector<int>> CNA) { CNA_events.erase(CNA); }
        //void remove_CNAs_in_region(int region);
        double exchange_Loss_CNLOH(std::vector<int> candidate_regions);
        void change_alleles_CNA();
        void change_alleles_CNA_locus(int locus,bool heterozygous);
        double exchange_Loss_Double_loss(std::vector<int> candidate_regions);
        
        
        // Various accessors
        bool is_empty(){ return (mutations.size()==0 && CNA_events.size()==0);}

        //  Accessors for Mutations
        int get_n_ref_allele(int locus){return n_ref_allele[locus];}
        int get_n_alt_allele(int locus){return n_alt_allele[locus];}
        const std::vector<int>& get_mutations() {return mutations;}
        std::size_t get_number_mutations() {return mutations.size();}
        int get_number_non_germline_mutations() {
            int n = 0;
            for (std::size_t i = 0; i < mutations.size(); ++i) {
                if (data.locus_to_variant_type[mutations[i]] != VariantType::VT_GERMLINE) {
                    ++n;
                }
            }
            return n;
        }

        void move_germline_mutations_and_cnvs(Node* pNewRoot);

        //  Accessors for CNAs
        int get_cn_region(int region){return cn_regions[region];}
        std::set<std::tuple<int,int,std::vector<int>>> get_CNA_events() {return CNA_events;}
        std::set<int> get_affected_regions(){return affected_regions;}
        std::size_t get_number_CNA() {return CNA_events.size();}
        int get_number_CNA_noncopyneutral() {
            int count=0;
            for (auto CNA: CNA_events){
                if (std::get<1>(CNA)!=0) count++;
            }
            return count;
        }
        int get_number_CNA_mut() {
            int count=0;
            for (auto CNA: CNA_events){
                if (std::get<2>(CNA).size()>0) count+=1;
            }
            return count;
        }
        int get_number_CN_losses(){
            int n_CN_losses = 0;
            for (auto CNA: CNA_events){
                if (std::get<1>(CNA)<0) n_CN_losses++; 
            }
            return n_CN_losses;
        }
        int get_number_LOH(){
            int n_LOH = 0;
            for (auto CNA: CNA_events){
                if (std::get<1>(CNA)<=0 && std::get<2>(CNA).size()>0) n_LOH++; 
            }
            return n_LOH;
        }
        int get_number_effective_LOH(Node* parent){ //Count the number of losses and CNLOH, where the parent was heterozygous for at least one allele.
            int n_LOH = 0;
            for (auto CNA: CNA_events){
                if (std::get<1>(CNA)<=0 && std::get<2>(CNA).size()>0){
                    bool parent_het = false;
                    for (int i=0;i<std::get<2>(CNA).size();i++){
                        if (parent->get_n_alt_allele(data.region_to_loci[std::get<0>(CNA)][i])>0 && parent->get_n_ref_allele(data.region_to_loci[std::get<0>(CNA)][i])>0){
                            parent_het=true;
                        }
                    }
                    if (parent_het) n_LOH++;
                }
            }
            return n_LOH;
        }
        int get_number_disjoint_CNA(std::vector<int> regions_successor);
        int get_number_disjoint_LOH(std::vector<int> regions_successor); // only count copy number loss of regions containing a variant
        
        
        // Label of the node, used for plotting.
        std::string get_label();
        std::string get_label_simple(std::set<int> excluded_mutations);

};

#endif