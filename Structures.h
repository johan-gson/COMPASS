#ifndef DEF_STRUCTURES
#define DEF_STRUCTURES
#include <vector>
#include <map>
#include <array>
#include <cstdint>

enum CellType {
    CT_U = 0, //unknown
    CT_N = 1, //normal, non-malignant
    CT_M = 2 //malignant
};

enum VariantType {
    VT_UNKNOWN = 0, //unknown
    VT_GERMLINE = 1, //normal, non-malignant
    VT_SOMATIC = 2 //malignant
};

enum VariantImpact {
    VI_UNKNOWN = 0,
    VI_DRIVER = 1,
    VI_PASSENGER = 2,
    VI_GERMLINE = 3
};

enum CNType {
    CNT_UNKNOWN = 0, //unknown, if we have no WGS data
    CNT_COPY_NEUTRAL = 1, //for example from WGS data
    CNT_CNV_DETECTED = 2 //for example from WGS data
};

enum CNAllelePrior {
    CNAP_UNKNOWN = 0,
    CNAP_INCREASED = 1, //This means that the AF is higher in malignant cells, which can be used to initialize the CNA variant vector
    CNAP_DECREASED = 2 //Like INCREASED, but AF is lower in malignant cells
};

struct Cell{
    std::vector<int> ref_counts; // number of ref reads for each variable locus
    std::vector<int> alt_counts; // number of alt reads for each variable locus
    std::vector<int> region_counts; // number of reads in each region
    std::vector<int> amplicon_counts; //number of reads in each amplicon
    std::vector<int> genotypes;
    std::vector<int> GQ;
    CellType cell_type;
    std::string name;
    int total_counts; // sum of the read counts in each region, i.e., total reads counts for the cell
    std::size_t sample = 0;
};

struct Data{
    std::string sex; // "male" (1 chr X) or "female" (2 chr X)
    std::vector<std::string> locus_to_chromosome;
    std::vector<int> locus_to_position;
    std::vector<std::string> locus_to_reference;
    std::vector<std::string> locus_to_alternative;
    std::vector<bool> variant_is_SNV; // true if the variant is a SNV, false otherwise (indel)
    std::vector<double> locus_to_freq; // frequency of the variant in the population
    std::vector<std::string> locus_to_name; //name of the variant, which describes the effect on the protein (if any)
    std::vector<int> locus_to_region;
    std::vector <VariantType> locus_to_variant_type;
    std::vector <VariantImpact> locus_to_variant_impact;
    std::vector <CNAllelePrior> locus_to_cna_allele_prior;
    std::vector<std::string> sample_ids;
    
    std::vector<std::vector<int>> region_to_loci; //list of loci on this amplicon (possibly empty)
    std::vector<std::string> region_to_chromosome;
    std::vector<std::string> region_to_name;
    std::vector<bool> region_is_reliable; // true for the amplicons that are used for computing the CNV score.
    std::vector<std::vector<double>> predetermined_region_weights; // [region][sample]probability for a read to fall in each of the regions, when these values are given as input and not inferred (otherwise they will be inferred by using the cells attached at the root)
    std::vector<CNType> region_to_cn_type;
    std::vector<std::vector<double>> region_to_theta; //[region][sample]
    std::vector<std::string> amplicon_to_name;
    std::vector<std::vector<double>> amplicon_to_theta;//[region][sample]
    std::vector<std::vector<double>> amplicon_weights;//[region][sample]
    std::vector<int> amplicon_to_region;
    std::vector<std::vector<std::size_t>> region_to_amplicons;
};

struct Params{
    double sequencing_error_rate;
    double omega_hom; //concentration parameter for beta-binomial in the SNV likelihood, when homozygous
    double omega_het; //concentration parameter for beta-binomial in the SNV likelihood, when heterozygous

    double sequencing_error_rate_indel; // more errors for indels than SNVs
    double omega_hom_indel;
    double omega_het_indel;

    double prior_dropoutrate_mean; // beta prior for the dropout rates. the mean is alpha/(alpha+beta)
    double prior_dropoutrate_omega; // beta prior for the dropout rates. the inverse overdispersion parameter omega is alpha+beta

    //double theta; // inverse dispersion of the negative binomial for the CNV likelihood

    double doublet_rate;
    
    bool use_doublets;
    bool filter_regions;
    bool filter_regions_CNLOH;
    bool verbose;

    // Penalties in the tree prior
    double node_cost;
    double CNA_cost;
    double LOH_cost;
    double mut_notAtRoot_cost;
    double mut_notAtRoot_freq_cost;
    double germline_root_bonus = 1000; //this shouldn't really be needed
    double somatic_non_root_bonus = 1000;
    double normal_not_at_root_penalty = 100;
    double malignant_at_root_penalty = 50;
    double cnloh_removing_unknown_somatic_mut_penalty = 500;
    double loss_removing_driver_somatic_mut_penalty = 20000;
    double cna_in_multiple_branches_penalty = 1000;
    double two_CNA_in_lineage_penalty = 100; //we want these to be well supported for them to exist
    bool allow_double_allele_loss = true;
    double min_exp_cn = 0.03; //used for double CN loss
    double small_node_penalty = 600; //penalty for nodes that have very few cells and do not represent branches (have more than one child)
};

//allele_changes describes the changes to allele 1 and 2. A CNLOH can thus be -1,1 or 1,-1, 
//a biallelic loss -1,-1, etc. A gain is typically 1,0 or 0,1, whole genome doubling would in theory be 1,1.
//variant_alleles describes for each variant the index in allele_changes representing the allele where the variant is.
//variant_alleles typically stay fixed when flipping which allele is lost, that is changed in allele_changes.
//for variants that are uncertain (i.e., when we don't know how they are correlated with the other variants within the segment),
//the variant_alleles may be shuffled as part of an mcmc move.
//Note that the alleles are not comparable across CNAs, allele 1 in one may physically be allele 2 in another CNA, we just don't know this.
struct CNADesc {
    CNADesc() {}
    CNADesc(std::array<int, 2> allele_changes_) {
        allele_changes = allele_changes_;
    }
    std::array<int, 2> allele_changes{ 0, 0 };
};

struct NodeDefinition {
    std::vector<std::size_t> variants; //indices to the variants in data. 
    std::map<std::size_t, CNADesc> CNAs; //only contains the CNAs that are changed
};

struct TreeDefinition {
    std::vector<NodeDefinition> nodes;
    std::vector<std::size_t> parents;
    std::vector<std::size_t> segment_variant_alleles;//This is a vector over all variants
};



#endif