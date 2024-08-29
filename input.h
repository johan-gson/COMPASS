#ifndef DEF_INPUT
#define DEF_INPUT

#include <string>
#include <vector>

void load_CSV(std::string base_name, std::string regionweights_file,bool use_CNA=true);
void filter_regions();
void init_params();
void fit_region_thetas();

#endif