#ifndef __YAML_SERIALIZATION_H
#define __YAML_SERIALIZATION_H
#include <string>
#include "structures.h"
//#include "types.h"

 
//Import/export of trees to graphviz format

TreeDefinition read_yaml(std::string filename, const Data& data);
void write_yaml(std::string filename, TreeDefinition tree, const Data& data);

 
#endif
