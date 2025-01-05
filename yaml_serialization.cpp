#include "yaml_serialization.h"
#include <fkYAML/node.hpp>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>


// overloads must be defined in the same namespace as user-defined types.
void from_node(const fkyaml::node& node, TreeDefinition& tree) {
    //assume the tree definition is empty
    
    //parents
    for (const auto& val : node["parents"]) {
        auto v = val.get_value<int>();
        tree.parents.push_back(std::size_t(v));
    }

    //nodes
    for (const auto& n : node["nodes"]) {
        NodeDefinition nd;
        auto x = n["CNAs"];
        for (const auto& CNA : n["CNAs"]) {
            CNADesc cd;
            cd.allele_changes[0] = std::size_t(CNA["allele_changes_1"].get_value<int>());
            cd.allele_changes[1] = std::size_t(CNA["allele_changes_2"].get_value<int>());
            auto ind = std::size_t(CNA["segment_index"].get_value<int>());
            nd.CNAs[ind] = cd;
        }
        for (const auto& v : n["variants"]) {
            nd.variants.push_back(std::size_t(v.get_value<int>()));
        }
        tree.nodes.push_back(nd);
    }


    //segment variant alleles
    for (const auto& val : node["segment_variant_alleles"]) {
        auto v = val.get_value<int>();
        tree.segment_variant_alleles.push_back(std::size_t(v));
    }
}


void to_node(fkyaml::node& node, const TreeDefinition& tree) {
    node = fkyaml::node{ {"nodes", fkyaml::node::sequence()},{"parents", fkyaml::node::sequence()},{"segment_variant_alleles", fkyaml::node::sequence()} };
    auto& parents = node["parents"].get_value_ref<fkyaml::node::sequence_type&>();
    auto& nodes = node["nodes"].get_value_ref<fkyaml::node::sequence_type&>();
    auto& segment_variant_alleles = node["segment_variant_alleles"].get_value_ref<fkyaml::node::sequence_type&>();

    for (auto parent : tree.parents) {
        parents.push_back(parent);
    }

    for (auto n : tree.nodes) {
        nodes.push_back(fkyaml::node::mapping());
        auto& node_ref = nodes.back().get_value_ref<fkyaml::node::mapping_type&>();
        
        node_ref["variants"] = fkyaml::node::sequence();
        auto& variants_ref = node_ref["variants"].get_value_ref<fkyaml::node::sequence_type&>();
        for (auto v : n.variants) {
            variants_ref.push_back(v);
        }
        
        node_ref["CNAs"] = fkyaml::node::sequence();
        auto& CNAs_ref = node_ref["CNAs"].get_value_ref<fkyaml::node::sequence_type&>();
        for (auto& c : n.CNAs) {
            CNAs_ref.push_back(fkyaml::node::mapping());
            auto& CNA_ref = CNAs_ref.back().get_value_ref<fkyaml::node::mapping_type&>();
            CNA_ref["segment_index"] = c.first;
            CNA_ref["allele_changes_1"] = c.second.allele_changes[0];
            CNA_ref["allele_changes_2"] = c.second.allele_changes[1];
        }
    }
    
    for (auto a : tree.segment_variant_alleles) {
        segment_variant_alleles.push_back(a);
    }
}


void write_yaml(std::string filename, TreeDefinition tree, const Data& data) {
    fkyaml::node yaml_tree{ {"tree", tree } };
    std::string yaml_str = fkyaml::node::serialize(yaml_tree);

    std::ofstream os(filename, std::ios::binary); //binary: avoid createing \r\n line endings on Windows
    if (!os.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    os << yaml_str;
    os.close();
    int ii = 13;
}

TreeDefinition read_yaml(std::string filename, const Data& data) {
    std::ifstream is(filename, std::ios::binary);
    fkyaml::node node = fkyaml::node::deserialize(is);
    fkyaml::node tree_node = node["tree"];
    TreeDefinition tree = tree_node.get_value<TreeDefinition>();
    int ii = 13;
    return tree;
}


/*

// Serialization and deserialization helpers
void serializeCNADesc(fkyaml::node & yaml, const CNADesc& cna) {
    yaml["allele_changes"] = fkyaml::array(cna.allele_changes.begin(), cna.allele_changes.end());
    yaml["variant_alleles"] = fkyaml::array(cna.variant_alleles.begin(), cna.variant_alleles.end());
}

CNADesc deserializeCNADesc(const fk::Yaml& yaml) {
    CNADesc cna;
    cna.allele_changes = {
        yaml["allele_changes"][0].as<int>(),
        yaml["allele_changes"][1].as<int>()
    };
    for (const auto& val : yaml["variant_alleles"]) {
        cna.variant_alleles.push_back(val.as<std::size_t>());
    }
    return cna;
}

void serializeNodeDefinition(fk::Yaml& yaml, const NodeDefinition& node) {
    yaml["variants"] = fk::Yaml::array(node.variants.begin(), node.variants.end());

    fk::Yaml cnaMap;
    for (const auto& [key, value] : node.CNAs) {
        fk::Yaml cnaYaml;
        serializeCNADesc(cnaYaml, value);
        cnaMap[std::to_string(key)] = cnaYaml;
    }
    yaml["CNAs"] = cnaMap;
}

NodeDefinition deserializeNodeDefinition(const fk::Yaml& yaml) {
    NodeDefinition node;

    for (const auto& val : yaml["variants"]) {
        node.variants.push_back(val.as<std::size_t>());
    }

    for (const auto& [key, value] : yaml["CNAs"]) {
        std::size_t cnaKey = std::stoul(key);
        node.CNAs[cnaKey] = deserializeCNADesc(value);
    }

    return node;
}

void serializeTreeDefinition(fk::Yaml& yaml, const TreeDefinition& tree) {
    fk::Yaml nodesYaml;
    for (const auto& node : tree.nodes) {
        fk::Yaml nodeYaml;
        serializeNodeDefinition(nodeYaml, node);
        nodesYaml.push_back(nodeYaml);
    }

    yaml["nodes"] = nodesYaml;
    yaml["parents"] = fk::Yaml::array(tree.parents.begin(), tree.parents.end());
}

TreeDefinition deserializeTreeDefinition(const fk::Yaml& yaml) {
    TreeDefinition tree;

    for (const auto& nodeYaml : yaml["nodes"]) {
        tree.nodes.push_back(deserializeNodeDefinition(nodeYaml));
    }

    for (const auto& val : yaml["parents"]) {
        tree.parents.push_back(val.as<std::size_t>());
    }

    return tree;
}

// Read TreeDefinition from a YAML file
TreeDefinition readTreeFromYAML(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    file.close();

    fk::Yaml root = fk::Yaml::parse(buffer.str());
    return deserializeTreeDefinition(root);
}

// Write TreeDefinition to a YAML file
void writeTreeToYAML(const TreeDefinition& tree, const std::string& filename) {
}

*/






