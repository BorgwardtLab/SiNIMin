/*
* @Author: guanja
* @Date:   2019-07-04 17:19:08
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-09 09:50:05
*/

// Include standard libs.
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>

// Include C++ data classes.
#include <vector>
#include <set>


#include "../../C/struct/edges.cpp"
#include "../../C/struct/mapping.cpp"
#include "../../C/utils/utils.cpp"



int main(int argc, char** argv)
{

  std::string edge_file = "../../examples/data/edges.txt";
  std::string map_file = "../../examples/data/mapping_file.txt";

  std::vector<std::string> snp_order = {"rs1", "rs2", "rs3", "rs4", "rs5",
                                        "rs6", "rs7", "rs8", "rs9", "rs10"};

  std::cout << "-----------------------------------" << std::endl;                                        
  std::cout << "EDGES" << std::endl << std::endl;

  Edges edge_class(edge_file);
  std::cout << "# nodes: " << edge_class.n_nodes << std::endl;
  std::cout << "# edges: " << edge_class.n_edges << std::endl;
  std::cout << "1. Egdes with string names ";
  std::cout << "(no duplicates/symmetries, no self-loops)." << std::endl;
  print_nested_vector(edge_class.edges_str);
  std::cout << "2. Egdes with integer names " << std::endl;;
  print_nested_vector(edge_class.edges_int);

  std::cout << "-----------------------------------" << std::endl;
  std::cout << "MAPPING:" << std::endl << std::endl;
  Mapping mapping(map_file);
  std::cout << "1. Mapping with SNP-IDs." << std::endl;
  print_map_vector(mapping.geneview_map_str);

  // translate the mapping.
  mapping.translate_snp_ids(snp_order);
  std::cout << "2. SNP ordering in data set." << std::endl;
  print_vector(snp_order);
  std::cout << "3. Mapping with SNP-indices in ordering vector." << std::endl;
  print_map_vector(mapping.geneview_map_idx);

}