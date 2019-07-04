/*
* @Author: guanja
* @Date:   2019-07-04 17:19:08
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-04 19:14:29
*/

// Include standard libs.
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>

// Include C++ data classes.
#include <vector>
#include <set>


#include "../struct/edges.cpp"
#include "../struct/mapping.cpp"
#include "../utils/utils.cpp"



int main(int argc, char** argv)
{

  std::string edge_file = "../examples/data/edges.txt";
  std::string map_file = "../examples/data/mapping_file.txt";

  std::vector<std::string> snp_order = {"rs1", "rs2", "rs3", "rs4", "rs5",
                                        "rs6", "rs7", "rs8", "rs9", "rs10"};

  Edges edge_class(edge_file);
  Mapping mapping(map_file);

  std::cout << edge_class.n_nodes << std::endl;
  std::cout << edge_class.n_edges << std::endl;

  print_nested_vector(edge_class.edges_str);
  print_nested_vector(edge_class.edges_int);

  print_vector(snp_order);

  // translate the mapping.
  mapping.translate_snp_ids(snp_order);

  print_map_vector(mapping.geneview_map_idx);



}