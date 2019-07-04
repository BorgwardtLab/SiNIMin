/*
* @Author: guanja
* @Date:   2019-07-04 17:19:08
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-04 18:29:17
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
#include "../utils/utils.cpp"



int main(int argc, char** argv)
{

  std::string edge_file = "../examples/data/edges.txt";

  Edges edge_class(edge_file);

  std::cout << edge_class.n_nodes << std::endl;
  std::cout << edge_class.n_edges << std::endl;

  print_nested_vector(edge_class.edges_str);
  print_nested_vector(edge_class.edges_int);



}