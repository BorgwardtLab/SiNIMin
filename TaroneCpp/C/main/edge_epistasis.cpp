/*
* @Author: guanja
* @Date:   2019-07-04 17:19:08
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-08 19:02:55
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

  // Read the edges.                                      
  Edges edge_class(edge_file);

  // Read the mapping from SNPs to genes.
  Mapping mapping(map_file);

  // translate the mapping to the SNP-IDs that are provided with the data set.
  mapping.translate_snp_ids(snp_order);


  // Iterate the edges.

  /* For every edge, enumerate the intervals in a breadth-first, i.e.

    1. Enumerate all intervals of length '1' and check, if they are testable.
    2. Grow those intervals that were testable and compute their supports by
       combination, i.e. for interval (tau, length), combine its support with 
       the support of interval (tau-1, length). (Not sure whether this is the
       best approach though, maybe brute-force would also work). Using an
       unordered map with tuple (start, length) as key, support (complete, not 
       per table, as this would not work for WY permutations) as value
       would be good.
    3. For every edge, check whether intervals were already enumerated.
    4. Do a pairwise test of the supports of all possible interval combinations.
  */




}