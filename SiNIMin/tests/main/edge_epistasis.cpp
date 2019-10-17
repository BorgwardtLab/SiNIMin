/*
* @Author: guanja
* @Date:   2019-07-04 17:19:08
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-31 19:59:44
*/

// Include standard libs.
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>

// Include C++ data classes.
#include <deque>
#include <vector>
#include <set>
#include <tuple>
#include <unordered_map>



#include "../../C/struct/data.cpp"
#include "../../C/struct/edges.cpp"
#include "../../C/struct/mapping.cpp"
// #include "../../C/struct/tarone_cmh.cpp"
// #include "../../C/struct/tarone_cmh_wy.cpp"
#include "../../C/utils/utils.cpp"

#include "../../C/methods/edge_epistasis_nowy.cpp"



int main(int argc, char** argv)
{

  // input files.
  std::string edg_file = "../../examples/data/edges.txt";
  std::string map_file = "../../examples/data/mapping_file.txt";
  std::string dat_file = "../../examples/data/data.txt";
  std::string lab_file = "../../examples/data/labels.txt";
  std::string cov_file = "../../examples/data/covar.txt";
  std::string snp_file = "../../examples/data/snp_ids.txt";

  int maxlen = 0;
  std::string out_file="test.txt";
  bool encode_or=false;

  // Read the data.
  Data dataset(dat_file, lab_file, snp_file, cov_file);
  
  // Read the edges.                                      
  Edges edges(edg_file);

  // Read the mapping from SNPs to genes.
  Mapping mapping(map_file);

  // translate the mapping to the SNP-IDs that are provided with the data set.
  mapping.translate_snp_ids(dataset.snp_ids);

  // Init the edge-epistasis method.
  EdgeEpistasis edge_epistasis(dataset, edges, mapping, 0.05, maxlen,
                               out_file, encode_or);

  // Process the edges.
  edge_epistasis.process_edges();

}
