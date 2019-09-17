/*
* @Author: guanja
* @Date:   2019-07-04 17:19:08
* @Last Modified by:   guanja
* @Last Modified time: 2019-09-17 11:43:08
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

#include "../../C/utils/utils.cpp"

#include "../../C/methods/sinimin_fwer.cpp"



int main(int argc, char** argv)
{

  // input files.
  std::string edg_file = "../../examples/sim_data/sim_ps_0.05_pcon_0.05_simID_0_edges.txt";
  std::string map_file = "../../examples/sim_data/sim_ps_0.05_pcon_0.05_simID_0_snp_map.txt";
  std::string dat_file = "../../examples/sim_data/sim_ps_0.05_pcon_0.05_simID_0_X.txt";
  std::string lab_file = "../../examples/sim_data/sim_ps_0.05_pcon_0.05_simID_0_Y.txt";
  std::string cov_file = "../../examples/sim_data/sim_ps_0.05_pcon_0.05_simID_0_C.txt";
  std::string snp_file = "../../examples/sim_data/sim_ps_0.05_pcon_0.05_simID_0_snpID.txt";

  int maxlen = 1;
  int N_PERM = 1000;
  std::string out_file="sim_ps_0.05_pcon_0.05_simID_0_fwer";
  std::string tarone_file = out_file + std::string("_tarone.txt");
  std::string frequency_file = out_file + std::string("_minpv_frequencies.txt");

  bool encode_or = true;
  double thresh = 2.53571e-05;

  // Read the data.
  Data dataset(dat_file, lab_file, snp_file, cov_file);
  
  // Read the edges.                                      
  Edges edges(edg_file);

  // Read the mapping from SNPs to genes.
  Mapping mapping(map_file);

  // translate the mapping to the SNP-IDs that are provided with the data set.
  mapping.translate_snp_ids(dataset.snp_ids);

  // Init the edge-epistasis method.
  siniminFWER sinimin_fwer(dataset, edges, mapping, thresh, maxlen, 
                           N_PERM, encode_or);
  // Process the edges.
  sinimin.process_edges();
  
  // Write the tarone summary.
  sinimin.tarone.write_summary(out_file);

}
