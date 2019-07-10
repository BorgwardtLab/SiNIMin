/*
* @Author: guanja
* @Date:   2019-07-04 17:19:08
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-10 16:26:03
*/

// Include standard libs.
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>

// Include C++ data classes.
#include <vector>
#include <set>

// This has to be included, otherwise getopt cannot be found.
#include <getopt.h>


#include "../struct/data.cpp"
#include "../struct/edges.cpp"
#include "../struct/mapping.cpp"
#include "../struct/tarone_cmh.cpp"
#include "../struct/tarone_cmh_wy.cpp"
#include "../struct/timer.cpp"
#include "../utils/utils.cpp"
#include "../utils/output.cpp"

#include "../methods/edge_epistasis_nowy.cpp"



int main(int argc, char** argv)
{
  /*
    Read command-line input.
  */
  char opt;

  // output filenames.
  std::string pvalue_file;
  std::string profiling_file;
  std::string tarone_file;
  std::string significant_file;

  // input filenames.
  std::string data_file;
  std::string snp_file;
  std::string covar_file;
  std::string labels_file;
  std::string map_file;
  std::string edge_file;
  std::string out_prefix;

  // Boolean flags (to check whether all input files have been specified).
  bool flag_input_file = false;
  bool flag_snp_file = false;
  bool flag_covar_file = false;
  bool flag_labels_file = false;
  bool flag_map_file = false;
  bool flag_edge_file = false;
  bool flag_fwer = false;

  // Default values for FWER and number of threads to use for analysis.
  double target_fwer = 0.05;
  int n_threads = 1;

  while( (opt = getopt(argc, argv, "i:s:c:l:m:e:f:n:o:")) != -1){
    switch(opt){
      case 'i':
        data_file = std::string(optarg); 
        flag_input_file = true;
        break;
      case 's':
        snp_file = std::string(optarg); 
        flag_snp_file = true;
        break;
      case 'c':
        covar_file = std::string(optarg); 
        flag_covar_file = true;
        break;
      case 'l':
        labels_file = std::string(optarg); 
        flag_labels_file = true;
        break;
      case 'm':
        map_file = std::string(optarg); 
        flag_map_file = true;
        break;
      case 'e':
        edge_file = std::string(optarg); 
        flag_edge_file = true;
        break;
      case 'f':
        target_fwer = atof(optarg); 
        flag_fwer = true;
        break;
      case 'n':
        n_threads = atof(optarg); 
        break;
      case 'o':
        out_prefix = std::string(optarg);

        pvalue_file = out_prefix + std::string("_pvalues.csv");
        profiling_file = out_prefix + std::string("_profiling.txt");
        tarone_file = out_prefix + std::string("_tarone.txt");
        significant_file = out_prefix + std::string("_significant.txt");
        break;
      default:
        std::cerr << "Error @ main: Incorrect arguments!!!" << std::endl;
        exit(-1);
    }
  }

  if (!flag_input_file || !flag_labels_file || !flag_edge_file || 
      !flag_snp_file || !flag_map_file)
  {
    std::cerr << "Error @ main: Minimum input missing." << std::endl;
    // print_options();
    exit(-1);
  }


  // Check if the output files can be generated.
  check_out_file(pvalue_file);
  check_out_file(profiling_file);
  check_out_file(tarone_file);


  /*
    Set number of threads. 
    REMARK: For small data sets there is an overhead of using high numbers of 
    threads.
  */
  omp_set_num_threads(n_threads);
  std::cout << "Using " << n_threads << " thread(s)." << std::endl << std::endl;

  // Init timer.
  Timer timer(n_threads);

  // start measuring the complete execution time, and the initialization time.
  double tic_global = measureTime();
  double tic_init = measureTime();
  std::cout << "Reading data. " << std::endl; 


  // ---------------------------------------------------------------------------
  // Initialize the data to an empty constructor, load it with the given files,
  // either with or without covariate files.
  Data dataset;
  if (flag_covar_file) 
  {
    Data tmp_data(data_file, labels_file, snp_file, covar_file);
    dataset = tmp_data;
  }else{
    Data tmp_data(data_file, labels_file, snp_file);
    dataset = tmp_data;
  }
    
  // Read the edges.                                      
  Edges edges(edge_file);

  // Read the mapping from SNPs to genes.
  Mapping mapping(map_file);

  // translate the mapping to the SNP-IDs that are provided with the data set.
  mapping.translate_snp_ids(dataset.snp_ids);

  // Stop timing the input.
  timer.io = measureTime() - tic_init;
  std::cout << "Time elapsed: " << timer.io << " sec." << std::endl; 
  std::cout << std::endl;


  // ---------------------------------------------------------------------------
  // Init the edge-epistasis method.
  EdgeEpistasis edge_epistasis(dataset, edges, mapping, target_fwer, 0, 
                               pvalue_file);

  // Process the edges.
  double tic_mine = measureTime();
  std::cout << "Running EdgeEpistasis." << std::endl;
  edge_epistasis.process_edges();
  timer.proc_edge = measureTime() - tic_mine;
  std::cout << "Time elapsed: " << timer.proc_edge;
  std::cout << " sec." << std::endl << std::endl;

  // Filter the significant edges.
  int n_significant = \
      write_significant_edge_epi(edge_epistasis.tarone.corr_threshold(),
                                 pvalue_file, significant_file);

  // Write the tarone summary.
  edge_epistasis.tarone.write_summary(tarone_file, n_significant);

  timer.total = measureTime() - tic_global;
  timer.write_profiling(profiling_file);

  
  
}