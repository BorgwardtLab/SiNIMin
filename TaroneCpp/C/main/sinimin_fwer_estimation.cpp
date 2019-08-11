/*
* @Author: guanja
* @Date:   2019-07-04 17:19:08
* @Last Modified by:   guanja
* @Last Modified time: 2019-08-11 17:05:32
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
#include "../struct/timer.cpp"
#include "../utils/utils.cpp"
#include "../utils/output.cpp"

#include "../methods/sinimin_fwer.cpp"



int main(int argc, char** argv)
{
  /*
    Read command-line input.
  */
  char opt;

  // output filenames.
  std::string fwer_file;
  std::string profiling_file;

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

  // Default values for FWER and number of threads to use for analysis.
  double threshold = 0.05;
  int n_threads = 1;
  bool encode_or = true;

  // Number of permutations, default: 0 (no WY permutations).
  int n_perm = 1000;

  // Maximal length (dimension) of intervals, default: 0.
  int max_interval_dim = 0;

  while( (opt = getopt(argc, argv, "i:s:c:l:m:e:n:p:d:a:o:t:r")) != -1){
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
      case 't':
        threshold = atof(optarg); 
        break;
      case 'n':
        n_threads = atof(optarg); 
        break;
      case 'p':
        n_perm = atof(optarg); 
        break;
      case 'd':
        max_interval_dim = atof(optarg); 
        break;
      case 'r':
        encode_or = false; 
        break;
      case 'o':
        out_prefix = std::string(optarg);

        fwer_file = out_prefix + std::string("_fwer.txt");
        profiling_file = out_prefix + std::string("_fwer_profiling.txt");
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
  check_out_file(fwer_file);
  check_out_file(profiling_file);

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
  
  // ---------------------------------------------------------------------------
  // Start data input.

  double tic_init = measureTime();
  std::cout << "Reading data. " << std::endl; 
  // Read the edges.                                      
  Edges edges(edge_file);

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
    
  // Read the mapping from SNPs to genes.
  Mapping mapping(map_file);

  // translate the mapping to the SNP-IDs that are provided with the data set.
  mapping.translate_snp_ids(dataset.snp_ids);

  // Stop timing the input.
  timer.io = measureTime() - tic_init;
  std::cout << "Time elapsed: " << timer.io << " sec." << std::endl; 
  std::cout << std::endl;


  // ---------------------------------------------------------------------------
  // Run the edge-epistasis method.

  double tic_mine = measureTime();
  std::cout << "Estimation of FWER from " << n_perm << " permutations.";
  std::cout << std::endl;

    
  if (n_perm < 0)
  {
    std::cerr << "Invalid argument for number of permutations.";
    std::exit(-1);
  }else{
    // If n_perm is set to a value larger than 0, WY-permutations will be run.
    siniminFWER sinimin(dataset, edges, mapping, threshold, 
                        max_interval_dim, n_perm, encode_or);
    sinimin.process_edges();
    sinimin.tarone.write_summary(fwer_file);

    timer.proc_edge = measureTime() - tic_mine;
    std::cout << "Time elapsed: " << timer.proc_edge;
    std::cout << " sec." << std::endl << std::endl;
  }

  
  

  timer.total = measureTime() - tic_global;
  timer.write_profiling(profiling_file);

}