/*
* @Author: guanja
* @Date:   2019-07-09 14:13:20
* @Last Modified by:   guanja
* @Last Modified time: 2019-09-17 11:29:45
*/


#ifndef _sinimin_wy_func_cpp_
#define _sinimin_wy_func_cpp_

#include <iostream>


// Include C++ data classes.
#include <deque>
#include <vector>
#include <set>
#include <tuple>
#include <unordered_map>

// data structures.
#include "../struct/data.cpp"
#include "../struct/edges.cpp"
#include "../struct/mapping.cpp"
#include "../struct/interval.hpp"
#include "../struct/support.hpp"
#include "../struct/tarone_cmh_wy.cpp"
#include "../utils/utils.cpp"


// Interval map, with keys indicating the starting position and length, values
// correspond to the full support of the interval.
typedef std::unordered_map< Interval, Eigen::MatrixXd> interval_map;


/*
  Main class to do the edge-epistasis analysis.
*/
class SiniminWY
{

  public:

    // maximum size of intervals to consider.
    int max_length;

    // The data set.
    Data data;

    // The edges representing the network;
    Edges edges;

    // The mapping, containing the SNP-to-gene mapping.
    Mapping mapping;

    // The Tarone object.
    TaroneCMHwy tarone;

    // The output file.
    std::string output_file;

    // Whether or not interaction should be computed as the binary OR
    bool encode_or;

    // -------------------------------------------------------------------------
    // unordered map that contains the minimum p-value, whether or not
    // the support is prunable, and the p-value.
    // vector: [min-pv, prunable (bool), pval]
    std::unordered_map<std::vector<bool>, Support> support_map;

    // count the number of support-evaluations.
    int supp_eval_count = 0;
    int n_patterns = 0;
    // -------------------------------------------------------------------------

    // Default constructor.
    SiniminWY() = default;

    // Constructor with output filename.
    SiniminWY(Data data_obj, Edges edge_obj, Mapping map_obj, 
              double alpha, int maxlen, int n_perm, std::string filename,
              bool do_encode_or);

    // Functions to run the mining.
    void process_edges();

  private:

    // Function to call generation of all intervals, using only gene name.
    interval_map make_gene_intervals(std::string gene_name);

    // Compute the supports of all intervals within a gene in depth-first 
    // search.
    interval_map intervals_depth_first(const Eigen::MatrixXd& mat);

    // Test combinations of intervals between two adjacent genes.
    void test_interval_combinations(interval_map gene_0_itvl,
                                    interval_map gene_1_itvl,
                                    std::string gene_0_name,
                                    std::string gene_1_name,
                                    std::ofstream& out_file);
};


/*
  Constructor.
*/
SiniminWY::SiniminWY(Data data_obj, Edges edge_obj, Mapping map_obj, 
                     double alpha, int maxlen=0, int n_perm=1000,
                     std::string out_file="./edge_epistasis.csv",
                     bool do_encode_or=true) 
{

  // Set input objects.
  max_length = maxlen;
  data = data_obj;
  edges = edge_obj;
  mapping = map_obj;
  encode_or = do_encode_or;

  // initialize the Tarone object.
  TaroneCMHwy tmp_tarone(alpha, data.pt_samples, data.pt_cases, n_perm);
  tarone = tmp_tarone;

  // init the output filestream.
  output_file = out_file;

  // write output message to indicate encoding of interactions.
  if (encode_or)
  {
    std::cout << "Encoding interactions with OR" << std::endl;  
  }else{
    std::cout << "Encoding interactions with AND" << std::endl;  
  }
}


/*
  Main function, processes edge-by-edge and tests all intervals in two 
  interacting genes.
*/
void SiniminWY::process_edges()
{

  // Unordered map to store the intervals. Each key will correspond to a gene,
  // each value will contain an object 'interval_map', which is again an
  // unordered map, that contains as keys the intervals (string-formatted: 
  // "start_length"), and a matrix representing the binary meta-vector of that
  // interval.
  // std::unordered_map<int, interval_map> gene_intervals;

  // This can be replaced by a vector of invervals, and we can access them by
  // index.
  std::vector<interval_map> gene_intervals(edges.n_nodes);

  // Initialize a vector of zeros that will keep track of whether or not
  // intervals at that gene have been enumerated or not.
  std::vector<int> intervals_enumerated(edges.n_nodes);

  // Generate the output file.
  std::ofstream out_stream(output_file);
  out_stream << "p-value,min_pval,gene0,start0_len0,gene1,start1_len1"; 
  out_stream << std::endl;

  std::cout << "here 1" << std::endl;

  for (int i=0; i<edges.n_edges; i++)
  {

    std::cout << "@ edge " << i+1 << " of " << edges.n_edges << "\r";

    // Get the gene string-names of the genes adjacent to the edge.
    std::string gene_0_str = edges.edges_str[i][0];
    std::string gene_1_str = edges.edges_str[i][1];

    // Get the gene integer-names of the genes adjacent to the edge.
    int gene_0_int = edges.edges_int[i][0];
    int gene_1_int = edges.edges_int[i][1];

    // check if both, source and sink nodes have SNPs overlapping with them.
    if (mapping.geneview_map_idx.find(gene_0_str) == \
        mapping.geneview_map_idx.end())
    {
      intervals_enumerated[gene_0_int] = 9;
      continue;
    }
    if (mapping.geneview_map_idx.find(gene_1_str) == \
        mapping.geneview_map_idx.end())
    {
      intervals_enumerated[gene_1_int] = 9;
      continue;
    }

    // If the intervals for gene0 have not been enumerated yet, do so.
    if (intervals_enumerated[gene_0_int] == 0)
    {
      gene_intervals[gene_0_int] = make_gene_intervals(gene_0_str);
      intervals_enumerated[gene_0_int] = 1;
    }

    // If the intervals for gene1 have not been enumerated yet, do so.
    if (intervals_enumerated[gene_1_int] == 0)
    {
      gene_intervals[gene_1_int] = make_gene_intervals(gene_1_str);
      intervals_enumerated[gene_1_int] = 1;
    }

    std::cout << "here 2" << std::endl;

    // Do the pairwise testing of all intervals in source/sink.
    test_interval_combinations(gene_intervals[gene_0_int],
                               gene_intervals[gene_1_int],
                               gene_0_str, gene_1_str,
                               out_stream);

  }

  std::cout << edges.n_edges << " edges processed. Finishing. " << std::endl;
  std::cout << " Number of processed supports =  " << supp_eval_count << std::endl;
  std::cout << " Number of all supports =  " << n_patterns << std::endl;
  std::cout << " Size of supp_map =  " << support_map.size() << std::endl;
}




/*
  Creates all intervals in gene given by 'gene_name'.
  Starts by extracting the corresponding SNPs and extracting them from the
  data matrix, then calls function for enumeration on the subsampled data 
  matrix.

  TODO: translate the keys to correspond to the actual SNP-pointers. Right
  now, they are relative to the start of the tmp_matrix. (Problem: I have to
  check if this works without problems when combining intervals later on, when
  increasing the length of intervals, as I rely on the distance being '1').
*/
interval_map SiniminWY::make_gene_intervals(std::string gene_name)
{

  // The SNP-ids.
  std::vector<int> snp_ids = mapping.geneview_map_idx[gene_name];

  // The data-matrix corresponding to the SNPs (see comments below for update).
  Eigen::MatrixXd tmp_matrix =  \
      data.matrix.block(snp_ids.front(), 0, snp_ids.size(), data.n_samples); 

  interval_map gene_supports = intervals_depth_first(tmp_matrix);

  return gene_supports;
}


/*
  Enumerates all intervals in matrix mat in a depth-first order, i.e. for
  every starting position, all possible length-intervals are generated.
  Using a depth-first enumeration of the intervals we can stop growing
  intervals once they become untestable and we can enumerate only closed
  intervals, i.e. those that actually change support when grown.
*/
interval_map SiniminWY::intervals_depth_first(const Eigen::MatrixXd& mat)
{

  interval_map intervals;

  // Iterate over all starting positions tau.
  for (int tau=0; tau<mat.rows(); tau++)
  {
    // Init the support to 0.
    Eigen::MatrixXd support = Eigen::MatrixXd::Zero(1, mat.row(tau).size());

    // Init the length of the first interval to 0.
    int len = 0;
    
    // Init the sum of the support.
    int supp_count = support.sum();

    while (tau+len < mat.rows())
    {

      // If the maximum length is reached, stop enumeration of intervals.
      if (max_length > 0 && len == max_length)
      {
        break;
      }
      
      support = binary_or(mat.row(tau+len), support);
      len += 1;
      
      // check if the support has changed from one length to the next. Only 
      // report, if it did.
      if (supp_count == support.sum())
      {
        continue;
      }else{
        supp_count = support.sum();
      }

      // Compute the per-table support and minimum p-value.
      Eigen::VectorXd pt_support = \
        tarone.compute_per_table_support(support);
      // double min_pv = tarone.compute_minpval(pt_support);

      // Check if the support is prunable. If it is, we do not enumerate any
      // super-intervals.
      if (!tarone.is_prunable(pt_support)) 
      {
        Interval tmp_int = Interval(tau, len);
        intervals[tmp_int] = support;
      }
    }
  }

  return intervals;
}


/*
  Tests all combinations of intervals between gene0 and gene1.
  We explore the intervals in a depth-first manner, as this allows us to 
  exploit monotonicities in the support-space.
*/
void SiniminWY::test_interval_combinations(interval_map gene_0_itvl,
                                           interval_map gene_1_itvl,
                                           std::string gene_0_name,
                                           std::string gene_1_name,
                                           std::ofstream& out_file)
{

  // Define the supports of the interval_0 and interval_1 and their 
  // combination.
  Eigen::MatrixXd supp0;
  Eigen::MatrixXd supp1;
  Eigen::MatrixXd support;

  // define the per-table support.
  Eigen::VectorXd pt_support;

  // Get the number of SNPs per gene
  int n_snps_0 = mapping.geneview_map_idx[gene_0_name].size();
  int n_snps_1 = mapping.geneview_map_idx[gene_1_name].size();

  // Start to iterate the intervals in a depth-first search.
  for (int tau_0=0; tau_0<=n_snps_0; tau_0 ++ )
  {
    // Init vector to keep track at which length the previous length-interval
    // was prunable with the next interval.
    std::vector<int> prunable(n_snps_1, n_snps_1+1);

    for (int len_0=1; len_0<=n_snps_0-tau_0; len_0 ++)
    {
      // Create the current interval.
      Interval interval_0 = Interval(tau_0, len_0);

      // Continue if it is not present in the current list of intervals.
      if (gene_0_itvl.find(interval_0) == gene_0_itvl.end())
      {
        continue;
      }

      // Get the support of interval_0 from the map.
      supp0 = gene_0_itvl[interval_0];
      
      // Iterate over the intervals in the second gene.
      for (int tau_1=0; tau_1<n_snps_1; tau_1++)
      {
        for (int len_1=1; len_1<=n_snps_1-tau_1; len_1++)
        {

          // Check if a sub-interval of (tau_0, len_0) has been found prunable
          // with the current interval (tau_1, len_1). Sub-interval means with
          // same starting point, but different length.
          if (len_1 >= prunable[tau_1])
          {
            break;
          }

          // Create the current interval.
          Interval interval_1 = Interval(tau_1, len_1);

          // Continue if it is not present in the current list of intervals.
          if (gene_1_itvl.find(interval_1) == gene_1_itvl.end())
          {
            continue;
          }

          supp1 = gene_1_itvl[interval_1];

          // Compute the binary support of the pattern:
          if (encode_or == true)
          {
            support = binary_or(supp0, supp1);  
          }else{
            support = binary_and(supp0, supp1);  
          }

          // Create a boolean vector (faster hashing in map).
          std::vector<bool> supp_bool = make_bool_vec(support);

          n_patterns += 1;

          // -------------------------------------------------------------------
          // If the current support has not been observed yet, process it.
          if (support_map.find(supp_bool) == support_map.end())
          {
            supp_eval_count += 1;
            // compute the per-table support.
            pt_support = tarone.compute_per_table_support(support);

            // compute the minimum p-value.
            double min_pv = tarone.compute_minpval(pt_support);

            // compute the envelope.
            double evlpe = tarone.envelope(pt_support);

            // Check if the current support is prunable.
            double pvalue = 1.0;
            if (!tarone.is_prunable(evlpe)) 
            {
              // Compute the p-value.
              int a = tarone.compute_supported_cases(support, data.labels);
              pvalue = tarone.compute_pval(a, pt_support);
            }

            // create the support object and store it in the support map.
            Support tmp_sup(min_pv, evlpe, pvalue);
            support_map[supp_bool] = tmp_sup;

            // Now, for the WY permutations, we only need to do the 
            // permutations for patterns we have not seen so far.
            if (tarone.is_testable(min_pv))
            {
              tarone.process_testable(pt_support, support);
            }
          }


          // -------------------------------------------------------------------
          // Compute the per-table support and corresponding minimum p-value.
          Support sup_summary = support_map[supp_bool];

          if (!tarone.is_prunable(sup_summary.envelope())) 
          {
          // Write the interval combination to the output directory.
          out_file << sup_summary.pvalue() << ",";
          out_file << sup_summary.min_pvalue() << "," << gene_0_name << ",";
          out_file << interval_0 << "," << gene_1_name << "," << interval_1;
          out_file << std::endl;
          }else{
          // If the pattern is prunable, we mark the length in the prunable
          // vector at that starting point for the current interval_0 with
          // starting point in tau_0.
          prunable[tau_1] = len_1;
          }
        }
      }
    }
  }
}


#endif


/*
  COMMENTS:

  * If the correct version of Eigen were installed, the slicing of the 
    tmp_matrix could be done like this:
    dataset.matrix(Eigen::seqN(snp_ids.front(), snp_ids.size()), Eigen::all ) 

*/