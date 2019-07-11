/*
* @Author: guanja
* @Date:   2019-07-09 14:13:20
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-11 13:15:37
*/


#ifndef _edge_epistasis_wy_cpp_
#define _edge_epistasis_wy_cpp_

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
#include "../struct/tarone_cmh_wy.cpp"
#include "../utils/utils.cpp"


// Interval map, with keys indicating the starting position and length, values
// correspond to the full support of the interval.
typedef std::unordered_map< std::string, Eigen::MatrixXd> interval_supports;

// Tuple of form (interval start, interval length), for faster enumeration of
// intervals.
typedef std::deque<std::tuple<int, int>> interval_queue;


/*
  Main class to do the edge-epistasis analysis.
*/
class EdgeEpistasisWY
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

    // Default constructor.
    EdgeEpistasisWY() = default;

    // Constructor with output filename.
    EdgeEpistasisWY(Data data_obj, Edges edge_obj, Mapping map_obj, 
                    double alpha, int maxlen, int n_perm, std::string filename);

    // Functions to run the mining.
    void process_edges();

  private:

    // Function to call generation of all intervals, using only gene name.
    interval_supports make_gene_intervals(std::string gene_name);

    // Function to call search for all testable intervals.
    interval_supports find_intervals(const Eigen::MatrixXd& mat);

    // Function to explore all length-1 intervals and their supports for a gene.
    void intervals_layer_1(const Eigen::MatrixXd& mat,
                           interval_supports& intervals, 
                           interval_queue& testable_queue);

    // Function to explore intervals of arbitrary length > 1 and their supports.
    void intervals_layer_x(interval_supports& intervals, 
                           interval_queue& testable_queue);

    // Test combinations of intervals between two adjacent genes.
    void test_interval_combinations(interval_supports gene_0_itvl,
                                    interval_supports gene_1_itvl,
                                    std::string gene_0_name,
                                    std::string gene_1_name,
                                    std::ofstream& out_file);
};


/*
  Constructor.
*/
EdgeEpistasisWY::EdgeEpistasisWY(Data data_obj, Edges edge_obj, Mapping map_obj, 
                                 double alpha, int maxlen=0, int n_perm=1000,
                                 std::string out_file="./edge_epistasis.csv") 
{
  // Set input objects.
  max_length = maxlen;
  data = data_obj;
  edges = edge_obj;
  mapping = map_obj;

  // initialize the Tarone object.
  TaroneCMHwy tmp_tarone(alpha, data.pt_samples, data.pt_cases, n_perm);
  tarone = tmp_tarone;

  // init the output filestream.
  output_file = out_file;
}


/*
  Main function, processes edge-by-edge and tests all intervals in two 
  interacting genes.
*/
void EdgeEpistasisWY::process_edges()
{

  // Unordered map to store the intervals. Each key will correspond to a gene,
  // each value will contain an object 'interval_supports', which is again an
  // unordered map, that contains as keys the intervals (string-formatted: 
  // "start_length"), and a matrix representing the binary meta-vector of that
  // interval.
  std::unordered_map<int, interval_supports> gene_intervals;

  // Generate the output file.
  std::ofstream out_stream(output_file);
  out_stream << "p-value,min_pval,gene0,start0_len0,gene1,start1_len1"; 
  out_stream << std::endl;

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
      continue;
    }
    if (mapping.geneview_map_idx.find(gene_1_str) == \
        mapping.geneview_map_idx.end())
    {
      continue;
    }

    // Check if the intervals for gene_0 and gene_1 have already been
    // created. If not, do so.
    if (gene_intervals.find(gene_0_int) == gene_intervals.end())
    {
      gene_intervals[gene_0_int] = make_gene_intervals(gene_0_str);
    }

    if (gene_intervals.find(gene_1_int) == gene_intervals.end())
    {
      gene_intervals[gene_1_int] = make_gene_intervals(gene_1_str);
    }

    // Do the pairwise testing of all intervals in source/sink.
    test_interval_combinations(gene_intervals[gene_0_int],
                               gene_intervals[gene_1_int],
                               gene_0_str, gene_1_str,
                               out_stream);

  }

  std::cout << edges.n_edges << " edges processed. Finishing. " << std::endl;
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
interval_supports EdgeEpistasisWY::make_gene_intervals(std::string gene_name)
{

  // The SNP-ids.
  std::vector<int> snp_ids = mapping.geneview_map_idx[gene_name];

  // The data-matrix corresponding to the SNPs (see comments below for update).
  Eigen::MatrixXd tmp_matrix =  \
      data.matrix.block(snp_ids.front(), 0, snp_ids.size(), data.n_samples); 

  // Find the intervals.
  interval_supports gene_supports = find_intervals(tmp_matrix);

  return gene_supports;
}


/*
  Enumerates all intervals in matrix mat, assuming that rows are ordered 
  according to genetic location of SNPs.
  Starts by first enumerating all intervals of length 1, and continues to 
  creater longer intervals by combination of intervals in previous length-
  layers.
*/
interval_supports EdgeEpistasisWY::find_intervals(const Eigen::MatrixXd& mat)
{
  // Init the containers to store the data.
  interval_supports intervals;
  interval_queue testable_queue;

  // Find the testable intervals in the first length-layer.
  intervals_layer_1(mat, intervals, testable_queue);

  // Find longer intervals.
  intervals_layer_x(intervals, testable_queue);

  // return the testable intervals only.
  return intervals;
}


/*
  Initialize all intervals of length 1 that are not prunable.
*/
void EdgeEpistasisWY::intervals_layer_1(const Eigen::MatrixXd& mat, 
                                        interval_supports& intervals, 
                                        interval_queue& testable_queue)
{

  // Add all intervals of length 1, if they are testable.
  for (int i=0; i<mat.rows(); i++)
  {

    // Compute the per-table support and corresponding minimum p-value.
    Eigen::VectorXd pt_support = \
        tarone.compute_per_table_support(mat.row(i));
    double min_pv = tarone.compute_minpval(pt_support);

    // Here, we do not process the individual intervals, since we only want
    // to account for pairwise interactions.
    if (!tarone.is_prunable(pt_support)) 
    {
      std::string key = std::to_string(i) + "_" + std::to_string(1);
      intervals[key] = mat.row(i);
      testable_queue.push_back(std::make_tuple(i, 1));
    }
  }

}


/*
  Function to find intervals of length>1 by combining lower-order intervals.
*/
void EdgeEpistasisWY::intervals_layer_x(interval_supports& intervals, 
                                        interval_queue& testable_queue)
{

  int length_layer = 2;

  while (!testable_queue.empty())
  {
    // get the next element.
    std::tuple<int, int> tmp_interval = testable_queue.front();
    testable_queue.pop_front();
    int tmp_start = std::get<0>(tmp_interval);
    int tmp_len = std::get<1>(tmp_interval);

    if (length_layer == tmp_len)
    {
      length_layer += 1;
    }

    // stop if the maximal length layer has been explored.
    if (max_length > 0 && length_layer > max_length)
    {
      break;
    }

    // grow the interval by combining the two previous ones.
    std::string prev_key = std::to_string(tmp_start-1) + "_" \
                         + std::to_string(tmp_len);
    std::string curr_key = std::to_string(tmp_start) + "_" \
                         + std::to_string(tmp_len);                

    // If the previous interval is in intervals, we can create the combination
    // (this is not the case if, for example, the previous interval was
    // prunable, or it exceeds the boundaries of the gene).                         
    if (intervals.find(prev_key) != intervals.end())
    {
      // Compute the joint support of the two intervals.
      Eigen::MatrixXd support = binary_or(intervals[prev_key], 
                                          intervals[curr_key]);

      // Compute the per-table support and corresponding minimum p-value.
      Eigen::VectorXd pt_support = tarone.compute_per_table_support(support);
      double min_pv = tarone.compute_minpval(pt_support);

      // Here, we do not process the individual intervals, since we only want
      // to account for pairwise interactions.
      if (!tarone.is_prunable(pt_support)) 
      {
        // create the key and store it in the intervals, as well as the queue.
        std::string key = std::to_string(tmp_start-1) + "_" \
                        + std::to_string(length_layer);
        intervals[key] = support;
        testable_queue.push_back(std::make_tuple(tmp_start-1, length_layer));
      }
    }
  }
}


/*
  Tests all combinations of intervals between gene0 and gene1.
*/
void EdgeEpistasisWY::test_interval_combinations(interval_supports gene_0_itvl,
                                                 interval_supports gene_1_itvl,
                                                 std::string gene_0_name,
                                                 std::string gene_1_name,
                                                 std::ofstream& out_file)
{

  for ( auto it0 = gene_0_itvl.begin(); it0 != gene_0_itvl.end(); ++it0 )
  {
    // Get the intervals in gene0.
    std::string key0 = it0->first;
    Eigen::MatrixXd supp0 = it0->second;

    for (auto it1=gene_1_itvl.begin(); it1!=gene_1_itvl.end(); ++it1)
    {
      // Get the intervals in gene1.
      std::string key1 = it1->first;
      Eigen::MatrixXd supp1 = it1->second;

      // Compute the binary or between the two interval-supports.
      Eigen::MatrixXd support = binary_or(supp0, supp1);

      // Compute the per-table support and corresponding minimum p-value.
      Eigen::VectorXd pt_support = tarone.compute_per_table_support(support);
      double min_pv = tarone.compute_minpval(pt_support);

      // Here we do adapt the threshold.
      if (tarone.is_testable(min_pv))
      {
        tarone.process_testable(pt_support, support);
      }
      if (!tarone.is_prunable(pt_support)) 
      {
        // Compute the p-value.
        int a = tarone.compute_supported_cases(support, data.labels);
        double pvalue = tarone.compute_pval(a, pt_support);

        // Write the interval combination to the output directory.
        out_file << pvalue << "," << min_pv << "," << gene_0_name << ",";
        out_file << key0 << "," << gene_1_name << "," << key1 << std::endl;
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