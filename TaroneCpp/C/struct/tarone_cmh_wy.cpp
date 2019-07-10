/*
* @Author: Anja Gumpinger
* @Date:   2018-11-13 19:32:54
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-10 21:25:11
*/

#ifndef _tarone_cmh_wy_cpp_
#define _tarone_cmh_wy_cpp_

#include <vector>
#include <ctime> // needed to generate random seed.
#include <unordered_map> 

#include <Eigen/Dense>

#include "../utils/EasyGWAS/chi2.h"
#include "../utils/utils.cpp"
#include "./tarone_cmh.cpp"


/*
  Main TaroneCMH.
*/
class TaroneCMHwy: public TaroneCMH
{

  public:
    /*
      Inherits all public memebers of TaroneCMH, i.e.:
      * double target_fwer;
      * int n_cov;
      * Eigen::VectorXd pt_samples;
      * Eigen::VectorXd pt_cases;
      * TaroneCMH (double, Eigen::VectorXd, Eigen::VectorXd);
      * double delta_t();
      * double corr_threshold();
      * long long n_testable();
      * double compute_minpval(Eigen::VectorXd x);
      * bool is_testable(double min_pv);
      * void process_testable(double min_pv);
      * bool is_prunable(Eigen::VectorXd x);
      * double compute_pval(int a, Eigen::VectorXd x);
      * void write_summary(std::string filename, long long n_enumerated, 
                           long long n_significant);
    */

    // The number of permutations.
    int n_perm;

    // The number of samples.
    int n_samples;

    // permutation labels, each row corresponds to the labels of one 
    // permutation.
    Eigen::MatrixXd perm_labels;
    
    // permutation p-values.
    std::vector<double> perm_pvalues;

    // Constructor.
    TaroneCMHwy() = default;
    TaroneCMHwy (double alpha, Eigen::VectorXd per_table_samples, 
                 Eigen::VectorXd per_table_cases, int n_permutations);

    // Function to compute the empirical FWER.
    double empirical_fwer();

    // Function to process the testables.
    void process_testable(Eigen::VectorXd pt_support,
                          const Eigen::MatrixXd& support);

    // Function to compute the corrected threshold.
    double corr_threshold();

    // Function to write summary file.
    void write_summary(std::string filename, int n_significant);

  protected:
    /*
      Inherits all protected members of TaroneCMH, i.e.:
      * int idx_t; 
      * std::vector<int> freq_counts;
      * double NGRID; 
      * double LOG10_MIN_PVAL;
      * double log10_p_step;
      * std::vector<double> gammat;
      * std::vector<double> gammabint;
      * std::vector<int> hypercorner_bnd;
      * std::vector<double> pgrid;
      * void init_pvalue_grid();
      * void init_cmh_helpers();
      * void init_freq_count();
      * int bucket_idx(double pval);
      * void decrease_threshold();
    */

  private:
    // Private members are accessible by the class itself.

    // Initializers.
    void init_permutation_pvalues();

    // Generation of label permutations.
    std::vector<int> generate_binary_labels_vector(int n_samples, int n_cases);
    void permute_vector(std::vector<int>& vec, long long seed);
    Eigen::MatrixXd permute_labels();

    // Functions to process a testable itemset.
    std::tuple<int, int> find_amin_amax(Eigen::VectorXd pt_support);
    void permutation_testing(Eigen::VectorXd pt_support, 
                             const Eigen::MatrixXd& support);
    void process_testable_spec(Eigen::VectorXd pt_support,
                               const Eigen::MatrixXd& support);
    


};


/*
  Constructor, initializes TaroneCMH and adds additional functionality.
*/
TaroneCMHwy::TaroneCMHwy(double alpha, Eigen::VectorXd per_table_samples, 
                         Eigen::VectorXd per_table_cases, int n_permutations) :
    TaroneCMH(alpha, per_table_samples, per_table_cases)
{

  // Number of permutations.
  n_perm = n_permutations;

  // Number of samples.
  n_samples = per_table_samples.sum();

  // Initialize the permutation p-values.
  init_permutation_pvalues();

  // Generate the label permutations.
  perm_labels = permute_labels();
}


/*
  Initialize a vector of 1's of length n_perm.
*/
void TaroneCMHwy::init_permutation_pvalues()
{
  for (int i=0; i<n_perm; i++)
  {
    perm_pvalues.push_back(1);
  }
}


/*
  Creates a binary vector with n_cases 1's and (n_samples-n_cases) 0's.
*/
std::vector<int> TaroneCMHwy::generate_binary_labels_vector(int n_samples, 
                                                            int n_cases)
{
  std::vector<int> vec;
  vec.insert(vec.end(), n_cases, 1);
  vec.insert(vec.end(), n_samples-n_cases, 0);
  return vec;
}


/*
  Randomly permute a binary vector.
*/
void TaroneCMHwy::permute_vector(std::vector<int>& vec, long long seed)
{
  srand(seed);
  std::random_shuffle(vec.begin(), vec.end());
}


/*
  Generates the label-permutations.
*/
Eigen::MatrixXd TaroneCMHwy::permute_labels()
{
  Eigen::MatrixXd permutations(n_perm, n_samples);

  // init the vector to store the current permuted labels.
  std::vector<int> perm_p_labels;

  for (int p=0; p<n_perm; p++)
  {
    // create the random seed for the current permutation.
    long long tmp_seed = time(NULL) * p;

    // create ONE randomized labels-vector PER covariate class and stack them 
    // together.
    for (int i=0; i<n_cov; i++)
    {
      // Generate the binary vector.
      std::vector<int> tmp_vec = \
          generate_binary_labels_vector(pt_samples[i], pt_cases[i]);

      // Shuffle the labels.
      permute_vector(tmp_vec, tmp_seed+i);

      // Add the permutations for the current covariate class to the permuted 
      // labels vector.
      perm_p_labels.insert(perm_p_labels.end(), tmp_vec.begin(), 
                             tmp_vec.end());
    }

    // Add the permuted label vector to the permutations matrix.
    for (int i=0; i<perm_p_labels.size(); i++)
    {
      permutations(p, i) = perm_p_labels[i];
    }
    perm_p_labels.clear();
  }

  return permutations;
}


/*
  Compute the empirical FWER.
*/
double TaroneCMHwy::empirical_fwer(){
  int count = 0;
  for (int p=0; p<n_perm; p++)
  {
    if (perm_pvalues[p] <= delta_t())
    {
      count += 1;
    } 
  }
  return (double)count / (double)n_perm;
}


/*
  Computes the corrected significance threshold as the alpha-quantile of 
  the permutation vectors.
*/
double TaroneCMHwy::corr_threshold(){
  return alpha_quantile(perm_pvalues, target_fwer);
}


/*  
  Wrapper calling the processing of a testable itemset. 
*/
void TaroneCMHwy::process_testable(Eigen::VectorXd pt_support, 
                                   const Eigen::MatrixXd& support)
{
  process_testable_spec(pt_support, support);
}


/*
  Process a testable itemset.
*/
void TaroneCMHwy::process_testable_spec(Eigen::VectorXd pt_support, 
                                        const Eigen::MatrixXd& support)
{
  // Do the permutation testing, including updating the perm_pvalues.
  permutation_testing(pt_support, support);

  // adapt the threshold.
  while(empirical_fwer() > target_fwer) 
  {
    decrease_threshold(); 
  }
}


/*
  For a itemset, compute the minimum and maximum value of 'a' that can be 
  obtained.
*/
std::tuple<int, int> TaroneCMHwy::find_amin_amax(Eigen::VectorXd pt_support)
{
  int a_min = 0;
  int a_max = 0;
  for (int i=0; i<n_cov; i++)
  {
    a_min += std::max((int)0, 
      (int)pt_support[i] - (int)pt_samples[i] +  (int)pt_cases[i]);
    a_max += std::min((int)pt_support[i], (int)pt_cases[i]);
  }

  return std::make_tuple(a_min, a_max);

}


/*
  Updates the vector of minimum attainable p-values given the current pattern.
  This is done in 2 steps:
  1. All p-values that could be obtained for all permutations can be precomputed
     by (i) finding all possible table configurations and (ii) precomputing
     their p-values.
  2. Updating the p-values for all permutations.
*/
void TaroneCMHwy::permutation_testing(Eigen::VectorXd pt_support, 
                                      const Eigen::MatrixXd& support)
{

  // Step 1(i): Find the values of a_min and a_max for the current itemset.
  std::tuple<int, int> a_min_max = find_amin_amax(pt_support);
  int a_min = std::get<0>(a_min_max);
  int a_max = std::get<1>(a_min_max);

  // Step 1(ii): Compute all possible p-values and store them in a vector.
  // This step can be parallelized.
  std::vector<double> tmp_pvals(a_max-a_min+1);
  #pragma omp parallel
  #pragma omp for 
  for (long long a=a_min; a<=a_max; a++)
  {
    int a_idx = a - a_min;
    tmp_pvals[a_idx] = compute_pval(a, pt_support);
  } 

  // Step 2: Compute the minimum p-value over all permutations.
  // This step can be parallelized.
  #pragma omp parallel
  #pragma omp for
  for (int p=0; p<n_perm; p++)
  {
    long long a = compute_supported_cases(support, perm_labels.row(p));
    int a_idx = a - a_min;
    double tmp = std::min(perm_pvalues[p], tmp_pvals[a_idx]);
    perm_pvalues[p] = tmp;
  }
}


void TaroneCMHwy::write_summary(std::string filename, int n_significant)
{

  std::ofstream file(filename);

  file << "Data set characteristics: "                      << std::endl;
  file << "Number samples: "            << pt_samples.sum() << std::endl;
  file << "Number cases: "              << pt_cases.sum()   << std::endl;
  file << "Number covariates: "         << n_cov            << std::endl; 
  file << std::endl;

  file << "Covariate classes:" << std::endl;
  for (int i=0; i<n_cov; i++)
  {
    file << "covar_" << i << ": ";
    file <<  pt_cases(i) << "/" << pt_samples(i);
    file << " (n_cases/n_samples)" << std::endl;
  }
  file << std::endl;

  // Tarone results.
  file << "Tarone results:"                                 << std::endl;
  file << "Number permutations: "     << n_perm             << std::endl;
  file << "Empirical FWER: "          << empirical_fwer()   << std::endl;
  file << "Testability threshold: "   << delta_t()          << std::endl;
  file << "Target fwer: "             << target_fwer        << std::endl;
  file << "Significance threshold: "  << corr_threshold()   << std::endl;
  file << "Number significant: "      << n_significant      << std::endl;

}






#endif
