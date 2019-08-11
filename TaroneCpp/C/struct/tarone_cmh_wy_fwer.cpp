/*
* @Author: Anja Gumpinger
* @Date:   2018-11-13 19:32:54
* @Last Modified by:   guanja
* @Last Modified time: 2019-08-11 16:50:32
*/

#ifndef _tarone_cmh_wy_fwer_cpp_
#define _tarone_cmh_wy_fwer_cpp_

#include <vector>
#include <ctime> // needed to generate random seed.
#include <unordered_map> 

#include <Eigen/Dense>

#include "../utils/EasyGWAS/chi2.h"
#include "../utils/utils.cpp"
#include "./tarone_cmh.cpp"
#include "./tarone_cmh_wy.cpp"




/*
  Main TaroneCMH.
*/
class TaroneCMHwyFWER: public TaroneCMHwy
{
  public:
    /*
      Inherits all public memebers of TaroneCMH, i.e.:
    */

    // 
    double threshold;

    // Constructor.
    TaroneCMHwyFWER() = default;
    TaroneCMHwyFWER (double thresh, Eigen::VectorXd per_table_samples, 
                     Eigen::VectorXd per_table_cases, int n_permutations);

    // Estimate the empirical FWER.
    double empirical_fwer();

    // Adapt the threshold.
    void permutation_testing(Eigen::VectorXd pt_support, 
                             const Eigen::MatrixXd& support); 


    // Function to write summary file.
    void write_summary(std::string filename);

  protected:
    /*
      Inherits all protected members of TaroneCMH, i.e.:
    */

    // Functions to process a testable itemset.
    std::tuple<int, int> find_amin_amax(Eigen::VectorXd pt_support);

};


/*
  Constructor, initializes TaroneCMH and adds additional functionality.
*/
TaroneCMHwyFWER::TaroneCMHwyFWER(double thresh, 
                                 Eigen::VectorXd per_table_samples, 
                                 Eigen::VectorXd per_table_cases, 
                                 int n_permutations) :
    TaroneCMHwy(0.05, per_table_samples, per_table_cases, n_permutations)
{

  // Set the threshold.
  threshold = thresh;

  // All other initializations are done when initializing parent class.

}



/*
  Compute the empirical FWER.
*/
double TaroneCMHwyFWER::empirical_fwer(){
  int count = 0;
  for (int p=0; p<n_perm; p++)
  {
    if (perm_pvalues[p] <= threshold)
    {
      count += 1;
    } 
  }
  return (double)count / (double)n_perm;
}



void TaroneCMHwyFWER::write_summary(std::string filename)
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
  file << "Tarone results:"                                       << std::endl;
  file << "Number permutations: "           << n_perm             << std::endl;
  file << "threshold: "                     << threshold          << std::endl;
  file << "Empirical FWER at threshold: "   << empirical_fwer()   << std::endl;
}



/* -----------------------------------------------------------------------------
  ALL OF THE SUBSEQUENT FUNCTIONS SHOULD MOVE TO PROTECTED IN WY-TARONE!!!
  Then we can inherit them.
*/


/*
  For a itemset, compute the minimum and maximum value of 'a' that can be 
  obtained.
*/
std::tuple<int, int> TaroneCMHwyFWER::find_amin_amax(Eigen::VectorXd pt_support)
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
void TaroneCMHwyFWER::permutation_testing(Eigen::VectorXd pt_support, 
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






#endif
