/*
* @Author: Anja Gumpinger
* @Date:   2018-11-13 19:32:54
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-11 12:58:12
*/

#ifndef _tarone_cmh_cpp_
#define _tarone_cmh_cpp_

#include <vector>

#include <Eigen/Dense>

#include "../utils/EasyGWAS/chi2.h"
#include "../utils/utils.cpp"


/*
  Main TaroneCMH.
*/
class TaroneCMH{

  public:
    // Public members are accessible by anyone calling the function.

    // target family-wise error rate.
    double target_fwer;

    // Number of covariates.
    int n_cov;

    // Vectors containing the number of samples/cases per covariate.
    Eigen::VectorXd pt_samples;
    Eigen::VectorXd pt_cases;
      
    // Default constructor.
    TaroneCMH() = default;

    // Constructor with data input.
    TaroneCMH (double alpha, Eigen::VectorXd per_table_samples, 
               Eigen::VectorXd per_table_cases);
      
    // Functions to compute the current thresholds.
    double delta_t();
    double corr_threshold();
    long long n_testable();
    
    // Functions to assess and process testability of a itemset.
    double compute_minpval(Eigen::VectorXd pt_support);
    bool is_testable(double min_pv);
    void process_testable(double min_pv);
    bool is_prunable(Eigen::VectorXd x);
    double compute_supported_cases(const Eigen::MatrixXd& support, 
                                   const Eigen::VectorXd& labels);
    Eigen::VectorXd compute_per_table_support(Eigen::MatrixXd support);
    double compute_pval(int a, Eigen::VectorXd pt_support);

    // Function to write a summary of the tarone run.
    void write_summary(std::string filename, int n_significant);

  protected:
    // Protected members are accessible by the class itself and any class 
    // inheriting.

    // The index pointing to the current threshold.
    int idx_t; 

    // vector to keep track of the frequency counts.
    std::vector<int> freq_counts;

    // vector containing the cumsum of the samples.
    Eigen::VectorXd pt_samples_cumsum;

    // Parameters to compute the p-value grid.
    double NGRID; 
    double LOG10_MIN_PVAL;
    double log10_p_step;

    // helpers to compute the CMH statistics.
    std::vector<double> gammat;
    std::vector<double> gammabint;
    std::vector<int> hypercorner_bnd;

    // The grid of p-values.
    std::vector<double> pgrid;

    void init_cumsum();
    void init_pvalue_grid();
    void init_cmh_helpers();
    void init_freq_count();
    int bucket_idx(double pval);
    void decrease_threshold();

  private:
    // Private members are accessible by the class itself.

    void process_testable_spec(double min_pv);

};


/*
  Constructor function. (Called at initialization of class, precomputes every-
  thing and sets all default values).
*/
TaroneCMH::TaroneCMH (double alpha, Eigen::VectorXd per_table_samples, 
                      Eigen::VectorXd per_table_cases) {

  // init the target_fwer.
  target_fwer = alpha;

  // set the number of covariates.
  n_cov = per_table_samples.rows();

  // set the samples/cases per covariate class.
  pt_samples = per_table_samples;
  pt_cases = per_table_cases;

  // init the cumsum of samples.
  init_cumsum();

  // init the p-value grid.
  init_pvalue_grid();

  // init the helpers to compute p-values.
  init_cmh_helpers();

  // init the vector counting how often a p-value has occurred.
  init_freq_count();

  // initialize the index pointing to the threshold.
  idx_t = 1; 
  
}


/*
  Returns the current testability threshold.
*/ 
double TaroneCMH::delta_t()
{
  return pgrid[idx_t];
}


/*
  Returns the corrected significance threshold.
*/ 
double TaroneCMH::corr_threshold()
{
  if (n_testable() == 0)
  {
    return target_fwer;
  }else{
    return target_fwer / n_testable();  
  }
}


/* 
  Computes the number of testable items.
*/
long long TaroneCMH::n_testable()
{
  long long count = 0;
  for (int i=idx_t+1; i<=NGRID; i++)
  {
    count += freq_counts[i];
  }
  return count;
}


/*
  Inits the vector that will keep track of the counts of itemsets that are
  testable at a given threshold.
*/
void TaroneCMH::init_freq_count()
{
  for (int i=0; i<=NGRID; i++) 
  {
    freq_counts.push_back(0);
  }
}


/*
  Precomputes the p-values that will be used as testability thresholds.
*/
void TaroneCMH::init_pvalue_grid()
{
  NGRID = 500;
  LOG10_MIN_PVAL = -30.0;
  log10_p_step =- LOG10_MIN_PVAL/NGRID;

  // init the p-value grid.
  double log10_p = 0;
  for(int j=0; j<=NGRID; j++) {
    log10_p -= log10_p_step;
    pgrid.push_back(pow(10, log10_p));
  }
}


/*
  Function to init the CMH helpers (to compute the p-values faster).
*/
void TaroneCMH::init_cmh_helpers()
{
  for (int i=0; i<n_cov; i++){
    long long tmp_c = pt_samples(i) - pt_cases(i);  // number of controls.
    
    gammat.push_back((double)pt_cases(i) / (double)pt_samples(i));
    
    gammabint.push_back(gammat[i]*(1-gammat[i])); 
    
    int maj = pt_cases(i) > tmp_c ? pt_cases(i) : tmp_c;
    hypercorner_bnd.push_back(maj);
  }
}


/*
  Create the cummulative sum of number of samples per covariate class.
*/
void TaroneCMH::init_cumsum()
{
  pt_samples_cumsum = Eigen::VectorXd::Zero(n_cov);
  pt_samples_cumsum(0) = pt_samples(0);
  for (int i=1; i<pt_samples_cumsum.rows(); i++)
  {
    pt_samples_cumsum(i) = pt_samples(i) + pt_samples_cumsum(i-1);
  }
}


/*
  Process a testable itemset.
*/
void TaroneCMH::process_testable(double min_pv)
{
  // Call the private function to do the method-specific processing of a 
  // testable itemset.
  process_testable_spec(min_pv);
}


/*
  Processes a testable itemset in the method-dependent fashion.
*/
void TaroneCMH::process_testable_spec(double min_pv)
{
  // Compute the bucket index.
  int idx = bucket_idx(min_pv);

  // Increase the count vector of testable items at the threshold.
  freq_counts[idx] += 1;

  // Check FWER criterion.
  while (n_testable() * delta_t() > target_fwer)
  {
    decrease_threshold();
  }
}



/*
  Function to decrease threshold and update the W-shape parameters.
*/
void TaroneCMH::decrease_threshold()
{
  // increase the pointer by one, thereby lowering the threshold.
  idx_t += 1; 
}


/*
  Testability of a support value.
  Compute the minimum p-value, and if it is smaller than delta_t, it is 
  testable.
*/
bool TaroneCMH::is_testable(double min_pv)
{
  if (min_pv <= delta_t()) 
    {
      return true;
    }
  return false;
}


/*
  given a p-value, find its 'bucket' (to count the number of testables).
  If the p-value is 0, the log10 is not defined, hence we automatically put 
  this p-value into the last bucket.
*/
int TaroneCMH::bucket_idx(double pval)
{
  if (pval == 0)
  {
    return NGRID;
  }

  int idx = (int)floor(-log10(pval)/log10_p_step);

  if(idx<0) 
    {
      idx = 0;
    }
  if(idx>NGRID) 
    {
      idx = NGRID;
    }
  return idx;
}


/*
  Compute the cell-count 'a' in the contingency table, i.e. those cases that
  have support = 1.
*/
double TaroneCMH::compute_supported_cases(const Eigen::MatrixXd& support, 
                                          const Eigen::VectorXd& labels)
{
  long long a = 0;
  for (int i=0; i<labels.rows(); i++) 
  {
    if (labels(i) && support(i)) 
      {
        a++;
      }
  }
  return a;
}


/*
  Compute the per-table support.
*/
Eigen::VectorXd TaroneCMH::compute_per_table_support(Eigen::MatrixXd support)
{
  int k = 0;
  Eigen::VectorXd pt_support = Eigen::VectorXd::Zero(n_cov);
  for (int i=0; i<pt_samples.sum(); i++)
  {
    if (support(0, i) == 1) 
    {
      pt_support(k) += 1;
    }
    if (i+1 == pt_samples_cumsum(k)) 
      {
        k += 1;
      }
  }
  return pt_support;
}

/*
  Compute the CMH-test p-value.
*/
double TaroneCMH::compute_pval(int a, Eigen::VectorXd x)
{
  double num = a;
  double den = 0;

  for (int i=0; i<x.size(); i++)
  {
    num -= x(i) * gammat[i];
    den += x(i) * (1 - ((double)x(i))/pt_samples(i)) * gammabint[i];
  }

  num *= num;
  if (den == 0) 
    {
      return 1;
    }else{
      return Chi2_sf(num/den, 1); // switch this for boost chi2
    }
}


/* 
  Compute the minimum p-value of the CMH test using only the table margins.
*/
double TaroneCMH::compute_minpval(Eigen::VectorXd x)
{
  // Holds the numerator of T_amin
  double left_tail_num = 0; 
  // Holds the numberator of T_amax
  double right_tail_num = 0; 
  // Holds the denominator.
  double den = 0;
  // Auxiliaries for computation of T-scores.
  double aux1, aux2;

  for(int k=0; k<x.rows(); k++)
  {
    aux1 = x(k) - (pt_samples(k) - pt_cases(k)); 
    aux2 = x(k) * gammat[k];
    left_tail_num += ((aux1 > 0) ? aux1 : 0) - aux2;
    right_tail_num += ((x(k) > pt_cases(k)) ? pt_cases(k) : x(k)) - aux2;
    den += x[k] * (1-((double)x(k))/pt_samples(k)) * gammabint[k];
  }
  left_tail_num *= left_tail_num; 
  right_tail_num *= right_tail_num;

  if(den==0)
  { 
    return 1;
  }else{ 
    double minpval = Chi2_sf(((left_tail_num > right_tail_num) ? 
      left_tail_num : right_tail_num)/den, 1);
    return minpval;
  }
}


/*
  Compute the pruning criterion based on the OR encoding, i.e. increasing x.
  Following supplementary Alg. 4 in Bioinformatics 2017 paper.
*/
bool TaroneCMH::is_prunable(Eigen::VectorXd x)
{

  // allocate memory.
  std::vector<double> f_vals;
  std::vector<double> g_vals;
  std::vector<double> beta;
  std::vector<int> idx_beta_sorted; 

  double f_sum = 0;
  double g_sum = 0;
  double Tcmh_aux_corner;
  double Tcmh_max_corner_r;
  double Tcmh_max_corner_l;


  // If for any of the tables, its margin x is smaller than the maximum of n 
  // and N-n, then we cannot prune the set (we are not in the "top-right" 
  // hypercorner)
  for(int k=0; k<x.rows(); k++) 
  {
    if(x[k] < hypercorner_bnd[k]) return false;
  }

  // compute the righthandside values.
  for (int k=0; k<x.rows(); k++)
  {
    if (x(k) < pt_samples(k))
    {
      double f = gammat[k] * (pt_samples(k)-x(k));
      double g = gammabint[k] * x(k) * (1-((double)x(k))/pt_samples(k));
      f_vals.push_back(f);
      g_vals.push_back(g);
      beta.push_back(g/f);
    }
  }

  // get permutation of indices of beta values.
  idx_beta_sorted = argsort(beta);

  // compute CMH test statistic.
  // Skip tables that only have one class, i.e. g_sum == 0.
  f_sum = 0;
  g_sum = 0;
  Tcmh_max_corner_r = 0;
  for(int k=0; k<beta.size(); k++)
  {
    f_sum += f_vals[idx_beta_sorted[k]];
    g_sum += g_vals[idx_beta_sorted[k]];
    if (g_sum == 0) 
      {
        continue;
      }
    Tcmh_aux_corner = (f_sum*f_sum)/g_sum;
    Tcmh_max_corner_r = (Tcmh_max_corner_r >= Tcmh_aux_corner) ? 
      Tcmh_max_corner_r : Tcmh_aux_corner; 
  }

  g_vals.clear();
  f_vals.clear();
  beta.clear();

  // compute the lefthandside values.
  for(int k=0; k<x.rows(); k++)
  {
    // Discard all dimensions for which x[k]==Nt[k], as they don't 
    // contribute to the function neither in the numerator nor the 
    // denominator.
    if(x(k) < pt_samples(k))
    {
      double f = (1-gammat[k]) * (pt_samples(k)-x(k));
      double g = gammabint[k] * x(k) * (1-((double)x(k))/pt_samples(k));
      f_vals.push_back(f);
      g_vals.push_back(g);  
      beta.push_back(g/f);
    }
  }

  // get permutation of indices of beta values.
  idx_beta_sorted = argsort(beta);

  f_sum = 0; 
  g_sum = 0; 
  Tcmh_max_corner_l = 0;
  for(int k=0; k<beta.size(); k++)
  {
    f_sum += f_vals[idx_beta_sorted[k]];
    g_sum += g_vals[idx_beta_sorted[k]];
    if (g_sum == 0)
      {
        continue;
      }
    Tcmh_aux_corner = (f_sum*f_sum)/g_sum;
    Tcmh_max_corner_l = (Tcmh_max_corner_l >= Tcmh_aux_corner) ? 
      Tcmh_max_corner_l : Tcmh_aux_corner; 
  }

  double min_pval = Chi2_sf(((Tcmh_max_corner_r >= Tcmh_max_corner_l) ? 
    Tcmh_max_corner_r : Tcmh_max_corner_l),1);

  return  min_pval > delta_t(); 
}


void TaroneCMH::write_summary(std::string filename, int n_significant)
{

  std::ofstream file(filename);

  file << "Data set characteristics: " << std::endl;
  file << "Number samples: " << pt_samples.sum() << std::endl;
  file << "Number cases: " << pt_cases.sum() << std::endl;
  file << "Number covariates: " << n_cov << std::endl << std::endl;

  file << "Covariate classes:" << std::endl;
  for (int i=0; i<n_cov; i++)
  {
    file << "covar_" << i << ": ";
    file <<  pt_cases(i) << "/" << pt_samples(i);
    file << " (n_cases/n_samples)" << std::endl;
  }
  file << std::endl;

  // Tarone results.
  file << "Tarone results:" << std::endl;
  file << "Number testable: " << n_testable() << std::endl;
  file << "Testability threshold: " << delta_t() << std::endl;
  file << "Target fwer: " << target_fwer << std::endl;
  file << "Significance threshold: " << corr_threshold() << std::endl;
  file << "Number significant: " << n_significant << std::endl;

}

#endif
