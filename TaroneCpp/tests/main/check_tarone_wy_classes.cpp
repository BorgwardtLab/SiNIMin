/*
* @Author: guanja
* @Date:   2019-07-04 17:19:08
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-08 18:41:40
*/

// Include standard libs.
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>

// Include C++ data classes.
#include <vector>
#include <set>
#include <Eigen/Dense>


#include "../../C/struct/tarone_cmh_wy.cpp"
#include "../../C/utils/utils.cpp"



int main(int argc, char** argv)
{

  // set the target FWER.
  double alpha = 0.05;
  int n_perm = 1000;


  /* ***************************************************************************
    Check the computation of the support/ pvalues etc.
  *************************************************************************** */

  // Init data vecs to test the functions.
  std::vector<int> samples_pt_vec = {6, 4};
  std::vector<int> cases_pt_vec = {3, 2};
  std::vector<int> label_vec = {1, 1, 1, 0, 0, 0, 1, 1, 0, 0};
  std::vector<int> support_vec = {1, 0, 1, 0, 1, 0, 0, 1, 0, 1};
  

  // Turn them into Eigen::VectorXd.
  Eigen::VectorXd samples_pt = vec_to_eigen(samples_pt_vec);
  Eigen::VectorXd cases_pt = vec_to_eigen(cases_pt_vec);
  Eigen::VectorXd labels = vec_to_eigen(label_vec);
  Eigen::MatrixXd support(1, 10);
  for (int i=0; i<support_vec.size(); i++)
  {
    support(0, i) = support_vec[i];
  }


  // Test 1: Check if the support_pt is computed correctly.
  TaroneCMHwy tarone_test1(alpha, samples_pt, cases_pt, n_perm);

  Eigen::VectorXd support_pt = tarone_test1.compute_per_table_support(support);

  if (support_pt[0] != 3)
  {
    std::cerr << "support_pt incorrect. (" << support_pt[0] << ", " ;
    std::cerr << support_pt[1] << "). Exiting ." << std::endl;
    std::exit(-1);
  }
  if (support_pt[1] != 2)
  {
    std::cerr << "support_pt incorrect. (" << support_pt[0] << ", ";
    std::cerr << support_pt[1] << "). Exiting ." << std::endl;
    std::exit(-1);
  }

  // Test 2: Check if the supported case count a is computed correctly.
  int a = tarone_test1.compute_supported_cases(support, labels);
  if (a != 3)
  {
    std::cerr << "error @ supported cases 'a' incorrect (" << a << "). " ;
    std::cerr << "Exiting ." << std::endl;
    std::exit(-1);
  }

  // Test 3: Compute the minimum p-value. (exact value computed by hand for
  // the above example)
  double min_pv = tarone_test1.compute_minpval(support_pt);

  if (min_pv - 0.0015654 > 1e-8)
  {
    std::cerr << "error @ min-pval. (" << min_pv << "). Exiting ." << std::endl;
    std::exit(-1);
  }
  if ( 0.0015654 - min_pv > 1e-8)
  {
    std::cerr << "error @ min-pval. (" << min_pv << "). Exiting ." << std::endl;
    std::exit(-1);
  }


  // Test 4: Compute the p-value.
  a = tarone_test1.compute_supported_cases(support, labels);
  double pval = tarone_test1.compute_pval(a, support_pt);

  if (pval - 0.5270892568655381 > 1e-8)
  {
    std::cerr << "error @ pvalue. (" << pval << "). Exiting ." << std::endl;
    std::exit(-1);
  }
  if ( 0.5270892568655381 - pval > 1e-8)
  {
    std::cerr << "error @ pvalue. (" << pval << "). Exiting ." << std::endl;
    std::exit(-1);
  }


  /* ***************************************************************************
    Check the computation of permutation p-values.
  *************************************************************************** */

  // Given the minimum attainable p-value from before, we now look into
  // the permutation based lowering of the threshold.

  tarone_test1.process_testable(support_pt, support);
  std::cout << "delta_t = " << tarone_test1.delta_t() << std::endl;
  std::cout << "emp. FWER = " << tarone_test1.empirical_fwer() << std::endl;

  // try next pattern.
  support_vec = {1, 1, 0, 0, 0, 0, 1, 1, 0, 0};
  for (int i=0; i<support_vec.size(); i++)
  {
    support(0, i) = support_vec[i];
  }
  support_pt = tarone_test1.compute_per_table_support(support);
  tarone_test1.process_testable(support_pt, support);
  std::cout << "delta_t = " << tarone_test1.delta_t() << std::endl;
  std::cout << "emp. FWER = " << tarone_test1.empirical_fwer() << std::endl;

  // try next pattern.
  support_vec = {0, 1, 0, 0, 0, 0, 1, 1, 0, 0};
  for (int i=0; i<support_vec.size(); i++)
  {
    support(0, i) = support_vec[i];
  }
  support_pt = tarone_test1.compute_per_table_support(support);
  tarone_test1.process_testable(support_pt, support);
  std::cout << "delta_t = " << tarone_test1.delta_t() << std::endl;
  std::cout << "emp. FWER = " << tarone_test1.empirical_fwer() << std::endl;

}