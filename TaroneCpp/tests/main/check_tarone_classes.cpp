/*
* @Author: guanja
* @Date:   2019-07-04 17:19:08
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-11 11:39:06
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


#include "../../C/struct/tarone_cmh.cpp"
#include "../../C/utils/utils.cpp"



int main(int argc, char** argv)
{

  // Set the pt_cases.
  Eigen::VectorXd cases = Eigen::VectorXd::Zero(2);
  cases(0)=10;
  cases(1)=5;

  Eigen::VectorXd samples = Eigen::VectorXd::Zero(2);
  samples(0) = 15;
  samples(1) = 20;
  // set the target FWER.
  double alpha = 0.05;

  // Init tarone.
  TaroneCMH tarone(alpha, samples, cases);

  std::vector<double> minpv_vec{ 0.001, 0.99, 0.15, 0.03, 0.02, 0.5, 0.75, 
                                 0.005, 0.0000069, 0.12, 0.002, 0.0045, 
                                 0.870964};

  for (int i=0; i<minpv_vec.size(); i++)
  {
    std::cout << "Processing min-pvalue = " << minpv_vec[i] << ". ";
    std::cout << "Thresh. before proc.: delta_t = " << tarone.delta_t() << ". ";
    tarone.process_testable(minpv_vec[i]);
    std::cout << "Thresh. after proc.: delta_t = " << tarone.delta_t() << ". ";
    std::cout << "n_testable = " << tarone.n_testable() << ". " << std::endl;
  }

  // Assert that a p-value of 0 is processed correctly, i.e. put into the 
  // last bin.

  std::cout << "++++++++++" << std::endl;
  std::cout << "Processing minpval=0" << std::endl;
  tarone.process_testable(0);
  std::cout << "n_testable = " << tarone.n_testable() << ". " << std::endl;

  std::cout << "++++++++++" << std::endl;
  std::cout << "Processing minpval=8.70964e-31" << std::endl;
  tarone.process_testable(8.70964e-31);
  std::cout << "n_testable = " << tarone.n_testable() << ". " << std::endl;

}