/*
* @Author: guanja
* @Date:   2019-07-04 17:19:08
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-08 13:18:56
*/

// Include standard libs.
#include <iostream>
#include <string>

// Include C++ data classes.
#include <Eigen/Dense>
#include "../../C/struct/data.cpp"


int main(int argc, char** argv)
{
  std::string data_file = "../../examples/data/data.txt";
  std::string labels_file = "../../examples/data/labels.txt";
  std::string covar_file = "../../examples/data/covar.txt";

  Data data(data_file, labels_file, covar_file);

  std::cout << data.matrix << std::endl;
  std::cout << data.labels << std::endl;
  std::cout << data.pt_samples << std::endl;
  std::cout << data.pt_cases << std::endl;
  std::cout << data.pt_samples_cumsum << std::endl;


}