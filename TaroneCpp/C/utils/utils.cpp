/*
* @Author: guanja
* @Date:   2019-07-04 17:39:39
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-10 20:36:26
*/

#ifndef _utils_cpp_
#define _utils_cpp_

#include <unordered_map> 
#include <vector>
#include <numeric> // iota

/* declarations. */

// print vector to std::cout.
template <typename T>
void print_vector(std::vector<T> vec);

// print nested vector to std::cout.
template <typename T>
void print_nested_vector(std::vector<std::vector<T>> nvec);

// function to print a map.
template <typename T>
void print_map(std::unordered_map<std::string, T> map);

// function to print a map.
template <typename T>
void print_map_vector(std::unordered_map<std::string, std::vector<T>> map);

// function to check whether a key is present in an unordered map.
bool check_key(std::unordered_map<std::string, int> m, std::string k);

// Argsort function.
template <class T>
std::vector<int> argsort(const std::vector<T> & unsorted);

// Transforms a std::vector to an Eigen::Vector;
template <class T>
Eigen::VectorXd vec_to_eigen(const std::vector<T>&  x);

// Transforms an Eigen::Vector to a std::vector;
std::vector<int> eigen_to_vec(const Eigen::VectorXd& x);

// Compute the binary or between two vectors (matrices).
Eigen::MatrixXd binary_or(Eigen::MatrixXd vec0, Eigen::MatrixXd vec1);

// Check if a file exists.
void check_file(std::string filename);

// Check if output file can be generated.
void check_out_file(std::string& filename);


/* Functions. */

template <typename T>
void print_vector(std::vector<T> vec){
  for(auto i: vec) 
  {
    std::cout << i << " ";
  }
  std::cout << std::endl;
}


template <typename T>
void print_nested_vector(std::vector<std::vector<T>> nvec){
  for(auto i: nvec) 
    {
      print_vector(i);
    }
  std::cout << std::endl;
}


template <typename T>
void print_map(std::unordered_map<std::string, T> map){
  for (std::pair<std::string, T> element : map)
  {
    std::cout << element.first << " :: " << element.second << std::endl;
  }
}


// function to print a map.
template <typename T>
void print_map_vector(std::unordered_map<std::string, std::vector<T>> map)
{
  for (std::pair<std::string, std::vector<T>> element : map)
  {
    std::cout << element.first << " :: ";
    for(auto i: element.second) 
    {
      std::cout << i << " ";
    }
    std::cout << std::endl;
  }
}



bool check_key(std::unordered_map<std::string, int> m, std::string k){
  if (m.find(k) == m.end()) 
  {
    return false; 
  }
  return true;
}


/*
  Argsort.
*/
template <class T>
std::vector<int> argsort(const std::vector<T> & unsorted){

  // Generate index vector.
  std::vector<int> index_sorted(unsorted.size());
  std::iota(index_sorted.begin(), index_sorted.end(), 0);
  // Sort the index vector, using unsorted for comparison
  std::sort(index_sorted.begin(), index_sorted.end(), 
    [&](int i1, int i2) { return unsorted[i1] < unsorted[i2]; });

  return index_sorted;
}


/*
  Transform an Eigen::VectorXd to a std::vector.
*/

std::vector<int> eigen_to_vec(const Eigen::VectorXd& x)
{
  std::vector<int> y;
  for(int i=0; i<x.size(); i++) 
  {
    y.push_back(x[i]);
  }
  return y;
}


/*
  Transform an std::vector to an Eigen::VectorXd.
*/
template <class T>
Eigen::VectorXd vec_to_eigen(const std::vector<T>&  x)
{
  Eigen::VectorXd y(x.size());
  for(int i=0; i<x.size(); i++) 
    {
      y(i) = x[i];
    }
  return y;
}


/*
 Compute the binary or between two Matrices of dimension (1, n_samples).
*/
Eigen::MatrixXd binary_or(Eigen::MatrixXd vec0, Eigen::MatrixXd vec1)
{
  // Add the two matrices.
  Eigen::MatrixXd sum_vec0_vec1 = vec0 + vec1;

  Eigen::MatrixXd vec_binary = Eigen::MatrixXd::Zero(vec0.rows(), vec0.cols());

  for (int i=0; i<vec_binary.cols(); i++)
  {
    if (sum_vec0_vec1(0, i) > 0)
    {
      vec_binary(0,i) = 1;
    }
  }
  return vec_binary;
}


/*
  Check if an input-file exists. Throw an error if not.
*/

void check_file(std::string filename){
  std::ifstream f(filename);
  if (!f.is_open()){
    std::cerr << "Error @ check_file: Unable to open file " << filename;
    std::cerr << std::endl;
    exit(-1);
  }
}


/*
  Checks if an output file can be generated.
*/
void check_out_file(std::string& filename){
  std::ofstream file_obj;
  file_obj.open(filename);
  if(!file_obj.is_open()) {
    std::cerr << "Error @ check_out_file: Unable to create output file ";
    std::cerr << filename << std::endl;
    exit(-1);
  }
}


/*
  Computes the alpha-quantile of a vector of doubles.
*/
double alpha_quantile(std::vector<double> vec, double alpha){

  // sort vector in increasing order.
  std::sort(vec.begin(), vec.end()); 

  // find the alpha-quantile and adjust for 0-based indexing.
  int idx = (int) (alpha * vec.size());
  idx -= 1;

  // if the value at index idx+1 is the same as the value at index idx, pick 
  // next smaller one (otherwise, quantile will be incorrect).
  if (vec[idx+1] != vec[idx])
  {
    return vec[idx];
  }else{
    int tmp_idx = idx;
    while(vec[tmp_idx] == vec[idx]) 
    {
      tmp_idx -= 1;
    }
    return vec[tmp_idx];
  }
}


#endif
