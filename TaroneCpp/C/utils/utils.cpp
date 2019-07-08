/*
* @Author: guanja
* @Date:   2019-07-04 17:39:39
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-08 15:32:52
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

#endif
