/*
* @Author: guanja
* @Date:   2019-07-04 17:39:39
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-04 18:30:05
*/

#ifndef _utils_cpp_
#define _utils_cpp_

#include <unordered_map> 
#include <vector>

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

// function to check whether a key is present in an unordered map.
bool check_key(std::unordered_map<std::string, int> m, std::string k);



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


bool check_key(std::unordered_map<std::string, int> m, std::string k){
  if (m.find(k) == m.end()) 
  {
    return false; 
  }
  return true;
}

#endif
