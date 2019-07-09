/*
* @Author: guanja
* @Date:   2019-07-04 16:56:06
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-09 18:15:29
*/


/* Edge class.

Reads the edges from an edge file, where each row contains the two adjacent 
vertices as string names.

The following pre-processing steps are done:
1. The edges are read, while self-loops are skipped.
2. Each edge is stored only one, i.e. in a non-symmetric fashion.
3. The vertices are mapped to integer IDs (for faster processing later), this
   mapping is stored in to unordered_maps, one containing the mapping from 
   string name to integer name, the second one from integer to string name.
4. The string-named edges are translated to integer-edges.

*/

#ifndef _edges_cpp_
#define _edges_cpp_

// sstringstream.
#include <sstream>

#include <algorithm>

#include <vector>
#include <tuple>
#include <unordered_map> 

#include "../utils/utils.cpp"


class Edges{
public:

  // Number of edges.
  long long n_edges;

  // Number of nodes.
  long long n_nodes;

  // representation of edges as tuple of strings.
  std::vector<std::vector<std::string>> edges_str;

  // representation of edges as tuple of ints.
  std::vector<std::vector<int>> edges_int;

  // mapping of char label to int label.
  std::unordered_map<std::string, int> vertex_str_to_int_map;

  // mapping of int label to string label.
  std::unordered_map<int, std::string> vertex_int_to_str_map;

  // Empty constructor for initialization.
  Edges() = default;
  // Constructor with edge file.
  Edges(std::string);

private:

  void read_edges(std::string edge_file);
  void remove_duplicate_edges();
  void create_str_int_maps();
  void translate_edges();

};


Edges::Edges(std::string edge_fn){

  check_file(edge_fn);

  read_edges(edge_fn);
  remove_duplicate_edges();
  
  // Create the mapping of integer and string labels.
  create_str_int_maps();

  // Translate the edges to the integer labels.
  translate_edges();

  // set the number of edges.
  n_edges = edges_str.size();
  n_nodes = vertex_str_to_int_map.size();
}


/*
  Reads the edges from a file.
*/
void Edges::read_edges(std::string edge_file){

  std::string line;
  std::stringstream ss_line;

  std::string vertex_name;
  std::vector<std::string> edge_vec;

  std::ifstream file(edge_file);
  while(std::getline(file, line)){
    ss_line << line;
    while(ss_line >> vertex_name) {
      edge_vec.push_back(vertex_name);
    }

    // sort the vertices in alphabetic order to later remove duplicates.
    std::sort(edge_vec.begin(), edge_vec.end());

    // If self loop, continue.
    if (!(edge_vec[0] == edge_vec[1]))
    {
      edges_str.push_back(edge_vec);
    }

    edge_vec.clear();
    ss_line.clear();
  }
  file.close();   
}


/* 
  Remove duplicate edges. 
*/
void Edges::remove_duplicate_edges(){
  // sort the nested vector.
  std::sort(edges_str.begin(), edges_str.end());
  // keep unique edges only.
  edges_str.erase(std::unique(edges_str.begin(), edges_str.end()), 
                  edges_str.end());
}


/* Create the edges from string to integer label and vice versa. */
void Edges::create_str_int_maps(){

  int current_count = 0;

  for(int i=0; i<edges_str.size(); i++)
  {
    for(int j=0; j<2; j++)
    {
      if (!check_key(vertex_str_to_int_map, edges_str[i][j])){
        vertex_str_to_int_map[edges_str[i][j]] = current_count;
        vertex_int_to_str_map[current_count] = edges_str[i][j];
        current_count += 1;
      }
    }
  }
}


/* 
  Translate vertices to int. 
*/
void Edges::translate_edges(){
  // sort the nested vector.

  std::vector<int> edge_int;

  for (int i=0; i<edges_str.size(); i++)
  {
    for (int j=0; j<2; j++)
    {
      edge_int.push_back(vertex_str_to_int_map[edges_str[i][j]]);
    }
    edges_int.push_back(edge_int);
    edge_int.clear();
  }
}


#endif