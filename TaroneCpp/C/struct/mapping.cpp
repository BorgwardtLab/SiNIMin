/*
* @Author: guanja
* @Date:   2019-07-04 18:40:57
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-09 18:15:24
*/

/* Mapping class.

Reads the mapping of SNPs to genes from a file.

The following pre-processing steps are done:
1. The edges are read, while self-loops are skipped.
2. Each edge is stored only one, i.e. in a non-symmetric fashion.
3. The vertices are mapped to integer IDs (for faster processing later), this
   mapping is stored in to unordered_maps, one containing the mapping from 
   string name to integer name, the second one from integer to string name.
4. The string-named edges are translated to integer-edges.

*/

#ifndef _maping_cpp_
#define _maping_cpp_

// sstringstream.
#include <sstream>

#include <algorithm>

#include <vector>
#include <tuple>
#include <unordered_map> 

#include "../utils/utils.cpp"


class Mapping
{

public:

  // Number of edges.
  long long n_genes;

  // mapping of char label to int label.
  std::unordered_map<std::string, std::vector<std::string>> geneview_map_str;

  // mapping of char label to int label.
  std::unordered_map<std::string, std::vector<int>> geneview_map_idx;


  // Empty constructor for initialization.
  Mapping() = default;
  // Constructor with edge file.
  Mapping(std::string);

  // Function to translate the SNP-IDs to their position in the data matrix.
  void translate_snp_ids(std::vector<std::string> snp_ordering);

private:

  void read_mapping(std::string map_file);

};


Mapping::Mapping(std::string map_file)
{

  check_file(map_file);
  read_mapping(map_file);
  
}


/*
  Reads the mapping from a file.
*/
void Mapping::read_mapping(std::string map_file)
{

  std::string line;
  std::stringstream ss_line;

  int col_idx;
  std::string value;
  std::string gene_name;
  std::vector<std::string> snp_vec;

  std::ifstream file(map_file);
  while(std::getline(file, line))
  {
    col_idx = 0;
    ss_line << line;
    while(ss_line >> value) 
    {
      if (col_idx==0)
      {
        gene_name = value;
      }else{
        snp_vec.push_back(value);
      }
      geneview_map_str[gene_name] = snp_vec;
      col_idx += 1;
    }
    snp_vec.clear();
    ss_line.clear();
  }
  file.close();   
}


/*
  Translates the SNP-IDs to their position given in the snp_ordering file.
*/
void Mapping::translate_snp_ids(std::vector<std::string> snp_ordering)
{

  // create the mapping of SNP <-> index.
  int position = 0;
  std::unordered_map<std::string, int> position_dict;
  for (int i=0; i<snp_ordering.size(); i++)
  {
    position_dict[snp_ordering[i]] = position;
    position += 1;
  }

  // translate the snp_names to row_idx in data matrix.
  for (std::pair<std::string, std::vector<std::string>> elm : geneview_map_str)
  {
    std::vector<int> snp_idx;
    for (int j=0; j<elm.second.size(); j++)
    {
      snp_idx.push_back(position_dict[elm.second[j]]);
    }
    geneview_map_idx[elm.first] = snp_idx;
    snp_idx.clear();
  } 
}


#endif