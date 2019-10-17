/*
* @Author: guanja
* @Date:   2019-07-10 15:28:17
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-10 16:09:29
*/

#ifndef _output_cpp_
#define _output_cpp_

#include <sstream>
#include <string>
#include <fstream>

// filter significant files from a pvalues.cvs file.
int write_significant_edge_epi(double threshold, 
                               std::string input_filename, 
                               std::string output_filename);



/*
  Reads the significant itemsets from the pvalue.csv file of a Tarone run.
*/
int write_significant_edge_epi(double threshold, 
                               std::string input_filename, 
                               std::string output_filename)
{

  std::ifstream in_file(input_filename);
  std::ofstream out_file(output_filename);

  out_file << "p-value,min_pval,gene0,start0_len0,gene1,start1_len1" ;
  out_file << std::endl;

  std::string line;
  std::string cell;
  std::getline(in_file,line);
  long long sig_count = 0;
  
  while(std::getline(in_file, line))
  {
    std::istringstream ss_line(line);
    while(std::getline(ss_line, cell, ',')) 
    {
      if (std::stod(cell) < threshold)
      {
        sig_count ++;
        out_file << line << std::endl;
      }
      break;
    }
    ss_line.clear();
  }
  return sig_count;
}

#endif