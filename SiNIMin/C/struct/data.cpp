/*
* @Author: Anja Gumpinger
* @Date:   2018-11-12 13:57:32
* @Last Modified by:   guanja
* @Last Modified time: 2019-07-10 21:38:36
*/

#ifndef _data_cpp_
#define _data_cpp_

#include <sstream>
#include <fstream>

#include <vector>
#include <tuple>
#include <set>

#include <Eigen/Dense>
#include "../utils/utils.cpp"


/*
  Sructure holding data.
*/
class Data{
public:

  // Metrics of the data set.
  long long n_features;
  long long n_samples;
  long long n_cases;
  long long n_covar;

  // The data set and the labels.
  Eigen::MatrixXd matrix;
  Eigen::VectorXd labels;

  // The snp-IDs.
  std::vector<std::string> snp_ids;

  // Vectors storing the number of samples per covariate class.
  Eigen::VectorXd pt_samples;

  // Vectors storing the number of cases per covariate class.
  Eigen::VectorXd pt_cases;
 
  // Empty constructor for initialization.
  Data() = default;

  // Constructor with covariate file.
  Data(std::string, std::string, std::string, std::string);

  // Constructor without covariate file.
  Data(std::string, std::string, std::string);

private:

  // Vector that contains the cumulative sum over the samples in the covariate
  // classes.
  Eigen::VectorXd pt_samples_cumsum;

  // Get initial dimensionality of a data set.
  std::tuple<int, int> get_dimensions(std::string filename);

  // Reading of the data set from text-files.
  Eigen::MatrixXd load_data(std::string filename, int n_rows, int n_cols,
                            std::vector<int> sort_vec);

  // Function to load the labels from a txt-file.
  Eigen::VectorXd load_labels(std::string filename, long long n_samples,
                              std::vector<int> sort_vec);

  // Function to load the covariates from a txt.file.
  std::vector<int> load_covar(std::string filename, long long n_samples);

  // Function to load the SNP-IDs from a txt-file.
  std::vector<std::string> load_snps(std::string filename);

  // Functions to initialize the different data metrics.
  void init_cases();
  void init_cumsum();
  void init_covar(std::vector<int> tmp_covar, std::vector<int> sort_vec);
  
};


/*
  Constructor function with covariate file.
*/
Data::Data(std::string data_fn, std::string labels_fn, std::string snp_id_fn,
           std::string covar_fn)
{

  // Check if the input files exist.
  check_file(data_fn);
  check_file(labels_fn);
  check_file(snp_id_fn);
  check_file(covar_fn);

  // get and set the data dimensions.
  std::tuple<int, int> data_dim = get_dimensions(data_fn);
  n_samples = std::get<1>(data_dim);
  n_features = std::get<0>(data_dim);

  // get the covariate dimensions.
  std::tuple<int, int> covar_dim = get_dimensions(covar_fn);

  // get the label dimensions.
  std::tuple<int, int> label_dim = get_dimensions(labels_fn);

  // Assert that the covariate dimensions match the ones in the data set.
  if (n_samples != std::get<0>(covar_dim)){
    std::cerr << "Error @ dimensions: Different number of samples in data-";
    std::cerr << "file (" << n_samples << ") and covariate-file (";
    std::cerr << std::get<0>(covar_dim) << ")." << std::endl;
    std::exit(-1);
  }

  // Assert that the label dimensions match the ones in the data set.
  if (n_samples != std::get<0>(label_dim)){
    std::cerr << "Error @ dimensions: Different number of samples in data-";
    std::cerr << "file (" << n_samples << ") and label-file (";
    std::cerr << std::get<0>(label_dim) << ")." << std::endl;
    std::exit(-1);
  }

  // Load the covariates.
  std::cout << "Loading covariates. ";
  std::vector<int> tmp_covar = load_covar(covar_fn, n_samples);
  std::cout << "Done" << std::endl;

  // Get the sorting indices of the covariates (ordering of samples).
  std::vector<int> y;
  for(int i=0; i<tmp_covar.size(); i++) 
  {
    y.push_back(tmp_covar[i]);
  }
  std::vector<int> sort_id = argsort(y);

  // Read the snps.
  std::cout << "Loading SNPs. ";
  snp_ids = load_snps(snp_id_fn);
  std::cout << "Done" << std::endl;

  // Read the labels and sort them according to the covariates.
  std::cout << "Loading labels. ";
  labels = load_labels(labels_fn, n_samples, sort_id);
  std::cout << "Done" << std::endl;

  // Read the data and sort the samples according to the covariates.
  std::cout << "Loading data. ";
  matrix = load_data(data_fn, n_features, n_samples, sort_id); 
  std::cout << "Done" << std::endl;

  // invert labels if the positive class is the majority class.
  if (labels.sum() > n_samples/2)
  {
    for (int i=0; i<labels.rows(); i++)
    {
      labels(i) = !labels(i);
    } 
  }

  n_cases = labels.sum();

  // get the per-table counts of cases and samples.
  init_covar(tmp_covar, sort_id);
  n_covar = pt_samples.size();
  init_cumsum();
  init_cases();

}


/*
  Constructor function without the covariate file. (Called at initialization of 
  class, precomputes everything and sets all default values).
*/
Data::Data (std::string data_fn, std::string labels_fn, std::string snp_id_fn)
{

  // Check if the input files exist.
  check_file(data_fn);
  check_file(labels_fn);
  check_file(snp_id_fn);

  // Get and set the data dimensions.
  std::tuple<int, int> dimensions = get_dimensions(data_fn);
  n_samples = std::get<1>(dimensions);
  n_features = std::get<0>(dimensions);

  // Init the covariates (all belong to the same class, i.e. class 0).
  std::vector<int> tmp_covar(n_samples, 0);

  // Get the sorting indices of the covariates (ordering of samples). In this,
  // case, the ordering does not change.
  std::vector<int> sort_id;
  for(int i=0; i<tmp_covar.size(); i++) 
    {
      sort_id.push_back(i);
    }

  // Read the snps.
  std::cout << "Loading SNPs. ";
  snp_ids = load_snps(snp_id_fn);
  std::cout << "Done" << std::endl;

  // Read the labels and sort them according to the covariates.
  std::cout << "Loading labels. ";
  labels = load_labels(labels_fn, n_samples, sort_id);
  std::cout << "Done" << std::endl;

  // Read the data and sort the samples according to the covariates.
  std::cout << "Loading data. ";
  matrix = load_data(data_fn, n_features, n_samples, sort_id); 
  std::cout << "Done" << std::endl;

  // invert labels if majority class.
  if (labels.sum() > n_samples/2)
  {
    for (int i=0; i<labels.rows(); i++) 
    {
      labels(i) = !labels(i);
    }
  }

  n_cases = labels.sum();

  // get the per-table counts of cases and samples.
  init_covar(tmp_covar, sort_id);
  n_covar = pt_samples.size();
  init_cumsum();
  init_cases();
}


/*
  Create the number of samples per covariate class.
*/
void Data::init_covar(std::vector<int> tmp_covar, 
      std::vector<int> sort_vec)
{
  
  // sort the tmp_covar according to the sort_vec.
  std::vector<int> sort_covar(tmp_covar.size());
  for(int i=0; i<tmp_covar.size(); i++) 
  {
    sort_covar[i] = tmp_covar[sort_vec[i]];
  }

  std::vector<int> classes;
  int tmp_count = 1;
  for (int i=1; i<tmp_covar.size(); i++)
  {
    if (sort_covar[i-1] != sort_covar[i] ) 
    {
      classes.push_back(tmp_count);
      tmp_count = 1;
    }else{
      tmp_count += 1;
    }
  }
  classes.push_back(tmp_count);
  pt_samples = Eigen::VectorXd::Zero(classes.size());
  for (int i=0; i<classes.size(); i++) 
  {
    pt_samples[i] = classes[i];
  }
}


/*
  Create the number of cases per covariate class.
*/
void Data::init_cases(){
  pt_cases = Eigen::VectorXd::Zero(n_covar);
  int k = 0;
  for (int i=0; i<labels.rows(); i++)
  {
    if (labels(i) == 1) 
      {
        pt_cases(k) += 1;
      }
    if (i+1 == pt_samples_cumsum(k))
    {
      k += 1;
    }
  }
}


/*
  Read the number of rows (=features) and columns (=samples) from file.
*/
std::tuple<int, int> Data::get_dimensions(std::string filename)
{
  
  int n_rows = 0;
  int n_cols = 0;
  double idx;
  std::ifstream f(filename);
  std::string line;
  std::stringstream ss_line;

  while(std::getline(f, line))
  {
    if (n_rows == 0) 
    {
      ss_line << line;
      while(ss_line >> idx) n_cols ++;
    }
    n_rows++;
  }
  f.close();

  return std::make_tuple(n_rows, n_cols);
}


/*
  Load the data matrix sort them according to the covariates.
*/
Eigen::MatrixXd Data::load_data(std::string filename, int n_rows, int n_cols,
    std::vector<int> sort_vec)
{

  Eigen::MatrixXd data_matrix(n_rows, n_cols);
  std::string line;
  std::stringstream ss_line;
  
  double value;
  int row_idx = 0;
  int col_idx = 0;

  std::ifstream f(filename);
  while(std::getline(f, line))
  {
    ss_line << line;
    while(ss_line >> value) 
    {
      data_matrix(row_idx, col_idx) = value;
      col_idx ++;
    }
    ss_line.clear();
    row_idx ++;
    col_idx = 0;
  }
  f.close();      

  if(row_idx != n_rows)
  {
    std::cerr << "Error @ load_data: File " << filename;
    std::cerr << " incorrectly formatted: ";
    std::cerr << row_idx << " vertices processed, expected " << n_rows << ".";
    std::cerr << std::endl;
    exit(-1);
  }  

  Eigen::MatrixXd data_matrix_sorted(n_rows, n_cols);
  for(int i=0; i<sort_vec.size(); i++)
  {
    data_matrix_sorted.col(i) = data_matrix.col(sort_vec[i]);
  }
  return data_matrix_sorted;
}


/*
  Load the labels and sort them according to the covariates.
*/
Eigen::VectorXd Data::load_labels(std::string filename, long long n_samples,
  std::vector<int> sort_vec)
{

  std::string line;
  std::stringstream ss_line;
  int idx = 0;
  bool value;
 
  std::vector<int> tmp_labels(n_samples);
  std::ifstream f(filename);
  while(std::getline(f, line))
  {
    ss_line << line;
    while(ss_line >> value) 
    {
      tmp_labels[idx] = value;
    }
    ss_line.clear();
    idx ++;
  }

  // sort the labels according to the sort_vec and return.
  Eigen::VectorXd labels = Eigen::VectorXd::Zero(n_samples);
  for(int i=0; i<n_samples; i++) labels(i) = tmp_labels[sort_vec[i]];
  return labels;
}


/*
  Read the covariates.
*/
std::vector<int> Data::load_covar(std::string filename, long long n_samples)
{

  std::string line;
  std::stringstream ss_line;
  int idx = 0;
  int value;
  std::vector<int> tmp_covar(n_samples);
 
  std::ifstream f(filename);
  while(std::getline(f, line))
  {
    ss_line << line;
    while(ss_line >> value) 
    {
      tmp_covar[idx] = value;
    }
    ss_line.clear();
    idx ++;
  }
  return tmp_covar;
}


/*
  Read the SNP-IDs from the file.
*/
std::vector<std::string> Data::load_snps(std::string filename)
{

  std::string line;
  std::stringstream ss_line;
  std::string value;
  std::vector<std::string> snps;
 
  std::ifstream f(filename);
  while(std::getline(f, line))
  {
    ss_line << line;
    while(ss_line >> value) 
    {
      snps.push_back(value);
    }
    ss_line.clear();
  }
  return snps;
}


/*
  Create the cummulative sum of number of samples per covariate class.
*/
void Data::init_cumsum()
{
  pt_samples_cumsum = Eigen::VectorXd::Zero(n_covar);
  pt_samples_cumsum(0) = pt_samples(0);
  for (int i=1; i<pt_samples_cumsum.rows(); i++)
  {
    pt_samples_cumsum(i) = pt_samples(i) + pt_samples_cumsum(i-1);
  }
}


#endif