# SiNIMin

This repository contains the code for the *significant Network Interval Mining* approach, short SiNIMin, described in _Network-guided detection of candidate intervals that exhibit genetic heterogeneity_ (under review).


## Data formatting.

Assuming we are given a data set of n samples with d binary features. An example of all files can be found in the folder _examples_.

The method requires the following input:
  1. __data file__ with d rows corresponding to features (important: the features are assumed to follow a natural ordering, such as genetic variants by their position on the DNA ) and n columns, corresponding to n samples. The values are supposed to be binary.
  2. __label file__ with n rows, that contains the binary phenotype of the n samples. Samples are assumed to be in same ordering as in data file.
  3. __feature file__ d rows, that contains the name of the d features. Samples are assumed to be in same ordering as in data file.
  4.  __edge file__, where each row contains the names of the nodes adjacent to the edge in tab-separated format.
  5. __mapping file__, linking the features to the nodes in the network. Each row contains the name of the node followed by a white-space separated list of feature names.
  6. __target FWER__, the target family wise error rate, default: 0.05.
  7. __covariate file__ (optional) with n rows. Each row contains the index of the class of the corresponding sample. Samples are assumed to be in same ordering as in data file.
  
  

## Usage information

### Compilation
Note that the package relies on the Eigen-library. This library has to be linked upon re-compilation of the method.
OpenMP is used for parallelization of permutation testing.

### Example usage.

An example on how to use the code can be found in the examples folder.
The executable for the method is called `sinimin` and can be found in the _SiNIMin/compiled_.

```bash
./sinimin \
  -i "${data_file}" \
  -l "${labels_file}" \
  -c "${covariate_file}" \
  -m "${mapping_file}" \
  -e "${edge_file}" \
  -s "${feature_file}" \
  -f 0.05 \
  -o "${output_prefix}" \
 ```
There exist additional flags that can be set, namely:
```bash
  -d ${maxlen} \
  -n ${number_threads} \
  -p ${number_permutations} 
```

The `-d` flag toggles the maximum length of intervals to be tested. For example, if `d` is set to 1, only interactions between single features are tested.
The `-p` flag toggles the number of permutations. If this flag is set, SiNIMin-WY is executed, i.e. Westfall-Young permutations are used to estimate family-wise error rates.
The `-n` flag sets the number of processes. This parameter only results in a speed-up for permutation testing. `sinimin` uses OMP to parallelize.

## Contact
anja.gumpinger@bsse.ethz.ch
