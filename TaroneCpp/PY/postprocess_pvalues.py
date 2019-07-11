#!/usr/bin/env python3
# @Author: guanja
# @Date:   2019-07-11 09:51:24
# @Last Modified by:   guanja
# @Last Modified time: 2019-07-11 09:59:43

import argparse
import logging
import os
import pdb

import numpy as np
import pandas as pd
import scipy.stats as st


logging.basicConfig(level='INFO', format='.. %(message)s')


def genomic_inflation_chi2(pv_lst):
  """Takes a list of pvalues from chi-squared statistic as input, and 
  computes the genomic inflation factor. 
  
  Remark: 
    survival function 'isf' is monotonically decreasing, i.e. we can 
    interchange median and survival function computation, thus speeding the 
    analysis up.

  Args:
    pv_lst: list of p-values (have to come from chi2 test, otherwise invalid
      genomic inflation returned).

  Returns:
    genomic inflation: observed median test-statistic / 
      expected median test statistic
  """
  obs_med =  st.chi2.isf(np.median(pv_lst), 1)
  exp_med =  st.chi2.isf(0.5, 1)

  return obs_med/exp_med


def main():

  parser = argparse.ArgumentParser(description='Postprocessing of pvalues' \
                                               'found in an edge-' \
                                               'epistasis run.')
  parser.add_argument('--pv_file', dest='pv_file', type=str, 
                      help='File containing pvalues.')
  parser.add_argument('--tar_file', dest='tar_file', type=str, 
                      help='File containing the testability threshold.')
  parser.add_argument('--out_file', dest='out_file', type=str, 
                      help='Output file.')
  args = parser.parse_args()


  # obtain the tarone threshold of a run.
  with open(args.tar_file, 'r') as fin:
    for line in fin:
      if line.startswith('Testability'):
        treshold = float(line.strip().split()[-1])

  # load the pvalues.
  df = pd.read_csv(args.pv_file, sep=',')
  pvalue_lst = df['p-value'].values

  lambda_gc = genomic_inflation_chi2(pvalue_lst)

  with open(args.out_file, 'w') as fout:
    fout.write(f'Genom. inflation: {lambda_gc}\n')

  pass

if __name__ == '__main__':
  main()