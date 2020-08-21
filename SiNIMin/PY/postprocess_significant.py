#!/usr/bin/env python3
# @Author: guanja
# @Date:   2019-07-11 09:22:29
# @Last Modified by:   Anja Gumpinger
# @Last Modified time: 2019-10-11 09:29:18

import argparse
import logging
import os
import ipdb

import numpy as np
import pandas as pd

from sinimin_utils import genes


logging.basicConfig(level='INFO', format='.. %(message)s')


def main():

  parser = argparse.ArgumentParser(description='Postprocessing of significant' \
                                               'hits found in an edge-' \
                                               'epistasis run.')
  parser.add_argument('--sig_file', dest='sig_file', type=str, 
                      help='File containing significant hits.')
  parser.add_argument('--map_file', dest='map_file', type=str, 
                      help='File containing gene2snp mapping.')
  parser.add_argument('--out_file', dest='out_file', type=str, 
                      help='Output file.')
  args = parser.parse_args()


  # read the file containing the significant hits into a pd.DataFrame
  df = pd.read_csv(args.sig_file, sep=',')

  df.rename(index=str, columns={"start0_len0": "interval0",
                                "start1_len1": "interval1"}, inplace=True)

  # read the mapping.
  mapping_geneview = genes.read_snp2gene_map_geneview(args.map_file, 
                                                      position_col=False)

  

  df_filt = pd.DataFrame(columns=df.columns)


  for idx, row in df.iterrows():

    # get the starting SNP of the interval in the first gene.
    int0 = [int(x) for x in row.interval0.split('_')]
    start0 = mapping_geneview[row.gene0][int0[0]]
    len0 = int0[1]

    # get the starting SNP of the interval in the second gene.
    int1 = [int(x) for x in row.interval1.split('_')]
    start1 = mapping_geneview[row.gene1][int1[0]]
    len1 = int1[1]

    row.interval0 = (start0, len0)
    row.interval1 = (start1, len1)
    
    df_filt = df_filt.append(row, ignore_index=True)

  df_filt.sort_values(by='p-value', inplace=True)
  df_filt.to_csv(args.out_file, index=False, sep='\t')



  pass

if __name__ == '__main__':
  main()
