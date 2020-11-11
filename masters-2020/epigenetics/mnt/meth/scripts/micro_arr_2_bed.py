#!/usr/bin/env python
# Author roman.cherniatchik@jetbrains.com.com

import pandas as pd
import numpy as np
import sys
import argparse

def _cli(args):
  if not args:
    print("Usage: micro_arr_2_bed.py INPUT MAPPING OUTPUT")
    exit(0)

  input_file = args[0]
  mapping_file = args[1]
  output = args[2]
  parser = argparse.ArgumentParser("Convert Microarray data to BED like file")
  parser.add_argument('input', help="Micro array data (2 columns)")
  df_mapping = pd.read_csv(mapping_file, sep="\t").dropna()
  df_sample = pd.read_csv(input_file, sep="\t")
  df = pd.merge(df_mapping, df_sample, left_on="probeID", right_on="ID_REF").drop("ID_REF", axis=1)
  df['CpG_beg'] = df['CpG_beg'].astype('int')
  df['CpG_end'] = df['CpG_end'].astype('int')
  df.rename({'CpG_chrm': 'chr', 'CpG_beg': 'start', 'CpG_end': 'end', df_sample.columns[1]: 'Mvalue'}, axis=1, inplace=True)

  x=np.power(2, df['Mvalue'])
  df['Beta']=x/(x+1)
  df.loc[:, ('chr', 'start', 'end', 'probeID', 'Beta', 'probe_strand', 'Mvalue')].to_csv(output, sep="\t", index=None, header=False)

if __name__ == "__main__":
  _cli(sys.argv[1:])