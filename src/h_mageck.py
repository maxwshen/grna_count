import _config
import sys, os, fnmatch, datetime, subprocess, math
import numpy as np
import pandas as pd
from mylib import util
from collections import defaultdict

# Default params
inp_dir = _config.OUT_PLACE + 'c2_combine_counts/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
raw_out_dir = _config.OUT_PLACE + NAME + '/raw/'
util.ensure_dir_exists(out_dir)
util.ensure_dir_exists(raw_out_dir)

spec = pd.read_csv(_config.DATA_DIR + 'spec_mageck.csv')

PSEUDOCOUNT = 1

# Functions
def call_mageck(lib_nm, counts_fn, control_nm, treatment_nm):
  tool = '/cluster/mshen/tools/mageck-0.5.6/bin/mageck'

  vs_nm = '%s___%s__vs__%s' % (lib_nm, control_nm, treatment_nm)
  command = '%s test -k %s -c %s -t %s -n %s' % (tool, counts_fn, control_nm, treatment_nm, vs_nm)
  print command
  subprocess.check_output(command, shell = True)
  
  working_dir = os.getcwd() + '/'
  subprocess.check_output('mv %s%s* %s' % (working_dir, vs_nm, raw_out_dir), shell = True)

  # Rename columns
  def rename_gene_summary(s):
    if s in ['id', 'num']:
      return s
    s = s.replace('neg', control_nm)

    s = s.replace('pos', treatment_nm)
    return 'enriched_in_%s' % (s)

  def rename_sgrna_summary(s):
    if s in ['sgrna', 'Gene', 'score', 'p.low', 'p.high', 'p.twosided', 'FDR', 'adj_var']:
      return s
    if s == 'LFC':
      return 'LFC_%s_over_%s' % (treatment_nm, control_nm)

    s = s.replace('control', control_nm)

    s = s.replace('treatment', treatment_nm)
    s = s.replace('treat', treatment_nm)

    return 'enriched_in_%s' % (s)

  df = pd.read_csv(raw_out_dir + vs_nm + '.gene_summary.txt', delimiter = '\t')
  df.columns = [rename_gene_summary(s) for s in list(df.columns)]
  df.to_csv(out_dir + vs_nm + '.gene_summary.csv')

  df = pd.read_csv(raw_out_dir + vs_nm + '.sgrna_summary.txt', delimiter = '\t')
  df.columns = [rename_sgrna_summary(s) for s in list(df.columns)]
  df.to_csv(out_dir + vs_nm + '.sgrna_summary.csv')

  return





@util.time_dec
def main():
  print NAME  

  # Function calls
  for idx, row in spec.iterrows():
    control_nm = row['Control']
    treatment_nm = row['Treatment']
    lib_nm = row['Library']

    counts_fn = inp_dir + '%s.tsv' % (lib_nm)
    call_mageck(lib_nm, counts_fn, control_nm, treatment_nm)

  return 


if __name__ == '__main__':
  main()
