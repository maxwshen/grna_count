import _config
import sys, os, fnmatch, datetime, subprocess
sys.path.append('/cluster/mshen/')
import numpy as np
import pandas as pd
from mylib import util
from collections import defaultdict

# Default params
NAME = util.get_fn(__file__)
inp_dir = _config.OUT_PLACE + 'c_countgrna/'
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')
spec_sum = defaultdict(dict)
with open(_config.DATA_DIR + 'spec_sum.csv') as f:
  for i, line in enumerate(f):
    if i == 0:
      continue
    w = line.strip().split(',')
    lib_nm = w[0]
    combination_nm = w[1]
    components = [s for s in w[2:] if s != '']
    spec_sum[lib_nm][combination_nm] = components

# Functions
def combine_dfs(library_nm):
  print library_nm
  mdf = None

  timer = util.Timer(total = 60)
  for split in range(60):
    df = pd.read_csv(inp_dir + '%s_%s.csv' % (library_nm, split), index_col = 0)
    if mdf is None:
      mdf = df
    else:
      cols = [s for s in mdf.columns if s not in ['gRNA sequence', 'gRNA name', 'Gene']]
      for col in cols:
        mdf[col] += df[col]
    timer.update()

  # Sum columns matching library nm
  sss = spec_sum[library_nm]
  for combination_nm in sss.keys():
    print 'Summing for %s..' % (combination_nm)
    components = sss[combination_nm]
    mdf[combination_nm] = sum(mdf[s] for s in components)

  # Write CSV
  mdf.to_csv(out_dir + '%s.csv' % (library_nm))
  
  # Write TSV for MaGeCK
  renaming = {
    'gRNA name': 'sgRNA',
    'Gene': 'gene',
  }
  mdf = mdf.rename(renaming, axis = 'columns')
  mdf = mdf.drop('gRNA sequence', axis = 'columns')
  mdf = mdf.set_index('sgRNA')
  mdf.to_csv(out_dir + '%s.tsv' % (library_nm), sep = '\t')

  return


@util.time_dec
def main():
  print NAME  

  for library in set(exp_design['Library']):
    combine_dfs(library)

  return

if __name__ == '__main__':
  main()
