# counts the number of copies of each gRNA present in reads
#   Primer:
#     Allows 1bp insertion or deletion before gRNA
#     allows any mismatches before gRNA  
#   gRNA sequence
#     requires exact match to gRNA
#   

import _config
import sys, os, fnmatch, datetime, subprocess
sys.path.append('/cluster/mshen/')
import numpy as np
import pandas as pd
from mylib import util


# Default params
NAME = util.get_fn(__file__)
inp_dir = _config.OUT_PLACE + 'b_demultiplex/'
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)


exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')

##
# Functions
##
def get_library(library_nm):
  library_fold = '/cluster/mshen/prj/mera_tiling/libraries/'
  lib = pd.read_csv(library_fold + '%s.csv' % (library_nm))
  return lib

def find_grna(exp, grna_set, seq):
  # Finds which gRNA perfectly matches the sequence
  # Assumes that a sequence can only have 1 gRNA (no overlap)
  #
  # Ensure gRNAs is a set, so 'in' is O(1)
  # 'in' is O(n) for lists
  num_found = 0
  
  start_idx = len('TGTGGAAAGGACGAAACACC')
  bp_relax = 1
  starts = range(start_idx - bp_relax, start_idx + bp_relax + 1)
  lens = [21, 20, 19]
  seqs = []
  added_once = False
  for st in starts:
    for l in lens:
      s = seq[st : st + l]
      if s in grna_set:
        return s
  return False

def count_grna(exp, lib, split):


  reads_fn = inp_dir + exp + '/%s.fa' % (split)

  # Handle potential duplicates in designed gRNAs by placing counts only in the first occurrence 
  grna_set = set(lib['gRNA sequence'])
  grna_list = list(lib['gRNA sequence'])
  idxs = dict()
  for grna in grna_set:
    idxs[grna] = grna_list.index(grna)

  # Init list to be joined to lib dataframe in same order as gRNA sequence
  counts = [0] * len(grna_list)
  
  tot = 0
  num_reads_matched = 0
  timer = util.Timer(total = util.line_count(reads_fn))
  with open(reads_fn) as f:
    for i, line in enumerate(f):
      if i % 2 == 0:
        header = line.strip()
      else:
        read = line.strip()
        matched_grna = find_grna(exp, grna_set, read)
        if matched_grna is not False:
          counts[idxs[matched_grna]] += 1
          num_reads_matched += 1
        tot += 1
      timer.update()

  try:
    pct_reads_matched = float(num_reads_matched) / tot
  except ZeroDivisionError:
    pct_reads_matched = np.nan
  print num_reads_matched, '/', tot
  print pct_reads_matched * 100, '%', ' reads matched'

  return counts


##
# qsub
##
def gen_qsubs():
  # Generate qsub shell scripts and commands for easy parallelization
  print 'Generating qsub scripts...'
  qsubs_dir = _config.QSUBS_DIR + NAME + '/'
  util.ensure_dir_exists(qsubs_dir)
  qsub_commands = []

  num_scripts = 0
  for idx in range(0, 60):
    command = 'python %s.py %s' % (NAME, idx)
    script_id = NAME.split('_')[0]

    # Write shell scripts
    sh_fn = qsubs_dir + 'q_%s_%s.sh' % (script_id, idx)
    with open(sh_fn, 'w') as f:
      f.write('#!/bin/bash\n%s\n' % (command))
    num_scripts += 1

    # Write qsub commands
    qsub_commands.append('qsub -m e -wd %s %s' % (_config.SRC_DIR, sh_fn))

  # Save commands
  with open(qsubs_dir + '_commands.txt', 'w') as f:
    f.write('\n'.join(qsub_commands))

  print 'Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir)
  return


@util.time_dec
def main(split = ''):
  print NAME

  if split == '':
    gen_qsubs()
    return

  # Function calls
  for library_nm in set(exp_design['Library']):
    print library_nm
    exps = exp_design[exp_design['Library'] == library_nm]['Name']

    lib = get_library(library_nm)
    df = lib
    for exp in exps:
      print exp
      counts = count_grna(exp, lib, split)
      df[exp] = counts

    df.to_csv(out_dir + '%s_%s.csv' % (library_nm, split))

  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(split = sys.argv[1])
  else:
    main()