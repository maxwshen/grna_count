# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, imp
sys.path.append('/cluster/mshen/')
import numpy as np
from collections import defaultdict
from mylib import util, compbio
import pandas as pd


# Default params
inp_dir = _config.OUT_PLACE + 'a_split/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')

##
# Functions
##
def match(read, header):
  for idx, row in exp_design.iterrows():
    bc = row['Barcode'].upper()
    index = row['Index'].upper()
    bc_score = match_barcode(read, bc)
    idx_score = match_index(header, index)
    # print idx_score, bc_score, row['Name']
    if idx_score is True:
      if bc_score is None:
        return row['Name'], read
      elif bc_score <= 1:
        trimmed_read = read[len(bc):]
        return row['Name'], trimmed_read
  return 'other', read

def match_index(header, index):
  read_idx = header.split(':')[-1].strip()
  return bool(index == read_idx)

def match_barcode(read, bc):
  if bc != 'NONE':
    s = read[:len(bc)]
    score = sum([1 for i in range(len(bc)) if s[i] != bc[i]])
  else:
    score = None
  return score

##
# primary
##
def demultiplex(split):
  inp_fn = inp_dir + '%s.fq' % (split)
  for name in list(exp_design['Name']) + ['other']:
    util.ensure_dir_exists(out_dir + name)
    util.exists_empty_fn(out_dir + name + '/%s.fa' % (split))

  lc = util.line_count(inp_fn)
  num_bad_q, num_tot = 0, 0
  timer = util.Timer(total = lc)
  with open(inp_fn) as f:
    for i, line in enumerate(f):
      if i % 4 == 0:
        header = line.strip()
      if i % 4 == 1:
        read = line.strip()
      if i % 4 == 3:
        num_tot += 1
        qs = line.strip()
        quals = [ord(s)-33 for s in qs]
        if np.mean(quals) < 30:
          num_bad_q += 1
          continue

        demultiplex_id, trimmed_read = match(read, header)
        
        out_fn = out_dir +  '%s/%s.fa' % (demultiplex_id, split)
        with open(out_fn, 'a') as f:
          f.write('>' + header[1:] + '\n' + trimmed_read + '\n')
      
      timer.update()

  print 'Rejected %s fraction of reads' % (num_bad_q / num_tot)

  return

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

##
# Main
##
@util.time_dec
def main(split = ''):
  print NAME  

  if split == '':
    gen_qsubs()
    return

  demultiplex(split)

  return out_dir


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1])
  else:
    main()