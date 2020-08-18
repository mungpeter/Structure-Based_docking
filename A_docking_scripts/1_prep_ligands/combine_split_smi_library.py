#!/usr/bin/env python3

import sys,os

msg = '''\n  > {0}
    -l <>   [ list of all ZINC smiles (path + size in bytes) ]
    -p <>   [ prefix of output files ]\n
  Optional:
    -m <>   [ maximum smiles in each file (def: 100 000) ]
    -f <>   [ maximum file size in byte (def: 500 000 000 = 500 MB) ]\n
# use multiple cpu (6-12) to handle gzip operations in background
# collect smiles files+path  and size with Unix functions'''.format(sys.argv[0])
if len(sys.argv) < 5: sys.exit(msg)

import re
import subprocess
from argparse import ArgumentParser

##########################################################################
def main():

  args = UserInput()
  outpref = args.outpref
  if args.max_smi:
    max_smi = int(args.max_smi)
  else:
    max_smi     = 100000
  if args.max_size:
    max_size = int(args.max_size)
  else:
    max_size   = 500000000    # 500 MB max
  max_size_2 = max_size * 0.9

  ## read in list of smiles files, with file size (in bytes) next to it
  files = []
  with open(args.smi_list, 'r') as fi:
    for l in fi:
      tmp = l.rstrip().split()
      files.append([tmp[0], int(tmp[1])])

  ## i = current list number, num = combined file number, 
  ## start_num = separated out file number, total = current total file size
  ## stat = statement of smiles file names to combine
  i, num, start_num = 0, 1, 1  # changing linearly
  total, stat = 0, []    # constantly recycled
  
  while( i < len(files) ):

    print(' {0} -- remaining {1}'.format(i, len(files)-i) )

    ## if the next smiles file is big on its own, process it instead
    if files[i][1] > max_size_2:
      tmp_smi = '{0}.pre_{1}.smi'.format(outpref, num)
      start_num = ConcatFiles( files[i][0], tmp_smi, outpref,
                                max_smi, start_num )
      i += 1
      num += 1
      continue

    if total < max_size:
      stat.append(files[i][0])
      total += files[i][1]
      i += 1

    ## if total size is already maxed out, or reach end of input files
    if (total > max_size) or i == len(files)-1:
      tmp_smi = '{0}.pre_{1}.smi'.format(outpref, num)
      start_num = ConcatFiles( ' '.join(stat), tmp_smi, outpref,
                                max_smi, start_num )
      total, stat = 0, []
      num += 1


##########################################################################
def ConcatFiles( stat, tmp_smi, outpref, max_smi, start_num ):

  os.system('cat {0} | grep -v "smiles" | shuf > {1}'.format(stat, tmp_smi))
  print(tmp_smi)
  line_num = int(subprocess.check_output('cat {0} | wc -l'.format(tmp_smi), 
                shell=True ).decode('UTF-8'))
  print(line_num)

  i = line_num
  x = start_num
  while (i > 0):
    print(' x = {0}  --  i = {0}'.format(x, i))
    print('tail -{0} {1} | head -{2} > {3}.{4}.smi'.format(
          i, tmp_smi, max_smi, outpref, x ))
    subprocess.run('tail -{0} {1} | head -{2} > {3}.{4}.smi; gzip {3}.{4}.smi &'.format(
          i, tmp_smi, max_smi, outpref, x ), shell=True )
    i = i - max_smi
    x = x + 1

  os.system('rm '+tmp_smi)
  return x


##########################################################################
def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-l', dest='smi_list', required=True,
                  help='list of smiles file, with path and size in bytes')
  p.add_argument('-p', dest='outpref', required=True,
                  help='Output prefix')

  p.add_argument('-m', dest='max_smi', required=False,
                  help='Maximum number of SMILES string in each file (def: 100,000)')
  p.add_argument('-f', dest='max_size', required=False,
                  help='Maximum size of file (def: 500,000,000 = 500 MB)')

  args=p.parse_args()
  return args


##########################################################################
if __name__ == '__main__':
  main( )

##########################################################################
#
#  Peter M.U. Ung @ gRED
#
#  v1    20.07.16
#
#  concatenate database (mostly ZINC) smiles files up to certain maximum 
#  size, to generate multiple concatenated smiles files.
#  The concatenated files are each splitted into smaller trunks with a max.
#  number of lines in it. The lines are shuffled in the process.
#  The output sub files are numbered, and the number will continue with 
#  each splitted concatenated files
#
#  don't use "sort -R", it is hash-based and much slower than "shuf"
#  but "shuf" takes up momory, so make sure to have enough mem
#
#  concat_1.smi --> sub.1.smi, sub.2.smi, sub.3.smi ... sub.i.smi
#  concat_2.smi --> sub.i+1.smi, sub.i+2.smi. ... sub.i+j.smi
#  concat_3.smi --> sub.i+j+1.smi, ... sub.i+j+k.smi
#
#  > ll */*/*smi | tr -s ' ' | cut -d' ' -f9,5 | awk '{print $2 " " $1}' \
#                | cut -d' ' -f1,2  > allfiles.list
#
