#!/usr/bin/python

# Peter M.U. Ung @ MSSM
# v1.0  15.06.20
# Clean the Schrodinger Glide result format to
# simplified format

import sys,re

if len(sys.argv) != 3:
  sys.exit('\n   > {0}\n          [glide .rept] [output]\n'.format(sys.argv[0]))

read = False
Data = []
Data.append(['Title', 'GlideScore'])
with open(sys.argv[1], 'rh') as fi:
  for l in fi:
    if read:
      Itms = l.split()
      try:
        Data.append([Itms[1], Itms[3]])
      except IndexError:
        print l
        read = False
        continue

    if re.search(r'====', l):
      read = True

with open(sys.argv[2], 'wh') as fo:
  for Itms in Data:
    fo.write('{0}\t{1}\n'.format(Itms[0],Itms[1]))

