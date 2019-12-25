#!/usr/bin/python

# Peter M.U. Ung @ MSSM
# v1.0  15.06.20
# Clean the Schrodinger Glide result format to
# simplified format

import sys,re
read = False
Data = []
#Data.append(['Title', 'conf', 'GlideScore'])
Data.append(['Title', 'GlideScore'])

with open(sys.argv[1], 'rh') as fi:
  for l in fi:
    if read:
      Itms = l.split()
      try:
        if re.search(r'CHEM', Itms[2]):
          Data.append( [ "{0} {1}".format(Itms[1], Itms[2]), Itms[4] ] )
        else:
          Data.append( [ Itms[1], Itms[3] ] )
#        Data.append( [ Itms[1], Itms[3], Itms[5] ] )
#        Data.append([Itms[1], Itms[5]])
      except IndexError:
        print l
        read = False
        continue

    if re.search(r'====', l):
      read = True

with open(sys.argv[2], 'wh') as fo:
  for Itms in Data:
    fo.write('{0}\t{1}\n'.format(Itms[0],Itms[1]))
#    fo.write('{0} {1} {2}\n'.format(Itms[0], Itms[1], Itms[2]))

