#!/usr/bin/env python3

import sys,os,re
import matplotlib.pyplot as plt

def run(ref_lig, query, out_pref):

###########################################
  Ref_q = {}
  Ref_l = []
  with open(ref_lig, 'r') as fi:
    for idx, l in enumerate(fi):
      itms = l.split()
      Ref_q[itms[0]] = [ idx, itms[1] ]
      Ref_l.append([idx, itms[0], itms[1], itms[2] ])

  Query = []
  num   = 0
  with open(query, 'r') as fi:
    for l in fi :
      if not re.search(r'#', l):
        lig = l.rstrip().split()[0]
        if lig in Ref_q:
          Query.append( Ref_q[l.rstrip().split()[0]] )
        else:
          num += 1
          print('{0} not found in Refernece: {1}'.format(lig, num) )

  for x in range(len(Ref_l)-len(Query)):
    Query.append([0])


##########################################
  X = [ m[0] for m in Ref_l ]
  Y = [ m[0] for m in Query ]
  C = [ m[2] for m in Ref_l ]

  print(C)
  Scatter_Plot( X, Y, C, x_lab=ref_lig, y_lab=query, title=out_pref,
                          png_name=out_pref+'.scatter.png' )


def Scatter_Plot( X_In, Y_In, Color=None, x_lab=None, y_lab=None,
                              title=None, png_name=None ):
  DPI   = 150
  SIZ_X = 6
  SIZ_Y = 6

  fig = plt.figure(num=1, figsize=(SIZ_X,SIZ_Y), dpi=DPI)
  
  print(len(X_In))
  
  print(len(Y_In))

  plt.scatter(X_In, Y_In, c=Color, s=10, edgecolor='black' )

  plt.title(title)
  plt.grid(True)
  
  plt.savefig(png_name, dpi=DPI)
  plt.close()


###############################################
if __name__ == '__main__':
  run(sys.argv[1], sys.argv[2], sys.argv[3])