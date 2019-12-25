#!/usr/bin/env python3

##########################################################################
##
##	Peter M.U. Ung @ MSSM
##
##	v1.0	- 13.12.25
##
##	Compare two result lists by scatter plot and calculate the correlation
##	between the two lists
##
##########################################################################

import sys,re
if len(sys.argv) != 5:
  MSG = '''
          ## Usage: x.py [Score file 1: txt] [Score file 2: txt]
                         [Top MOL to compare: int]
	 		 [Prefix of output plot files]
        '''
  sys.exit(MSG)

from CommonUtility import *

##########################################################################
def doit():
  File_Name     = []
  max_num       = int(sys.argv[3])-1
  prefix        = sys.argv[4]
  File_Name.append(sys.argv[1])
  File_Name.append(sys.argv[2])
  All_Data = []

  ## Read in the score.txt, store the names and scores, then rank them
  for file_name in File_Name:
    Hash = {}

    with file_handle(file_name) as f:
      Data = []
      for line in f:
        if not line.strip(): continue
        if re.search(r'Title', line): print("Reading...")
        else:
          if re.search(r'_', line):
	    name = str(line.split()[0].split('_')[0])
          else:
            name = str(line.split()[0])
          score = float(line.split()[1])
          Data.append([name, score])
      Data.sort(lambda x,y: cmp(float(x[1]), float(y[1])))
      print("  # Found {0} molecules".format(len(Data)))
      f = open(file_name+'.sort.txt','w') ###
      for idx, Item in enumerate(Data):
        f.write(Item[0]+"\t"+str(Item[1])+"\n") ###
#        Hash[Item[0]] = str(idx)	# Use Ranking = Spearman?
	Hash[Item[0]] = str(Item[1])	# Use Score value
      All_Data.append(Hash)

  k = open('scatter-hash.txt', 'w')
  Scatter = []
  Failed  = []		# Molecules not found in both lists
  count   = 0		# Number of Molecule not found in both lists
  ## Compile the comparison list
  for key in All_Data[0].keys():

    if All_Data[1].get(key):
      Scatter.append([All_Data[0][key], All_Data[1][key]])
      k.write(key+"\t"+str(All_Data[0][key])+"\t"+str(All_Data[1][key])+"\n") ###
      del All_Data[1][key]
    else:
#      print("Item > "+key+" < is not found in > "+File_Name[1]+" <")
      Failed.append([key, '1'])

  ## Plot the scatter plot data with 1st item as x-axis
  # Calculate all data
  X = [M[0] for M in Scatter]
  Y = [M[1] for M in Scatter]
  Scatter_Plot(X, Y, X_Lim=(-18,-2), Y_Lim=(-18,-2), 
               x_lab=File_Name[0], y_lab=File_Name[1], 
               title=prefix+".all", png_name=prefix+".all.scatter.png")
  # Sort the data and Calculate only the Top-selected data
  Scatter.sort(lambda x,y: cmp(float(x[0]), float(y[0])))
  I = [M[0] for M in Scatter[:max_num]]
  J = [M[1] for M in Scatter[:max_num]]
  Scatter_Plot(I, J, X_Lim=None, Y_Lim=None, 
               x_lab=File_Name[0], y_lab=File_Name[1], 
               title=prefix+'.top_'+str(max_num+1),
               png_name=prefix+".top_"+str(max_num+1)+".scatter.png")

  ## Check the reminders in the 2nd list and output the failed molecules
  for key in All_Data[1].keys():
    if key is not None:
      Failed.append([key, '0'])

  out = open(prefix+".failed.txt", 'w')
  for idx, Item in enumerate(Failed):
    out.write(str(idx)+"\tMissing: "+Item[0]+"\t"+Item[1])
    print(str(idx)+"\tMissing: "+Item[0]+"\t"+Item[1])
  out.close()
  print("  ## Molecules not found in one of the lists: "+str(len(Failed)))


##########################################################################
## Do scatter plot, take in lists of X and Y
 # X        : a list of input for the X-axis
 # Y        : a list of input for the y-axis
 # Z        : a list of input for the z-axis (color)
 # X_Lim    : the limits of x-axis, in the format of [min(X), max(X)]
 # Y_Lim    : the limits of x-axis, in the format of [min(Y), max(Y)]
 # x_lab    : label of the x-axis
 # y_lab    : label of the y-axis
 # title    : title of the scatter plot
 # png_name : name of the scatter plot PNG figure
def Scatter_Plot(X_In, Y_In, Z=None, X_Lim=None, Y_Lim=None, 
                 x_lab=None, y_lab=None, title=None, png_name=None):
  ## Default setting for the figure
  DPI   = 150
  SIZ_X = 6
  SIZ_Y = 6

  X = [float(x) for x in X_In]
  Y = [float(y) for y in Y_In]
  x_min = min(X)
  x_max = max(X)
  y_min = min(Y)
  y_max = max(Y)

  ## Calculate Pearson and Spearman correlations
  import numpy as np
  import scipy.stats as sp
  R = sp.pearsonr  (X, Y)	# (r_value,  p_value)
  S = sp.spearmanr (X, Y)	# (spearman, p_value)
  slope, intercept, r_value, p_value, std_err = sp.linregress(X, Y)
  print("Pearson:   "+str(R))
  print("Spearman:  "+str(S))
  r2 = R[0] * R[0]	# Coefficient of Determination

  ## Calculate the standard deviation of each point to the best-fit line
  fit    = np.polyfit(X, Y, 1)	# Do 1st-order polynomial fitting
  fit_fn = np.poly1d(fit)	# Create 1-D plot function
  Y_Dist = [(abs(Y[idx] - fit_fn(x)))**2 for idx, x in enumerate(X)]
  stddev = np.sqrt( np.mean(Y_Dist, dtype=np.float64) )
  print(" stdev --> "+str(stddev))
  
  ## Create Scatter Plot
  import matplotlib.pyplot as plt
  from matplotlib import rc, rcParams

  print(" ## Generating Scatter Plot ##")
  fig = plt.figure(num=1, figsize=(SIZ_X,SIZ_Y), dpi=DPI)

  ## Plot 1-D linear regressed data, the 2*stddev lines, then scatter plot
  msg_1 = (r"$y$   = "+str("%.4f" % slope)+"\n"+
           r"int = "  +str("%.4f" % intercept)+"\n"+
           r"$R^2$ = "+str("%.4f" % r2  )+"\n"+
           r"$R$   = "+str("%.4f" % R[0])+"\n"+
	   r"$S$   = "+str("%.4f" % S[0])+"\n"+
	   r"$\sigma$   = "+str("%.4f" % stddev))
#  plt.plot(X, fit_fn(X), label=msg_1,
  plt.plot(X, fit_fn(X), 
           color='purple', linewidth=1., linestyle='-')
   # 2x-Simga lines of the linear regression line
  plt.plot(X, fit_fn(X)+(2*stddev),  
           color='violet', linewidth=1.,  linestyle='-')
  plt.plot(X, fit_fn(X)-(2*stddev),
           color='violet', linewidth=.75, linestyle='-')

  if X_Lim is not None:
    x_min, x_max = X_Lim
    plt.xlim(x_min, x_max)
  if Y_Lim is not None:
    y_min, y_max = Y_Lim
    plt.ylim(y_min, y_max)

  ## Plot a diagonal line for comparison
#  plt.plot([min(X), max(X)], [min(Y), max(Y)]
  plt.plot([x_min, x_max], [y_min, y_max],
           color='black', linewidth=2., linestyle='-')

  ## If Z has no input, use the Y_Dist as the Z-data instead
  if Z is None: Z = Y_Dist
  ## Do Scatter plot, with order of data is colored according to X
  plt.scatter(X, Y, c=Z, color='blue', s=10, 
              edgecolor='black', linewidth=0.2, cmap=plt.cm.jet)

  plt.xlabel(x_lab)
  plt.ylabel(y_lab)
  plt.title (title)
  
  plt.grid(True)
  plt.legend(loc='lower right')

#  plt.show()
  plt.savefig(png_name, dpi=DPI)
  plt.close()


##########################################################################
if __name__ == '__main__':
  doit()
