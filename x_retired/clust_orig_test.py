#!/usr/bin/python

##########################################################################
#
#   Peter M.U. Ung @ MSSM
#   v1.0    - 2014 ?Feb/March/April?
#   
#   Original script of dataset clustering with original Histrogram lines.
#   Clustering with ECFP_6 with Tanimoto/Dice.
#   Clustered result will be presented in html/table
#
##########################################################################

import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import Image
import HTML


input_SDF = sys.argv[1]
CUT = 1-float(sys.argv[2])

ECFP = 6    # ECFP level: default ECFP_6

##########################################################################
def ClusterFps(Fp, cutoff):
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina

    # first generate the distance matrix:
    dists = []
    nfps = len(Fp)
    for i in range(1,nfps):
#        sims = DataStructs.BulkTanimotoSimilarity(Fp[i], Fp[:i])
        sims = DataStructs.BulkDiceSimilarity(Fp[i], Fp[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs

##########################################################################

M_sdf = [x for x in Chem.ForwardSDMolSupplier(input_SDF,removeHs=False)
         if x is not None]
print len(M_sdf)

Fp = [AllChem.GetMorganFingerprintAsBitVect(x, ECFP/2, 2048) for x in M_sdf]

Cluster_List = ClusterFps(Fp,CUT)

M_Clusters = []
for List in Cluster_List:
  if len(List) > 1:
    M_Temp = [ M_sdf[clust_id] for clust_id in List[1:] ]
    M_Clusters.append(M_Temp)
  else:
    M_Clusters.append([M_sdf[List[0]]])
print "no. of cluster: "+str(len(M_Clusters))

##########################################################################
## Generate table of clustered molecules
Img_Data    = []
Multi_list  = []        # List of all clusters
Single_list = []        # Mol that are scored as a single cluster

for Cluster in M_Clusters:
  if len(Cluster) == 1:
    m = Cluster[0]
    i = m.GetProp("_Name").split()[0]
    s = i.split('::')
    Single_list.append([m,s[0]])      #(Mol, Name, Rank, Score)
  else:
    clust_list = []
    for m in Cluster:
      i = m.GetProp("_Name").split()[0]
      s = i.split('::')
      ## [mol_data, Name, Rank, Score]
      clust_list.append([m,s[0]])
    Multi_list.append(clust_list)

Multi_list.append(Single_list)

##########################################################################
## Generate html/table of clusters of ligands
for idx, List in enumerate(Multi_list):
  Img = []
  for M in List:

    img_name = 'temp.'+M[1]+'.png'
    AllChem.Compute2DCoords(M[0])
    Draw.MolToFile(M[0], img_name, size=(200,200))
    img_link = '<img src="'+img_name+'">'
        # Img = (image_link, Name, Rank, Score, Mol_file)
    Img.append([img_link, M[1]])
  Img_Data.append(Img)

PAGE = open(input_SDF+'.html', 'w')
for idx, C in enumerate(Img_Data):
  t = HTML.Table()
  c_temp = []   # Compound Row
  i_temp = []   # ID number Row
  for num, img in enumerate(C):
    if num == 0:
      c_temp.append(img[0])
      i_temp.append('<center>'+img[1]+'</center>')
    else:
      if num % 5 != 0:
        c_temp.append(img[0])
        i_temp.append('<center>'+img[1]+'</center>')
      else:
        t.rows.append(c_temp)
        t.rows.append(i_temp)
        c_temp = [img[0]]
        i_temp = ['<center>'+img[1]+'</center>']

      if num == len(C)-1:
        t.rows.append(c_temp)
        t.rows.append(i_temp)

  htmlcode = str(t)
  if idx == len(Img_Data)-1: print PAGE.write('<p><b>  ###  No Cluster  ###</b></p>')
  print PAGE.write("<p><b>  ## Cluster: "+str(idx+1)+" : "+str(len(C))+"Hits ##</b></p>")
  print PAGE.write(htmlcode)
print PAGE.write("<p><b>  ## no. of cluster: "+str(len(M_Clusters))+" ##</b></p>")
PAGE.close()
import os
os.system("google-chrome "+input_SDF+".html")
