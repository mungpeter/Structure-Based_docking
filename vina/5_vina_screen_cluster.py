#!/usr/bin/env python

###########################################################################
##
##	Peter M.U. Ung @ MSSM
##
##	v1.0	- 13.11.02
##
##	Purpose: Organize molecules by Chemical similarity, then Rank
##
##	Read in a SDF file with multiple structures generated from 
##	VINA docking processing (molecule name will contain ZINC ID, rank, 
##	score).
##	1) For each molecule calculate the 2D fingerprint
##	2) Cluster molecules based on Tanimoto similarity.
##	3) Sort each cluster by size of cluster (except the last group; each
##	   mol in this group has no cluster partner)
##	4) 2 Output: a) reordered Mol file according to cluster
##		     b) 2D Image of the clustered compounds
##
##	Required:	1_vina_screen_preprocess.py
##			2_vina_screen_get_top.py
##			3_vina_screen_extract_mol.pl
##			  *.top$i.txt
##
##		** The name of the Molecules in the input SDF file must have
##		   this format:
##		     -- (Mol Name, Rank, Score) separated by '::'
##		   	MOLECULE_NAME::RANK::SCORE
##
##	Options: Daylight-like fingerprint and MACCS keys can be used, just 
##		 need to change the code.
##
##	Preset:	ECFP: default radius is 4. FCFP can be used
##		(Martin YC, Muchmore S. DOI: 10.1002/qsar.200810176)
##		Similarity Choice: Tanimoto (Default, Dice Available)
##		Clustering method: Butina Clustering (JCICS 39 747 (1999))
##

##	This script can only handle VINA result
##
##	This script can be replaced by 6_vina_fred_cluster.py,
##	which can handle both VINA and FRED result

###########################################################################

import os,sys
msg = """\n  ## Usage: x.py [input mol: sdf] [Tanimoto Cutoff: float]
	      [top hits: int] [output prefix] 
              [receptor structure: pdb]

              # and output a list of molecule png (temp.$name.png)
              [png of clustered ligands]

              ### This Script Handle VINA results ONLY ###\n\n"""
if len(sys.argv) != 6: sys.exit(msg)

import re,gzip,bz2,rdkit,Image,HTML
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

input_SDF    = sys.argv[1]
output_SDF   = sys.argv[4]
top          = int(sys.argv[3])
CUT          = 1-float(sys.argv[2])
ref_pdb      = sys.argv[5]

#######################################################################
def ClusterFps(Fp, cutoff):
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina

    # first generate the distance matrix:
    dists = []
    nfps = len(Fp)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(Fp[i], Fp[:i])
#        sims = DataStructs.BulkDiceSimilarity(Fp[i], Fp[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs


######################################################################
sdf_handle = input_SDF
if   re.search(r'.gz',  input_SDF): sdf_handle = gzip.open(input_SDF,   'r')
elif re.search(r'.bz2', input_SDF): sdf_handle = bz2.BZ2File(input_SDF, 'r')
t_sdf = [x for x in Chem.ForwardSDMolSupplier(sdf_handle,removeHs=False) 
         if x is not None]
M_sdf = [x for x in t_sdf[:top]]
print len(M_sdf)

#### Daylight-like fingerprint
#from rdkit.Chem.Fingerprints import FingerprintMols
#Fp = [FingerprintMols.FingerprintMol(x) for x in M_sdf]

#### Morgan (ECFP-like) fingerprint
Fp = [AllChem.GetMorganFingerprintAsBitVect(x, 2, 2048) for x in M_sdf]

#### MACCS Keys
#from rdkit.Chem import MACCSkeys
#Fp = [MACCSkeys.GenMACCSKeys(x) for x in M_sdf]

Cluster_List = ClusterFps(Fp, cutoff=CUT)


#####################################################################
## Build a new Molecule list according to the clusters
M_Clusters = []
for List in Cluster_List:
  if len(List) > 1:
#    M_Temp = []
#    for num in range(1, len(List)):
#      M_Temp.append(M_sdf[List[num]])
    M_Temp = [ M_sdf[clust_id] for clust_id in List[1:] ]
    M_Clusters.append(M_Temp)
  else:
    M_Clusters.append([M_sdf[List[0]]])
print "no. of cluster: "+str(len(M_Clusters))


## Generate table of clustered molecules
Img_Data    = []
Multi_list  = []	# List of all clusters
Single_list = []	# Mol that are scored as a single cluster

for Cluster in M_Clusters:
  if len(Cluster) == 1:
    m = Cluster[0]
    i = m.GetProp("_Name")
    s = i.split('::')
    Single_list.append([m,s[0],s[1],s[2]])	#(Mol, Name, Rank, Score)
  else:
    clust_list = []
    for m in Cluster:
      i = m.GetProp("_Name")
      s = i.split('::')
      ## [mol_data, Name, Rank, Score]
      clust_list.append([m,s[0],s[1],s[2]])
    Multi_list.append(clust_list)

Multi_list.append(Single_list)


## Generate PyMOL session of clusters of ligands
pymol_nme = "temp.pml"
pymol_pml = open(pymol_nme, 'w')
pymol_pml.write("load "+ref_pdb+"\nshow cartoon, poly\nhide lines\n")
for idx, List in enumerate(Multi_list):
  List.sort(key=lambda rank: int(rank[2]))
  Img = []

  pse_sdf = Chem.SDWriter("TEMP.clust."+str(idx)+".sdf")
  for M in List:
    pse_sdf.write(M[0])

    img_name = 'TEMP.'+M[1]+'.png'
    AllChem.Compute2DCoords(M[0])
    Draw.MolToFile(M[0], img_name, size=(200,200))
    img_link = '<img src="'+img_name+'">'
        # Img = (image_link, Name, Rank, Score, Mol_file)
    Img.append([img_link, M[1], M[2], "%.1f" % float(M[3])])
  pse_sdf.flush()
  pse_sdf.close()
  pymol_pml.write("load TEMP.clust."+str(idx)+".sdf, clust."+str(idx)+"\n")
  pymol_pml.write("dist HB."+str(idx)+", poly, clust."+str(idx)+", mode=2\n")
  Img_Data.append(Img)
pymol_pml.write("show sticks, org\nshow lines, "+ref_pdb+
                " within 6 of org\n"+"save "+output_SDF+".pse\nquit\n")
pymol_pml.close()
os.system("pymol -c "+pymol_nme)
os.system("rm ./TEMP.clust.*.sdf")
os.system("rm temp.pml")

#######################################################################
## Print out a HTML page, in which every row has a maximum of 5 compound png.
## Every major cluster of compounds is grouped together.
## List the Name of the compound, then the Rank and Score.

PAGE = open(output_SDF+'.html', 'w')
for idx, C in enumerate(Img_Data):
  t = HTML.Table()
  c_temp = []	# Compound Row
  i_temp = []	# ID number Row
  s_temp = []	# Rank/Score Row
  for num, img in enumerate(C):
    if num == 0: 
      c_temp.append(img[0])
      i_temp.append('<center>'+img[1]+'</center>')
      s_temp.append('<center>Rank: '+str(img[2])+' | '+str(img[3])+'</center>')
    else:
      if num % 5 != 0: 
        c_temp.append(img[0])
        i_temp.append('<center>'+img[1]+'</center>')
        s_temp.append('<center>Rank: '+str(img[2])+' | '+str(img[3])+'</center>')
      else:
        t.rows.append(c_temp)
        t.rows.append(i_temp)
        t.rows.append(s_temp)
        c_temp = [img[0]]
        i_temp = ['<center>'+img[1]+'</center>']
        s_temp = ['<center>Rank: '+str(img[2])+' | '+str(img[3])+'</center>']

      if num == len(C)-1:
        t.rows.append(c_temp)
        t.rows.append(i_temp)
        t.rows.append(s_temp)
  
  
  htmlcode = str(t)
  if idx == len(Img_Data)-1: print PAGE.write('<p><b>  ###  No Cluster  ###</b></p>')
  print PAGE.write("<p><b>  ## Cluster: "+str(idx+1)+" : "+str(len(C))+"Hits ##</b></p>")
  print PAGE.write(htmlcode)
print PAGE.write("<p><b>  ## no. of cluster: "+str(len(M_Clusters))+" ##</b></p>")
PAGE.close()


import os
os.system("google-chrome "+output_SDF+".html")
