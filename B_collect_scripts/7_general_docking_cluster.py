#!/usr/bin/env python3

###########################################################################
##
##	Peter M.U. Ung @ MSSM
##
##  v1.0  - 13.11.21
##  v2.0  - 14.11.06 - corrected a bug that generate 1 fewer ligand
##                     in each multi-ligand cluster. Now generate a
##                     SDF file with 'cluster tag'.
##  v3.0  - 14.12.03 - revamped the flow of the script so that some
##                     components can be reused
##  v3.1  - 15.11.24 - remove explicit H for figure generation
##  v4.0  - 18.10.31 - add option to use only the top X molecules
##  v4.1  - 19.05.16   updated html grid library, add isomer diplay in figure
##                      go around a bug with python3 rdkit
##  v5.0  - 20.12.05   setting if generate image/pdf file
##  v6.0  - 21.11.29   add .xz compression capability
##  v7.0  - 22.08.09   pi-pi/cation interaction; update for py3
##  v8.0  - 22.09.15   add argparse
##  v9.0  - 23.03.06   append score to title
##  v10   - 23.04.12   use precalculated clusterIdx for clustering
##
##	Purpose: Summerize docking result with Chemical similarity, Rank
##		 the results by clusters and the first set of docking inputs
##
##		 Should be able to take in all SDF from different docking 
##		 programs if processed to give the naming format.
##
##	Read in a SDF file with multiple structures generated from 
##	Docking processing (molecule name will contain ZINC ID, Rank, 
##	Score (1 decimal place), Type).
##	1) For each molecule calculate the 2D fingerprint
##	2) Cluster molecules based on Tanimoto similarity.
##	3) Sort each cluster by size of cluster (except the last group; each
##	   mol in this group has no cluster partner)
##	4) 3 Output: a) reordered Mol file according to cluster
##		     b) 2D Image of the clustered compounds
##		     c) HTML page that contact links to the 2D images
##
##	Required:	1_vina_screen_preprocess.py
##			2_vina_screen_get_top.py
##			3_vina_screen_extract_mol.pl
##			other docking processing scripts
##			  *.top$i.txt
##
##	** The name of the Molecules in the input SDF file must have
##	   this format:
##	     -- (Mol Name, Rank, Score, Type) separated by '::'
##	   	MOLECULE_NAME::RANK::SCORE::DOCK_TYPE
##
##	   e.g. ZINC464618334::22::-12.2::FRED
##
##
##	Options: Daylight-like fingerprint and MACCS keys can be used.
##
##
##	Preset:	ECFP: default radius is 4. FCFP can be used
##		(Martin YC, Muchmore S. DOI: 10.1002/qsar.200810176)
##		Similarity Choice: Tanimoto (Default)
##		Clustering method: Butina Clustering (JCICS 39 747 (1999))
##
###########################################################################

import os,sys
msg = """\n## Usage: x.py 
      -in  <+>  [ Docking inputs: sdf ] gz,bz2,xz  - accept multiple files "*.sdf,x.sdf"
      -op  < >  [ output prefix ]
  Optional:              
      -cut < >  [ Tanimoto Cutoff: float (def: 0.4) ]
      -fp  < >  [ Fingerprint: ECFP_4: 0 | Daylight: 1 | MACCS 166-bit: 2 (def: 0) ]
                [ Override FP-based clustering if non-integer value (SDF Tag) is used ]
      -pdb < >  [ To generate PyMOL session, add receptor structure: pdb ]
      -top < >  [ no. of top molecule for clustering : int ]
      -sort < > [ Sort input molecules by <tag> ]
      -app  < > [ Append Score in <tag> to Mol Title ]
      -pdf      [ Generate mol cluster pdf ]

    # and output a list of molecule png (temp.$name.png)
      [png of clustered ligands]

    e.g.>  x.py -in  vina.sdf fred.sdf ehits.sdf 
                -cut 0.4 -op output.clust-04 -pdb receptor.pdb -fp 0\n\n"""
if len(sys.argv) < 5: sys.exit(msg)

import re,gzip,bz2,lzma,glob
#import PIL.Image
#import html
#import cairo
import numpy as np

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import DrawingOptions
from rdkit.Chem.rdmolops import AssignStereochemistryFrom3D

from argparse import ArgumentParser
from rdkit_grid_print import grid_print
from tqdm import tqdm


##########################################################################
## Sort and compare to eliminate duplicate
## For duplicates, the Type, Rank and Score are combined
def main( ):

  args = UserInput()

  ## Get the docking inputs
  Dock_Files = [ glob.glob(name)[0] for name in args.filenames ]

  All_Hit, Multi_Hit = SortDockFiles(Dock_Files, args.top_mol, args.sort_tag)

  ## Make clusters of all imported results
  output_name = args.out_pref
  Mol_Cluster = FPClusterList(All_Hit, float(args.CUT), args.fingerprint )

  ## Show clusters of imported results that have at least 1 overlap from multiple files; rarely used
#  if len(Multi_Hit) > 0:
#    output_name = args.out_pref+'.multi-hit'
#    Mol_Cluster = Do_Fp(Multi_Hit, float(args.CUT), int(args.fingerprint))

  if args.pmol_pdb:
    GenPyMOLClust(Mol_Cluster, output_name, args.pmol_pdb, Dock_Files)
  if args.pdf:
    GenClustTable(Mol_Cluster, output_name, column=5)



##########################################################################
def SortDockFiles( Dock_Files, top_mol, sort_tag ):
  Mol_Hash = {}	# Hash Table of all available molecule info

  for dock_file in Dock_Files:
    dock = file_handle(dock_file)
    Dock = [ x for x in (Chem.ForwardSDMolSupplier(dock, removeHs=False))
              if x is not None ]
    if top_mol:
      if sort_tag:
        Temp = Dock
        Dock = sorted(Temp, key=lambda m: m.GetProp(sort_tag))
      temp = Dock[0:int(top_mol)]
      Dock = temp
    print(" --> Read {0}: {1}".format(dock_file,len(Dock)))

#    from rdkit.Chem import MolStandardize
#    Dock2 = Dock  # temp fix for sdf input without explicit H
#    Dock = [MolStandardize.rdMolStandardize.Normalize(Chem.AddHs(m,addCoords=True)) for m in Dock2]

    Temp = {}   # Hash Table of current $dock_file 
    for m in Dock: Temp[m.GetProp('_Name').split('::')[0]] = m

    ## Compare $Mol_Hash and $dock_file hash table
    ## Find duplicate, append the title, then remove one of the copies
    for hash_key in Mol_Hash.keys():
      if not hash_key: break
      hash_m    = Mol_Hash[hash_key]
      Hash_Info = hash_m.GetProp('_Name').split('::') # (ZINC, Rank, Score, Type)

      ## If mol is found in the Master, update the existing entry with new name
      if Temp.get(hash_key):
        temp_m    = Temp[hash_key]
        Temp_Info = temp_m.GetProp('_Name').split('::')
        new_name  = hash_key+'::'+str(Hash_Info[1]+'/'+Temp_Info[1])+'::'+str(Hash_Info[2]+'/'+Temp_Info[2])+'::'+str(Hash_Info[3]+'/'+Temp_Info[3])
        temp_m.SetProp('_Name', new_name)
        Mol_Hash[hash_key] = temp_m
        del Temp[hash_key]	# Remove the duplicate from $dock_file

    ## Unique molecules are added to the Master Hash Table as new entries
#    new_name = hash_key+'::'+str(Hash_Info[1]+'/')+'::'+str(Hash_Info[2])+'::'+str(Hash_Info[3])
#    hash_m.SetProp('_Name', new_nm)
    for temp_key in Temp.keys():
      Mol_Hash[temp_key] = Temp[temp_key]
    print("  ## Added "+str(len(Temp))+" Unique Entries from "+dock_file+" ##")

## Convert Hash back into List. Isolate entries with multiple top hits in separate input files
  All_Hit   = []
  Multi_Hit = []
  Singl_Hit = []
  for key in Mol_Hash.keys():
    m = Mol_Hash[key]
    All_Hit.append(m)
    if re.search(r'/', m.GetProp('_Name')):
      Multi_Hit.append(m)
    else:
      Singl_Hit.append(m)
  print("  ### Unique Entries Found: "+str(len(Multi_Hit)+len(Singl_Hit))+" ###\n")
  print("  ### "+str(len(Multi_Hit))+" ligands with < 1+ > occurance in input ###")
  print("  ### "+str(len(Singl_Hit))+" ligands with < 1  > occurance in input ###\n")

  return All_Hit, Multi_Hit


########################################################################
########################################################################

def FPClusterList( Mol, cutoff, fingerprint ):

  if type(fingerprint) == int:
    if int(fingerprint) == 0:
      #### Morgan (ECFP-like) fingerprint ####
      fp = 'ECPF_4'
      print('  -- Calculating with Morgan (ECFP_4) Fingerprint --')
      FP = [AllChem.GetMorganFingerprintAsBitVect(x, 2, 2048) for x in Mol]

    if int(fingerprint) == 1:
      #### Daylight-like fingerprint ####
      fp = 'Daylight'
      from rdkit.Chem.Fingerprints import FingerprintMols
      print('  -- Calculating with DayLight-like Fingerprint --')
      FP = [FingerprintMols.FingerprintMol(x) for x in Mol]

    if int(fingerprint) == 2:
      #### MACCS 166-bit keys ####
      fp = 'MACCS'
      from rdkit.Chem import MACCSkeys
      print('  -- Calculating with MACCS 166-bit Keys --')
      FP = [MACCSkeys.GenMACCSKeys(x) for x in Mol]

    # Cluster based on fingerprint
    Cluster_List = ClusterFps(FP, cutoff)
  else:
    ## Use ClusterIdx tag in SDF for clustering
    fp = fingerprint
    Cluster_List = ConvertClusTag(Mol, fingerprint)

#############
  ## Use list_position_id stored in Cluster_List to build a List of mol clusters: 
  ## Clusters = [[clust1], [clust2], [clust3], ...]
  M_Clusters = []
  for cluster in Cluster_List:
    if len(cluster) > 1:
      M_Temp = []
      for clust_id in cluster:
        M_Temp.append(Mol[clust_id])
      M_Clusters.append(M_Temp)
    else:
      M_Clusters.append([Mol[cluster[0]]])

  if type(fingerprint) == int:
    print("  ## No. of cluster by <{0}> at Tc={1}: {2}\n".format(fp, cutoff, len(M_Clusters)) )
  else:
    print("  ## No. of cluster by SDF Tag <{0}>: {1}\n".format(fp, len(M_Clusters)) )

  Mol_List = RankClusters(M_Clusters, fp)
  return Mol_List


##########################################################################
## Generate table of clustered molecules with ranking of cluster list based
## on each cluster's averaged mol-docking ranking, and add some properties to each mol
def RankClusters( M_Clusters, fingerprint ):
  Multi_list  = []		# List of all clusters
  Single_list = []	# Mol that are scored as a single cluster

  for cid, Cluster in enumerate(M_Clusters):
    # If the cluster has multiple members, rank them based on average of all ranks
    if len(Cluster) == 1:
      m = Cluster[0]
      i = m.GetProp("_Name")
      if re.search(r'::', i):
        s = i.split('::')
        r = 0.0
        if re.search(r'/', s[1]): # If Rank has '/', that means it has multiple hits
          Y = [float(x) for x in s[1].split('/')]
          r = str(np.average(Y))
        else:			# otherwise, divide the rank by 2
          r = str(float(s[1])/2.0)
        mol = MolSetProp(m, s[0],s[1],s[2],s[3],r, cid, fingerprint)
        Single_list.append([mol,s[0],s[1],s[2],s[3],r, cid, fingerprint])	#(Mol, Name, Rank, Score, Type, rank_avg, cluster, fp)
      else:
        mol = MolSetProp(m, i,'','','','',cid,fingerprint)
        Single_list.append([mol, i,'','','',0,cid,fingerprint])
    else:
      clust_list = []
      for m in Cluster:
        i = m.GetProp("_Name")
        if re.search(r'::', i):
          s = i.split('::')
          r = 0.0
          if re.search(r'/', s[1]):
            Y = [float(x) for x in s[1].split('/')]
            r = str(np.average(Y))
          else:
            r = str(float(s[1])/2.0)
          ## [mol_data, Name, Rank, Score, Type, rank_avg]
          mol = MolSetProp(m, s[0],s[1],s[2],s[3],r,cid,fingerprint)
          clust_list.append([mol,s[0],s[1],s[2],s[3],r,cid,fingerprint])
        else:
          mol = MolSetProp(m, i, '', '', '', '', cid,fingerprint)
          clust_list.append([mol, i,'','','',0,cid,fingerprint])
      Multi_list.append(clust_list)

  Multi_list.append(Single_list)

  # Sort each cluster and save only the molecules, not additional info
  Mol_List = []
  for List in Multi_list:
    List.sort(key=lambda rank: float(rank[5]))
    Mol_List.append([Item[0] for Item in List])

  return Mol_List


#####################################################################
def ConvertClusTag( Mol, clus_tag ):

  ## Sort the Mol list into hash by external cluster index
  idxhash = {}
  for m_idx, m in enumerate(Mol):
    clustidx = m.GetProp(clus_tag)

    if clustidx in idxhash.keys():
      idxhash[clustidx].append(m_idx)
    else:
      idxhash[clustidx] = [m_idx]

  ## Convert hash into list of lists; only need the Values(mol_idx)
  dic_list = list(map(list, idxhash.values()))
  idx_list = sorted(dic_list, key=len, reverse=True)
  idx_tupl = tuple(map( tuple,idx_list))
  print(idx_tupl)

  return idx_list

###################################################################
def MolSetProp( mol, name, rank, score, s_type, r_avg, cid, fp ):
  mol.SetProp('Name', name)
  mol.SetProp('Rank', rank)
  mol.SetProp('Score',score)
  mol.SetProp('Type', s_type)
  mol.SetProp('RankAvg', r_avg)
  mol.SetProp('{0}_Cluster'.format(fp), str(cid))
  return mol


##########################################################################
## Takes in a molecular list that has list of lists (clusters of molecules)
## and generate a HTML-based table
def GenPyMOLClust( Mol_List, output_name, ref_pdb, Dock_Files ):
  pymol_pml = open(output_name+".pml", 'w')
  m_out     = Chem.SDWriter(output_name+'.sdf')

  ref_name = ref_pdb.split('/')[-1].split('.pdb')[0]
  pymol_pml.write('set cartoon_oval_length, .8\nset cartoon_oval_width, 0.2\nset cartoon_rect_length, 0.8\nset cartoon_rect_width, 0.2\nset_bond stick_radius, .13, poly\nset dash_gap, 0.3\n')
  pymol_pml.write("load "+ref_pdb+", "+ref_name+"\nshow cartoon, poly\nhide lines\ncolor white, poly\ncolor cyan, org\nshow sticks, org and not resn NMA+ACE\n")

  pymol_pml.write("set_bond stick_radius, .15, "+ref_name+" and org\n")
  pymol_pml.write("create ref_lig, "+ref_name+" and org and not resn NMA+ACE\n")
  pymol_pml.write("show lines, byres poly within 5 of ref_lig\n")
  pymol_pml.write("hide sticks, "+ref_name+" and org\n")

  ## load the unclustered original data
  for dock_file in Dock_Files:
    dock_name = dock_file.split('/')[-1].split('.sdf')[0]
    pymol_pml.write("load {0}, {1}\n".format(dock_file, dock_name))
    pymol_pml.write("dist HB.all, poly, {0}, mode=2\n".format(dock_name))
#    pymol_pml.write('dist pi.all, poly, {0}, mode=5\n'.format(dock_name))

  ## write out each cluster as temp sdf to load into pymol
  for idx, Mols in enumerate(Mol_List):
    pse_sdf = Chem.SDWriter('_TEMP.clust.{0}.sdf'.format(idx+1))
    for mol in Mols: 
      pse_sdf.write(mol)
      m_out.write(mol)
    pse_sdf.close()
    pymol_pml.write("load _TEMP.clust.{0}.sdf, clust.{0}\n".format(idx+1))
    pymol_pml.write("dist HB.{0}, poly, clust.{0}, mode=2\n".format(idx+1))
#    pymol_pml.write("dist pi.{0}, poly, clust.{0}, mode=5\n".format(idx+1))

  pymol_pml.write("show sticks, org\nset valence\n")
  pymol_pml.write("hide (h. and (e. c extend 1))\n")
  pymol_pml.write("util.cbas\n")
  pymol_pml.write("center org\nzoom org\nset dash_gap, 0.25\nhide labels\n")
  pymol_pml.write("show mesh, byres poly within 4 of ref_lig\n")
  pymol_pml.write("set mesh_width, 0.1\n")
  pymol_pml.write("set light_count, 1\nset ray_opaque_background, off\n")
  pymol_pml.write("color white, poly\ncolor cyan, "+ref_name+" and org\n")
#  pymol_pml.write("color cyan, ref_lig\ncolor cyan, pi.*\nutil.cnc\n")
  pymol_pml.write("disable clust.*\ndisable HB.*\ndisable ref_lig\nenable HB.all")
  pymol_pml.write("set ray_trace_mode, 1\nset ray_trace_gain, 0.008\n")
  pymol_pml.write("set ray_trace_color, black\n")
  pymol_pml.write('set pse_export_version, 1.70\n')
  pymol_pml.write("save "+output_name+".pse\nquit\n")
  pymol_pml.flush()
  pymol_pml.close()
  m_out.flush()
  m_out.close()

  os.system("pymol -c -q -Q {0}.pml".format(output_name))
  os.system('gzip {0}.*pse {0}.*sdf'.format(output_name))
  os.system("rm ./_TEMP.clust.*.sdf")
#  os.system("rm {0}.pml".format(output_name))

##########################################################################
def GenClustTable( Mol_List, output_name, column=5 ):
  Img_Data = []
  for idx, Mols in enumerate(Mol_List):
    Img = []

    for mol in Mols:
      # Get molecule info
      m1 = mol.GetProp('Name')
      m2 = mol.GetProp('Rank')
      m3 = mol.GetProp('Score')
      m4 = mol.GetProp('Type')

      # Create tag and write out to sdf file
      mol.SetProp('Cluster', str(idx+1))
      mol.SetProp('SMILES' , Chem.MolToSmiles(mol, isomericSmiles=True))
      AssignStereochemistryFrom3D(mol)

      # Create figure using SMILES instead of 3D structure
      svg_name = '_TEMP.'+m1+'.svg'
      mol = rdMolDraw2D.PrepareMolForDrawing(mol)
      mol = Chem.RemoveHs(mol)
      AllChem.Compute2DCoords(mol)
      DrawingOptions.atomLabelFontSize=18
      
      Draw.MolToFile(mol, svg_name, size=(225,225) )
      #cairosvg.svg2png( url=svg_name, write_to=png_name, dpi=240 )
      img_link = '<img src="'+svg_name+'">'
          # Img = (image_link, Name, Rank, Score, Type)
      Img.append([img_link, m1, m2, m3, m4])

    Img_Data.append(Img)

## Print out a HTML page, in which every row has a maximum of 5 compound png.
## Every major cluster of compounds is grouped together.
## List the Name of the compound, then the Rank and Score.

  grid_print(output_name, Img_Data, 'formatted', column=5)

#  os.system('rm ./_TEMP.*.svg ./{0}.html'.format(output_name))


#######################################################################
## to handle the raw SDF format, python3's rdkit has a documented bug and
## hasn't been fixed since 2016. https://github.com/rdkit/rdkit/issues/1065
## To avoid it, the input file cannot be an object handle of a regular file,
## i.e. handle = open('xxx.sdf','r') will fail but handle = 'xxx.sdf' is fine.
## It only happens to regular file but not to gzip.open() or bz2.BZ2File() in
## python3 rdkit but not in python2 rdkit...
## Fix it by replace handle = open('xxx.sdf','r') with handle = file('xxx.sdf')
## to make regular file to be import as an object.
##
## Handle gzip and bzip2 file if the extension is right. otherwise, just open
## outuput: file handle
def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    handle = gzip.open(file_name, 'rb')
  elif re.search(r'.xz$', file_name):
    handle = lzma.open(file_name, 'rb')
  elif re.search(r'.bz2$', file_name):
    handle = bz2.BZ2File(file_name, 'rb')
  else:
    handle = open(file_name, 'rb')

#  print("## Opening "+file_name)
  return handle

#######################################################################
## Cluster Molecules based on Fingerprint
def ClusterFps(Fp, cutoff):

    # first generate the distance matrix:
    dists = []
    nfps = len(Fp)
    for i in list(range(1, nfps)):
      sims = DataStructs.BulkTanimotoSimilarity(Fp[i], Fp[:i])
      dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
#    print(cs)  # show data table formating; tuple
    return cs


########################################################################
def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-in', dest='filenames', required=True, nargs='+',
                  help='Docking SDF inputs: gz,bz2,xz  - accept multiple files "*.sdf,x.sdf" ')
  p.add_argument('-op', dest='out_pref', required=True,
                  help='Prefix of Output')

  p.add_argument('-fp', dest='fingerprint', required=False, default=0,
                  help='Fingerprint: ECFP_4: 0 | Daylight: 1 | MACCS 166-bit: 2 (def: 0)')
  p.add_argument('-cut', dest='CUT', required=False, default=0.4,
                  help='Tanimoto Cutoff: float (def: 0.4)')
  p.add_argument('-top', dest='top_mol', required=False, default=False,
                  help='no. of top molecule for clustering : int')
  p.add_argument('-clus',dest='clus_tag', required=False, default=False,
                  help='Override FP-based clustering, use clusterIdx in SDF Tag instead')
  p.add_argument('-sort', dest='sort_tag', required=False, default=False,
                  help='Sort molecules by <tag>')

  p.add_argument('-pdb', dest='pmol_pdb', required=False, default=False,
                  help='To generate PyMOL session, add receptor structure: pdb')
  p.add_argument('-pdf', dest='pdf', required=False, action='store_true',
                  help='Generate mol cluster pdf')
  p.add_argument('-app', dest='app_score', required=False, default=False,
                  help='Append Score in <tag> to Mol Title')

  args = p.parse_args()
  return args


############################################################################
if __name__ == '__main__':
  main(  )
