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
##  v9.0    23.12.15   convert to use PandasTools, add Murcko clustering

##	Purpose: Summerize docking result with Chemical similarity, Rank
##		 the results by clusters and the first set of docking inputs
##
##		 Should be able to take in all SDF from different docking 
##		 programs if processed to give the naming format.
##
##	Read in a SDF file with multiple structures generated from 
##	Docking processing (molecule name will contain ZINC ID, Rank, 
##	Score (1 decimal place), Type).
##	1) For each molecule calculate the 2D fingerprint/generic Murcko
##	2) Cluster molecules based on Tanimoto similarity/Murcko matching
##	3) Sort each cluster by size of cluster (except the last group; each
##	   mol in this group has no cluster partner)
##	4) Output reordered Mol file according to cluster
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
##	   e.g. ZINC464618334::22::-12.2::FRED
##
##	Options: Daylight-like fingerprint and MACCS keys can be used.
##
##	Preset:	ECFP: default radius is 4. FCFP can be used
##		(Martin YC, Muchmore S. DOI: 10.1002/qsar.200810176)
##		Similarity Choice: Tanimoto (Default)
##		Clustering method: Butina Clustering (JCICS 39 747 (1999))
##
###########################################################################

import os,sys
from typing import Any
msg = """\n## Usage: x.py 
      -in  <+>  [ Docking inputs: sdf ] gz,bz2,xz  - accept multiple files "*.sdf,x.sdf"
      -pdb < >  [ Receptor structure for PyMOL session: pdb ]
      -op  < >  [ output prefix ]
  Optional:              
      -fp  < >  [ Fingerprint: ECFP_4: 0 | Daylight: 1 | MACCS 166-bit: 2 (def: 0) ]
      -cut < >  [ Tanimoto Cutoff: float (def: 0.4) ]\n
      -murcko   [ Generate Generic Murcko structure for clustering 'Murcko_SMILES' ]
      -murcko_col [ <SDF Tag> Generic Murcko structure used for Clustering (def: 'Murcko_SMILES') ]\n
      -top < >  [ no. of top molecule for clustering : int (def: None) ]
      -sort < > [ Sort molecules by <SDF tag> (def: None) ]\n
    FRED/HYBRID score tag: Chemgauss4
    GlideSP score tag:     r_i_glide_gscore\n
    e.g.>  x.py -in  vina.sdf fred.sdf ehits.sdf 
                -pdb receptor.pdb -op output.clust-04 -fp 0 -cut 0.4 \n\n"""
if len(sys.argv) < 4: sys.exit(msg)

import re
import gzip,bz2,lzma
import glob
import pandas as pd

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.ML.Cluster import Butina
from rdkit.Chem.Fingerprints import FingerprintMols

from argparse import ArgumentParser
from rdkit.Chem import PandasTools as rdpd


##########################################################################
def main( ):

  args = UserInput()

  ## Get the docking inputs
  Dock_Files = [ glob.glob(name)[0] for name in args.filenames ]
  sdf = GetMolFromSDF( Dock_Files, args.sort_col, args.top_mol )

  if args.murcko_gen:
    tmp_rept, tmp_sing = MurckoClustering( sdf, args.murcko_col )
  else:
    tmp_rept, tmp_sing = FPClustering( sdf, args.cutoff, args.fingerprint )

  ## Make clusters of all imported results
  all_clusters = SortClusters( tmp_rept, tmp_sing, args.sort_col )

  ## Show clusters of imported results that have at least 1 overlap from multiple files; rarely used
  GenPyMOLClust( all_clusters, args.out_pref, args.pmol_pdb, Dock_Files )


##########################################################################
def GetMolFromSDF( Dock_Files, sort_col, top_mol ):

  all_dock = []
  for dock_file in Dock_Files:
    df = rdpd.LoadSDF(file_handle(dock_file), molColName='ROMol', idName='ID', removeHs=False)
    df['ID'] = df['ID'].astype(str)		# deal with numeric molid
    all_dock.append(df)
  print(" --> Read {0}: {1}".format(dock_file, len(df)))
  sdf = pd.concat(all_dock)

  ## (Name::Rank::Score::Type)
  if re.search(df.ID[0], '::'):
    sdf['Rank']  = sdf.ID.apply(lambda m: m.split('::')[1]).astype(int)
    sdf['Score'] = sdf.ID.apply(lambda m: m.split('::')[2]).astype(float)

  if sort_col:
    sdf[sort_col] = sdf[sort_col].astype(float)
    temp = sdf.sort_values(by=sort_col).reset_index(drop=True)
    print(temp[sort_col])
    if top_mol:
      return temp[:int(top_mol)]
    else:
      return temp


########################################################################
class FPChoice(object):
  def __init__( self, choice=0 ):
    self.choice = choice

  def __call__( self, mol ):
    return self._fp_calculator(mol)

  def _fp_calculator( self, mol ):
    if self.choice == 0:
      #### Morgan (ECFP-like) fingerprint ####
      return AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
    if self.choice == 1:
      #### Daylight-like fingerprint ####
      return FingerprintMols.FingerprintMol(mol)
    if self.choice == 2:
      #### MACCS 166-bit keys ####
      return MACCSkeys.GenMACCSKeys(mol)

#############################
## Cluster Molecules based on Fingerprint
def CalcFPTanimoto( Fp, cutoff ):
  ## Generate the FP distance matrix:
  dists = []
  for i in list(range(1, len(Fp))):
    sims = DataStructs.BulkTanimotoSimilarity(Fp[i], Fp[:i])
    dists.extend([1-x for x in sims])

  ## Cluster the data into list of lists of positional index (df or array), 
  ## each list with members within the cutoff:
  cs = Butina.ClusterData(dists, len(Fp), float(cutoff), isDistData=True)
  return cs


##########################################################################
def FPClustering( sdf, cutoff, fingerprint ):
  ## Build list of cluster lists from FP clustering = [[clusts1], [clusts2], ...]
  ## positional index from cluster used to extract members into groups
  FP_calc  = FPChoice(choice=fingerprint)
  Clusters = CalcFPTanimoto( sdf['ROMol'].apply(FP_calc), cutoff )
  tmp_rept = [ sdf[sdf.index.isin(cluster)].reset_index(drop=True).drop(columns=['index'])
                for cluster in Clusters if len(cluster) >= 2 ]
  temp     = [ sdf[sdf.index.isin(cluster)] for cluster in Clusters if len(cluster) == 1 ]
  tmp_sing = pd.concat(temp).reset_index(drop=True).drop(columns=['index'])

  print("## No. of cluster: "+str(len(tmp_rept)+len(tmp_sing))+" for Tc = "+str(cutoff)+" ##\n")
  return tmp_rept, tmp_sing


##########################################################################
def MurckoClustering( sdf, murcko_col ):
  ## Summarize to unique Murcko and how many member in each unique Murcko
  group = sdf.groupby([murcko_col]).sum()
  group['murcko_num'] = sdf.groupby([murcko_col]).size()

  ## first identify singleton Murkco, then combine all singletons as one group, 
  ## slice these singletons from the original data, sort them by 'sort_col'
  singleton = group[group.murcko_num == 1].reset_index()
  scaffolds = [ item[1] for item in singleton[murcko_col].items() ]
  collects  = [ sdf[sdf[murcko_col] == murcko ].reset_index(drop=True) for murcko in scaffolds]
  tmp_sing  = pd.concat(collects).reset_index(drop=True).drop(columns=['index'])

  ## Slice Repeats into groups of clusters, Sort list of clusters by: 
  ## cluster size, then by cluster with better score first
  repeats   = group[group.murcko_num  > 1].reset_index()
  scaffolds = [ item[1] for item in repeats[murcko_col].items() ]
  tmp_rept  = [ sdf[sdf[murcko_col] == murcko ].reset_index(drop=True).drop(columns=['index'])
                for murcko in scaffolds ]

  return tmp_rept, tmp_sing


##########################################################################
def SortClusters( tmp_rept, tmp_sing, sort_col ):

  ## if "sort_col" is used, sort by "Cluster size"+"Score"; if not using 
  ## "sort_col", just sort by "Cluster size" alone
  if sort_col:
    sort_rept = sorted(tmp_rept, key=lambda x: (len(x), -x[sort_col].astype(float).min()), reverse=True)
    sort_sing = tmp_sing.sort_values(by=sort_col).reset_index(drop=True)
    all_clust = sort_rept + [sort_sing]
  else:
    sort_rept = sorted(tmp_rept, key=lambda x: len(x), reverse=True)
    all_clust = sort_rept + [tmp_sing]

  print("  ### "+str(len(tmp_rept))+" Clusters  with < 1+ > entry members in input ###")
  print("  ### "+str(len(tmp_sing))+" Singleton with < 1  > entry member  in input ###\n")

  return all_clust


##########################################################################
## Takes in a molecular list that has list of lists (clusters of molecules)
## and generate a HTML-based table
def GenPyMOLClust( Mol_List, output_name, ref_pdb, Dock_Files ):
  pymol_pml = open(output_name+".pml", 'w')

  ref_name = ref_pdb.split('/')[-1].split('.pdb')[0]
  pymol_pml.write("set cartoon_oval_length, .8\nset cartoon_oval_width, 0.2\nset cartoon_rect_length, 0.8\n")
  pymol_pml.write("set cartoon_rect_width, 0.2\nset_bond stick_radius, .13, poly\nset dash_gap, 0.3\n")
  pymol_pml.write("load "+ref_pdb+", "+ref_name+"\nshow cartoon, poly\nhide lines\nhide sticks\n")
  pymol_pml.write("show sticks, resname HOH+WAT\nhide sticks, "+ref_name+" and org\n")
  pymol_pml.write("create ref_lig, "+ref_name+" and org and not resn NMA+ACE\n")
  pymol_pml.write("show lines, ref_lig\nhide sticks, ref_lig\nutil.cnc\n")
  pymol_pml.write("set_bond stick_radius, .15\nshow lines, byres poly within 5 of ref_lig\n")

  ## load the unclustered original data
  for idx,dock_file in enumerate(Dock_Files):
    dock_name = dock_file.split('/')[-1].split('.sdf')[0]
    pymol_pml.write("load {0}, {1}\n".format(dock_file, dock_name))
    pymol_pml.write("show sticks, {0}\ncolor salmon, {0}\n".format(dock_name))
    pymol_pml.write("dist HB.all.{0}, poly, {1}, mode=2\n".format(idx+1, dock_name))

  ## write out each cluster as temp sdf to load into pymol
  for idx, Mols in enumerate(Mol_List):
    rdpd.WriteSDF(Mols, '_TEMP.clust.{0}.sdf'.format(idx+1), molColName='ROMol', idName='ID')
    pymol_pml.write("load _TEMP.clust.{0}.sdf, clust.{0}\n".format(idx+1))
    pymol_pml.write("show sticks, clust.{0}\ncolor salmon, clust.{0}\n".format(idx+1))
    pymol_pml.write("dist HB.{0}, poly, clust.{0}, mode=2\n".format(idx+1))

  pymol_pml.write("set valence\nhide (h. and (e. c extend 1))\n")
  pymol_pml.write("color white, poly\ncolor cyan, ref_lig\n")
  pymol_pml.write("color cyan, "+ref_name+" and org\nutil.cnc\n")
  pymol_pml.write("center ref_lig\nzoom ref_lig\nset dash_gap, 0.25\nhide labels\n")
  pymol_pml.write("show mesh, byres poly within 4 of ref_lig\nset mesh_width, 0.2\n")
  pymol_pml.write("disable clust.*\ndisable HB.*\ndisable ref_lig\nenable HB.all")

  pymol_pml.write("set light_count, 1\nset ray_opaque_background, off\n")
  pymol_pml.write("set ray_trace_mode, 1\nset ray_trace_gain, 0.008\nset ray_trace_color, black\n")
  pymol_pml.write("save "+output_name+".pse\nquit\n")
  pymol_pml.flush()
  pymol_pml.close()

  os.system("pymol -c -q -Q {0}.pml".format(output_name))
  os.system("rm ./_TEMP.clust.*.sdf {0}.pml".format(output_name))
  os.system('gzip -f {0}.pse'.format(output_name))


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

########################################################################
def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-in', dest='filenames', required=True, nargs='+',
                  help='Docking SDF inputs: gz,bz2,xz  - accept multiple files "*.sdf,x.sdf" ')
  p.add_argument('-pdb', dest='pmol_pdb', required=True,
                  help='To generate PyMOL session, add receptor structure: pdb')
  p.add_argument('-op', dest='out_pref', required=True,
                  help='Prefix of Output')

  p.add_argument('-fp', dest='fingerprint', required=False, default=0,
                  help='Fingerprint: ECFP_4: 0 | Daylight: 1 | MACCS 166-bit: 2 (def: 0)')
  p.add_argument('-cut', dest='cutoff', required=False, default=0.4,
                  help='Tanimoto Cutoff: float (def: 0.4)')
  p.add_argument('-top', dest='top_mol', required=False, default=False,
                  help='no. of top molecule for clustering : int')

  p.add_argument('-murcko', dest='murcko_gen', required=False, action='store_true',
                  help='Generate Generic Murcko structure "Murcko_SMILES"')
  p.add_argument('-murcko_col', dest='murcko_col', required=False, default='Murcko_SMILES',
                  help='<SDF Tag> Generic Murcko structure used for Clustering')

  p.add_argument('-sort', dest='sort_col', required=False, default=False,
                  help='Sort molecules by <SDF tag>')

  args = p.parse_args()
  return args


############################################################################
if __name__ == '__main__':
  main(  )
