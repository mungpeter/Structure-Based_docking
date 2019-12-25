#!/usr/bin/env python3

##########################################################################
##
##  Peter M.U. Ung @ MSSM
##
##  v1.0    --  12.05.2014
##
##  Purpose: Split a list of SDF files into sets of molecules larger/smaller
##           than the Molecular Weight cutoff.
##
##########################################################################
import os,sys,glob

msg = '''\n  Usage: {0}
           [SDF Files] [MW Cutoff]\n
      e.g.: {0} "*.sdf.bz2" 300\n'''.format(sys.argv[0])
if len(sys.argv) != 3: sys.exit(msg)

from rdkit_open import rdkit_open
from rdkit import Chem
from rdkit.Chem import Descriptors

##########################################################################
def main(in_files, cutoff):
  InFiles = glob.glob(in_files)
  print(InFiles)
  for file_name in InFiles:
    Small = []  # Smaller than cutoff MW
    Large = []  # Larger than cutoff MW
    prefix = file_name.split('.sdf')[0]

    Mol = rdkit_open([file_name])
    for mol in Mol:
      if float(Descriptors.MolWt(mol)) < cutoff:
        Small.append(mol)
      else:
        Large.append(mol)

    fs  = Chem.SDWriter(prefix+'.small.sdf')
    fss = open(prefix+'.small.txt', 'w')
    PrintFiles(fss, fs, Small)
    fs.flush()
    fss.flush()
    fs.close()
    fss.close()

    fl  = Chem.SDWriter(prefix+'.large.sdf')
    fls = open(prefix+'.large.txt', 'w')
    PrintFiles(fls, fl, Large)
    fl.close()
    fls.close()  


##########################################################################
def PrintFiles(list_handle, sdf_handle, Mol):
  for mol in Mol:
    score = mol.GetProp('FRED Chemgauss4 score')
    mname = mol.GetProp('_Name')
    mw    = Descriptors.MolWt(mol)
    list_handle.write('{0}\t{1}\t{2}\n'.format(mname, score, mw))
    sdf_handle.write(mol)


##########################################################################
if __name__ == '__main__':
  main(sys.argv[1], float(sys.argv[2]))
