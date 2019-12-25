#!/usr/bin/env python3

##########################################################################
##
##	Peter M.U. Ung @ MSSM
##	
##	v1.0	- 14.03.20
##      v2.0    - 19.12.19  updated
##
##	Search for similar compounds in a library using ECFP_4 (default)
##	or DayLight fingerprint. 
##	DayLight fingerprint seems to be slower and less accurate than
##	ECFP_4 tho.
##
##########################################################################

import sys
msg = '''
    > {0} 
        [list of compound: sdf|smi] 
        [library for search: sdf|smi]
        [output prefix]
        [Similarity Cutoff: float] 
        [DayLight: dl | ECFP_4: ec | MACCS: mc | total: all]\n
        '''.format(sys.argv[0])
if len(sys.argv) != 6: sys.exit(msg)

import gzip,bz2,re,glob
import pandas as pd

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops
from rdkit.Chem import MACCSkeys
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import DrawingOptions
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.rdmolops import AssignStereochemistryFrom3D

from rdkit_grid_print import grid_print

##########################################################################

def main():
  Cmpd_File = sys.argv[1].split(',')
  Lib_File  = sys.argv[2].split(',')
  out_pref  = sys.argv[3]
  cutoff    = float(sys.argv[4])
  fp_choice = sys.argv[5]

  Cmpd    = rdkit_open(Cmpd_File)
  Lib     = rdkit_open(Lib_File)
  Cmpd_FP = calculate_FP(Cmpd, fp_choice)
  Lib_FP  = calculate_FP(Lib,  fp_choice)

  Selection, Save = pick_similar_cmpd(Cmpd_FP, Lib_FP, cutoff, fp_choice)
  grid_print(out_pref, Selection, 'formatted')


#########
  df = pd.DataFrame(Save, columns=['name','mol']).drop_duplicates(subset='name',keep='last')

  fs = Chem.SmilesWriter(out_pref+'.smi')
  for mol in df.to_numpy():
    fs.write(mol[1])
  fs.close()


##########################################################################
## Calculate the FingerPrint with Daylight, ECFP_4, MACCS, and total
def calculate_FP(Mols, choice):

  FP = [[mol, 
         FingerprintMols.FingerprintMol(mol),
         AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048),
         MACCSkeys.GenMACCSKeys(mol)] for mol in Mols]

  if not re.search(r'dl|ec|mc|all', choice):
    sys.exit('\n  ## ERROR: Fingerprint selection: dl ec mc all\n')
  return FP


##########################################################################
## Compare the scaffold compound to the library compounds and add the 
## library compound to the output list if it passes the cutoff threshold
def pick_similar_cmpd(Cmpd_FP, Lib_FP, cutoff, choice):
  dl_r, ec_r, ms_r = 2,2,1
  Result, Save = [], []

  for Cmpd in Cmpd_FP:
    Print = []
    cmpd          = Cmpd[0]
    cmpd_name     = cmpd.GetProp('_Name').split('::')[0]
    cmpd_img_name = '_TEMP.{0}.svg'.format(cmpd_name)
    cmpd_img_link = '<img src="{0}">'.format(cmpd_img_name)

    AssignStereochemistryFrom3D(cmpd)
    rdMolDraw2D.PrepareMolForDrawing(cmpd)
    cmpd = Chem.RemoveHs(cmpd)
    AllChem.Compute2DCoords(cmpd)
    DrawingOptions.atomLabelFontSize=18

    Draw.MolToFile(cmpd, cmpd_img_name, size=(225,225))
    Print.append([cmpd_img_link, cmpd_name, '', '','Query'])
#    Save.append([cmpd_name,cmpd])  # query mol not output to analog selection

    for Lib in Lib_FP:
      Tc = [ DataStructs.FingerprintSimilarity(Cmpd[x], Lib[x]) 
               for x in range(1,4)]

      if choice == 'dl':
        coeff = Tc[0]
      if choice == 'ec':
        coeff = Tc[1]
      if choice == 'mc':
        coeff = Tc[2]
      if choice == 'all':
        coeff = ( Tc[0]*dl_r + Tc[1]*ec_r + Tc[2]*ms_r )/(dl_r+ec_r+ms_r)

      ## Create the list if compounds are similar enough
      if coeff >= cutoff:
        mol          = Lib[0]
        lib_name     = mol.GetProp('_Name').split('::')[0]
        lib_img_name = '_TEMP.{0}.svg'.format(lib_name)
        lib_img_link = '<img src="{0}">'.format(lib_img_name)

        AssignStereochemistryFrom3D(mol)
        rdMolDraw2D.PrepareMolForDrawing(mol)
        mol = Chem.RemoveHs(mol)
        AllChem.Compute2DCoords(mol)
        DrawingOptions.atomLabelFontSize=18

        Draw.MolToFile(mol, lib_img_name, size=(225,225))
        Print.append([lib_img_link, lib_name, choice, '{:.3f}'.format(coeff), 'analog'])
        Save.append([lib_name,mol])

    if len(Print): 
      Result.append(Print)

  
  print('## No. of Query found: {0}'.format(len(Result)))
  for idx, x in enumerate(Result):
    print('# Members in Query {0}: {1}'.format(idx+1, len(x)))
  print('')
  return Result, Save


##########################################################################

## Handle gzip and bzip2 file if the extension is right. otherwise, just open
## outuput: file handle
def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    handle = gzip.open(file_name, 'rb')
  elif re.search(r'.bz2$', file_name):
    handle = bz2.BZ2File(file_name, 'rb')
  else:
    handle = (file_name)
  return handle

#######################################################################
## new version of rdkit distinguish the input source of the file, treating
## regular utf-8 file as str input and bytes file (zipped) as object input.
## Forward_supplier only takes bytes files and Regular_supplier takes regulars.
## To get around this, use file('x.sdf') to make utf-8 file as an object.

def rdkit_open(File_Tuple):

  List = []

  for f in (File_Tuple):
    handle = file_handle(f)

    if re.search(r'.sdf', f):
      Mol = [x for x in Chem.ForwardSDMolSupplier(handle, removeHs=True)
             if x is not None]

    if re.search(r'.smi', f):
      with open(f, 'r') as fi:
        first_line = fi.readline()

      if re.search(r'smiles', first_line, re.IGNORECASE):
        Mol = [x for x in Chem.SmilesMolSupplier(handle, titleLine=True,
                 delimiter=' |\t|,') if x is not None]
      else:
        Mol = [x for x in Chem.SmilesMolSupplier(handle, titleLine=False,
                 delimiter=' |\t|,') if x is not None]

    ## not the official RDkit function, may fail
    if re.search(r'.mol2', f):
      Mol = [x for x in Mol2MolSupplier(f, removeHs=False) if x is not None]

    print( "\n# Found mol in {0}: {1}\n".format(f, len(Mol)))
    for mol in Mol: List.append(mol)

  return List


##########################################################################
if __name__ == '__main__':
    main()

