#!/usr/bin/env python3

##########################################################################
#
#   Peter M.U. Ung @ MSSM
#
#   v1.0    -- 13.??.??
#   v2.0    -- the script will sort the selected molecules by ranking/score
#
#
##########################################################################

import sys
usage = '''\n\n  ## Usage: {0}
              [list name: txt] 
              [sdf file(s): sdf | sep by ',']
              [Sorting: None|Name|Rank|Score]\n'''.format(sys.argv[0])
if len(sys.argv) != 4: sys.exit(usage)

import re,glob
import bz2,gzip
import pandas as pd

from rdkit import Chem
from rdkit_grid_print import grid_print


list_name = sys.argv[1]
Chemicals = sys.argv[2].split(',')
if re.search(r'None', sys.argv[3], re.IGNORECASE):
  option = None
elif re.search(r'^Name$', sys.argv[3], re.IGNORECASE): option = 'name'
elif re.search(r'^Rank$', sys.argv[3], re.IGNORECASE): option = 'rank'
elif re.search(r'^Score$', sys.argv[3], re.IGNORECASE):option = 'score'
else: option = sys.argv[3]

print(option)
##########################################################################
def main(list_name, Chemicals, option):

  ## Read in the list of selected ligand ID 
  df   = pd.read_csv(list_name, delimiter='\s+', header=None, comment='#').dropna()
  List = df.loc[:, 0].to_numpy()
  print('\n > Number of items in <{}>: {}\n'.format( list_name, len(List) ))


  ## Extract the selected ligands from the supplied SDF
  print(' > List of structure file(s) read: \n',Chemicals)
  temp = rdkit_open(Chemicals)

  sdf = dict()
  for m in temp:
    if re.search(r'::', m.GetProp('_Name')):
      name, rank, score, x = m.GetProp('_Name').split('::')
      sdf[name] = [m, name, rank, score]
    else:
      name = m.GetProp('_Name')
      sdf[name] = [m, name, 0, 0.0]
  Molecules = [sdf[chem] for chem in List if chem is not None]

  ## Sort data, if needed
  if option is not None:
    if   option == 'name': Molecules.sort(key=lambda tup: tup[1])
    elif option == 'rank': Molecules.sort(key=lambda tup: int(tup[2]))
    elif option == 'score':Molecules.sort(key=lambda tup: float(tup[3]))
    else:
      print(' ## Using SDF tag to sort ligand order: \033[31m{0}\033[0m\n'.format(option))
      Molecules.sort(key=lambda tup: float(tup[0].GetProp(option)))

  Mols = [mol[0] for mol in Molecules]
  out = Chem.SDWriter(list_name.split('.txt')[0] + '.sdf')
  for molecule in Mols: 
    out.write(molecule)
  out.flush()
  out.close()

  grid_print(list_name.split('.txt')[0], Mols, 'sdf')



##########################################################################
## Handle gzip and bzip2 file if the extension is right. otherwise, just open
## outuput: file handle
def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    handle = gzip.open(file_name, 'rb')
  elif re.search(r'.bz2$', file_name):
    handle = bz2.BZ2File(file_name, 'rb')
  else:
    handle = open(file_name, 'r')

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
      Mol = [x for x in Chem.ForwardSDMolSupplier(handle, removeHs=False)
              if x is not None]

    if re.search(r'.smi', f):
      with handle as fi:
        first_line = fi.readline()

      if re.search(r'smiles', first_line, re.IGNORECASE):
        Mol = [x for x in Chem.SmilesMolSupplier(f, titleLine=True,
                  delimiter=' |\t|,') if x is not None]
      else:
        Mol = [x for x in Chem.SmilesMolSupplier(f, titleLine=False,
                  delimiter=' |\t|,') if x is not None]

    print( "\n# Found mol in {0}: {1}\n".format(f, len(Mol)))
    for mol in Mol: List.append(mol)

  return List


###########################################################################
if __name__ == '__main__':
  main(list_name, Chemicals, option)
