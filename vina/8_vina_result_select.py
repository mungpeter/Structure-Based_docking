#!/usr/bin/python

##########################################################################
#
#   Peter M.U. Ung @ MSSM
#
#   v1.0    -- 13.??.??
#   v2.0    -- the script will sort the selected molecules by ranking/score
#
#
##########################################################################

import re,sys

from rdkit import Chem
from rdkit_open import *
from rdkit_grid_print import grid_print
from CommonUtility import *

usage = '''\n\n  ## Usage: {0}
             [list name: txt] 
             [sdf file: sdf | sep by ',']
             [Sorting: None|Name|Rank|Score]\n'''.format(sys.argv[0])
if len(sys.argv) != 4: sys.exit(usage)

list_name = sys.argv[1]
Chemicals = sys.argv[2].split(',')
if re.search(r'None', sys.argv[3], re.IGNORECASE):
  option = None
elif re.search(r'Name', sys.argv[3], re.IGNORECASE): option = 'name'
elif re.search(r'Rank', sys.argv[3], re.IGNORECASE): option = 'rank'
elif re.search(r'Score', sys.argv[3], re.IGNORECASE):option = 'score'
else:
  sys.exit('\n  ERROR: [Sorting:None|Name|Rank|Score] -- {0} \n'.format(sys.argv[3]))


##########################################################################
def main(list_name, Chemicals, option):
  ## Read in the list of selected ligand ID 
#  List = remove_remark(file_handle(list_name))
  List = [line.split()[0] for line in remove_remark(file_handle(list_name))]
  print len(List)

  ## Extract the selected ligands from the supplied SDF
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

  Mols = [mol[0] for mol in Molecules]
  out = Chem.SDWriter(list_name.split('.txt')[0] + '.sdf')
  for molecule in Mols: out.write(molecule)
  out.flush()
  out.close()

  grid_print(list_name.split('.txt')[0], Mols, 'sdf')


##########################################################################
if __name__ == '__main__':
  main(list_name, Chemicals, option)
