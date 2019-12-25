#!/usr/bin/env python3

import os,sys,glob
msg = '''\n  ## Usage: {0}
                [Input ligands :sdf] - accept multiple files "*.sdf"
                [SDF info keyword]   - keyword for cluster ID
                [Output prefix] 
                [receptor structure: pdb]\n\n'''.format(sys.argv[0])
if len(sys.argv) != 5: sys.exit(msg)

from rdkit_open import rdkit_open
from rdkit_grid_print import grid_print
from rdkit import Chem
from rdkit.Chem import AllChem


##########################################################################

def main(in_file, keyword, out_prefix, ref_pdb):
  Mol = rdkit_open(glob.glob(in_file))

  Collect = {}  
  for m in Mol:
    clust_id = m.GetProp(keyword)
    if Collect.has_key(clust_id):
      Collect[clust_id].append(m)
    else:
      Collect[clust_id] = [m]

  Temp  = []
  Singl = []
  for clust_id in Collect.keys():
    if len(Collect[clust_id]) > 1:
      Temp.append([int(clust_id), Collect[clust_id]])
    elif len(Collect[clust_id]) == 1:
      Singl.append(Collect[clust_id][0])
  Temp.append([len(Collect)+100, Singl])  
  
  Temp.sort(key=lambda rank: rank[0])
  Clusters = [Mol[1] for Mol in Temp]

  GenPyMOLClust(Clusters, out_prefix, ref_pdb)
  GenClustTable(Clusters, out_prefix, column=5)


########################################################################
if __name__ == '__main__':
  main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
