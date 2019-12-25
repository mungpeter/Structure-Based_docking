#!/usr/bin/env python3

import sys
msg = "\n\n  ## Usage: x.py [sdf file: sdf] [cluster cutoff: real] \nif cluster cutoff is 'x', no clustering occurs\n"
if len(sys.argv) != 3: sys.exit(msg)

import re,gzip,bz2
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit_grid_print import grid_print

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


input_name = sys.argv[1]

handle = input_name
if re.search(r'.gz',  input_name): handle = gz.open(input_name)
if re.search(r'.bz2', input_name): handle = bz2.BZ2File(input_name)

SDF = []
if re.search(r'.sdf', input_name):
  Temp = [x for x in Chem.ForwardSDMolSupplier(handle,removeHs=False)
          if x is not None]
  print("  ## Input is SDF format ##\n  # No. of molecules: "+str(len(Temp)))
  SDF = Temp
if re.search(r'.smi', input_name):
  Temp = [x for x in Chem.SmilesMolSupplier(handle,titleLine=False) if x is not None]
  print("  ## Input is SMI format ##\n  # No. of molecules: "+str(len(Temp)))
  SDF = Temp


## Build a new Molecule list according to the clusters
M_Clusters = []
if sys.argv[2] == 'x':
  print("  ## No clustering is specified ##")
  M_Clusters.append(SDF)
else:
  Fp = [AllChem.GetMorganFingerprintAsBitVect(x, 3, 2048) for x in SDF]
  Cluster_List = ClusterFps(Fp, cutoff=(1-float(sys.argv[2])))

  for List in Cluster_List:
    if len(List) > 1:
      M_Temp = [ SDF[clust_id] for clust_id in List[1:] ]
      M_Clusters.append(M_Temp)
    else:
      M_Clusters.append([SDF[List[0]]])
  print("no. of cluster: "+str(len(M_Clusters)))

## Formatting for print
Multi_list  = []        # List of all clusters
Single_list = []        # Mol that are scored as a single cluster

for Cluster in M_Clusters:
  if len(Cluster) == 1:
    m = Cluster[0]
    i = m.GetProp("_Name").split()[0]
    a = 'TEMP.'+i+'.svg'
    l = '<img src="'+a+'">'
    h = Chem.RemoveHs(m)
    mol = h
    rdMolDraw2D.PrepareMolForDrawing(mol)
    AllChem.Compute2DCoords(mol)
    Draw.MolToFile(m, a, size=(200,200))

    Single_list.append([l,i,"-","-",'-'])      #(Mol, Name, Rank, Score)
  else:
    clust_list = []
    for m in Cluster:
      i = m.GetProp("_Name").split()[0]
      a = 'TEMP.'+i+'.svg'
      l = '<img src="'+a+'">'
      h = Chem.RemoveHs(m)
      mol = h
      rdMolDraw2D.PrepareMolForDrawing(mol)
      AllChem.Compute2DCoords(mol)
      Draw.MolToFile(m, a, size=(200,200))

      ## [mol_data, Name, Rank, Score, Type]
      clust_list.append([l,i,"-",'-','-'])
    Multi_list.append(clust_list)

Multi_list.append(Single_list)


grid_print(sys.argv[1], Multi_list, 'formatted')
