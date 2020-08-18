#!/usr/bin/python

import sys,os
from rdkit import Chem
from CommonUtility import *

msg = '''   > {0}
        [Vina Score File: txt] [SDF Output Name]\n'''.format(sys.argv[0])
if len(sys.argv) != 3: sys.exit(msg)

vina_score = sys.argv[1]
sdf_name = sys.argv[2]

Extract = []
with file_handle(vina_score) as fi:
  fo   = Chem.SDWriter(sdf_name)
  fail = open(sdf_name+'.fail', 'w')
  for line in fi: 
    Items = line.split()
    
    conf_dir = Items[1]
    Temp     = Items[2].split('.')
    lig_dir  = Temp[0]+'.'+Temp[1]
    lig_name = Items[2]
    score    = Items[3]

    lig_path = conf_dir+'/'+lig_dir+'/'+lig_name
    print lig_path
    os.system('obabel -ipdbqt {0} -osdf -f 1 -l 1 > temp.sdf'.format(lig_path))

    mol = Chem.SDMolSupplier('temp.sdf',removeHs=False)[0]

    try:
      mol.SetProp('Directory', conf_dir+'/'+lig_dir+'/'+lig_name)
    except AttributeError:
      print 'Failed: '+lig_name
      fail.write(lig_name+'\n')
      continue
    mol.SetProp('Score', score)
    name = mol.GetProp('_Name')
    mol.SetProp('_Name', name+'::'+score+'::'+'vina')
    fo.write(mol)
  fo.flush()
  fo.close()
