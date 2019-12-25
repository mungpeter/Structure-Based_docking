#!/usr/bin/env python3

## 15.11.24
## compare 2 files to see which ligands are seen in the 
## compared files

import re,sys
from CommonUtility import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_open import *

msg = '''\n  > {0}
        [file 1] [file 2]
        [SDF of 1 or 2]
        [Output prefix]\n'''.format(sys.argv[0])
if len(sys.argv) != 5: sys.exit(msg)

#######################################################

with file_handle(sys.argv[1]) as fi:
  file1 = [l.split()[0] for l in fi]
with file_handle(sys.argv[2]) as fi:
  file2 = [l.split()[0] for l in fi]

match = set(file1).intersection(file2)

chems = {}
mols = rdkit_open([sys.argv[3]])
for m in mols:
  name = m.GetProp('_Name').split('::')[0]
  chems[name] = m

fo = open(sys.argv[4]+'.txt', 'wh')
fo.write('# {0}\t{1}\n# {2}\t{3}\n'.format(
  sys.argv[1], len(file1), sys.argv[2], len(file2)))
fo.write('# match: {0}\n\n'.format(len(match)))


fm = Chem.SDWriter(sys.argv[4]+'.sdf')
for m in match:
  if re.match(r'#', m): continue
  fo.write(m+'\n')
  fm.write(chems[m])

fo.close()
