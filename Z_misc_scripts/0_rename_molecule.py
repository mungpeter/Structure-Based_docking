#!/usr/bin/env python3

##########################################################################
##
##  Peter M.U. Ung @ MSSM
##  
##  v.1     -   15.03.02
##
##  Purpose: change the Title Name of the molecules (in SDF) to 
##           certain name stored in the data section.
##
##########################################################################

import os,sys
from CommonUtility import *
from rdkit_open import *
from rdkit import Chem
from rdkit.Chem import AllChem

msg = '''\n    > {0}
        [Mol File: sdf] ['Property' with the New Name]
        [Append name to 'Property']
        [Output Prefix] [Output Format: sdf, smi]\n
    e.g.> x.py file.sdf 'CAS' 'cas_' new_file sdf
'''.format(sys.argv[0])
if len(sys.argv) != 6: sys.exit(msg)

lig = rdkit_open([sys.argv[1]])

# get name from the property
for m in lig:
  m.SetProp('_Name', sys.argv[3]+m.GetProp(sys.argv[2]))

# specific which format writer to use
if sys.argv[5] == 'sdf':
  m_out = Chem.SDWriter(sys.argv[4]+'.sdf')
if sys.argv[5] == 'smi':
  m_out = Chem.SmilesWriter(sys.argv[4]+'.smi')

for m in lig: m_out.write(m)
m_out.flush()
m_out.close()
