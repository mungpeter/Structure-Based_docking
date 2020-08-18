#!/usr/bin/env python3

import sys
from rdkit import Chem
from rdkit.Chem import PandasTools as rdpd


df   = rdpd.LoadSDF(sys.argv[1], removeHs=False, molColName='ROMol')
s_id = df['SourceID']

for idx, row in df.iterrows():
  df['ROMol'][idx].SetProp('_Name', s_id[idx])

rdpd.WriteSDF(df, sys.argv[2], properties=list(df.columns))

##########################################################################
#
#  Peter M.U. Ung @ gRED
#
#  v1.0    20.08.10
#
#  Specifically used for gRED MOE-fastROCS search result -- 
#  Replace MolID in SDF MOL (Title) by <SourceID>
#
