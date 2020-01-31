#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import QED
from rdkit.Chem import PandasTools
from rdkit.Chem import MolStandardize
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import SaltRemover
from rdkit.Chem import rdchem

msg = '''  > {}
\t\t[sdf filename]\n
\t\tWarning: SDF file should not exceed 1.3GB; that will use ~ 25GB RAM\n'''.format(sys.argv[0])
if len(sys.argv) != 2: sys.exit(msg)

def main(prm_file):

  pref = prm_file.split('.sdf')[0]

  print('## Reading file...')
  prm_df = PandasTools.LoadSDF(prm_file, smilesName='SMILES',
            molColName='MOL', includeFingerprints=False)
  print(prm_df[:10])

  ## remove salts and rename the smiles
  print('## Cleaning moleucles...')
  remover = SaltRemover.SaltRemover()
  prm_df['mol'] = prm_df.MOL.apply(remover.StripMol)
  prm_df['smiles'] = prm_df.mol.apply(Chem.MolToSmiles)
  prm_df['ID'] = prm_df.idnumber
  print(prm_df[:10])

  ## recalculate molecular properties
  print('## Calculating properties...')
  prm_df['qed']  = prm_df.mol.apply(QED.properties)
  prm_df['MW']   = prm_df.qed.apply(lambda x: x.MW)
  prm_df['logP'] = prm_df.qed.apply(lambda x: x.ALOGP)
  prm_df['HBA']  = prm_df.qed.apply(lambda x: x.HBA)
  prm_df['HBD']  = prm_df.qed.apply(lambda x: x.HBD)
  prm_df['PSA']  = prm_df.qed.apply(lambda x: x.PSA)
  prm_df['ROTB'] = prm_df.qed.apply(lambda x: x.ROTB)
  prm_df['AROM'] = prm_df.qed.apply(lambda x: x.AROM)
  prm_df['HA']   = prm_df.mol.apply(rdchem.Mol.GetNumHeavyAtoms)
  print(prm_df[:10])
  print(' > number of molecules... ',len(prm_df))

  ## shuffle
  print('## Shuffling molecules...')
  prm_df = prm_df.reindex(np.random.permutation(prm_df.index))

  ## print out molecule properties and smiles (shuffled)
  print('## Writing results...')
  Cols_csv = ['ID','MW','HA','logP','LogS','HBA','HBD','PSA','ROTB','AROM','SaltType','smiles']
  Cols_smi = ['smiles','ID']

  prm_df.loc[(prm_df.MW > 150.) & (prm_df.MW <= 300.)].to_csv(pref+'.frag.csv.bz2',sep=',',float_format='%.2f',columns=Cols_csv, index=False)
  prm_df.loc[(prm_df.MW > 150.) & (prm_df.MW <= 300.)].to_csv(pref+'.frag.smi',sep='\t',columns=Cols_smi,index=False)

  prm_df.loc[(prm_df.MW > 300.) & (prm_df.MW <= 400.)].to_csv(pref+'.lead.csv.bz2',sep=',',float_format='%.2f',columns=Cols_csv, index=False)
  prm_df.loc[(prm_df.MW > 300.) & (prm_df.MW <= 400.)].to_csv(pref+'.lead.smi',sep='\t',columns=Cols_smi,index=False)

  prm_df.loc[prm_df.MW > 400.].to_csv(pref+'.drug.csv.bz2',sep=',',float_format='%.2f',columns=Cols_csv, index=False)
  prm_df.loc[prm_df.MW > 400.].to_csv(pref+'.drug.smi',sep='\t',columns=Cols_smi,index=False)

  prm_df.loc[prm_df.MW <= 150.].to_csv(pref+'.small.csv.bz2',sep=',',float_format='%.2f',columns=Cols_csv, index=False)
  prm_df.loc[prm_df.MW <= 150.].to_csv(pref+'.small.smi',sep='\t',columns=Cols_smi,index=False)

#############################
if __name__ == '__main__':
  main(sys.argv[1])

##########################################################################
#
#  Peter MU Ung @ MSSM/Yale
#
#  v1  20.01.30
#
#  take in the SDF molecule library and split it into fragment/lead/drug-like
#  sets according to their Molecular Weight:
#    small <= 150
#    150 < frag <= 300
#    300 < lead <= 400
#          drug >  400
#
