#!/usr/bin/env python3

import sys
import os,re
import gzip,bz2
import pandas as pd

from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import SaltRemover
from rdkit.Chem import MolStandardize
from rdkit.Chem import PandasTools as rdpd
from rdkit.Chem.MolStandardize import rdMolStandardize
rdBase.DisableLog('rdApp.*')

# rdkit.__version__   # 2019.09.03

from argparse import ArgumentParser

msg = '''
  > {0}
      -in  <file>        [ Input mol file: sdf|smi ]  * gzip/bzip2 accepted
      -fo  <format>      [ Output format:  sdf|smi ]
      -out <prefix>      [ Output mol file prefix  ]
      -id  <identifier>  [ Optional: Identifier used in file (Def: ID) ]
      -st  <start num>   [ Optional: Starting number of mol (Def: 0) ]
      -en  <end num>     [ Optional: Ending number of mol (Def: -1) ]
'''.format(sys.argv[0])
if len(sys.argv) == 1: sys.exit(msg)

##########################################################################

def main( ):

  args  = UserInput()
  if args.id is None:
    args.id = 'ID'
  if args.start is None:
    args.start = 0
  if args.end is None:
    args.end = -1

  df = RDkitRead(args.infile, args.id, removeHs=True, add_Hs=False)[int(args.start):int(args.end)].dropna()

  remover = SaltRemover.SaltRemover()
  normzer = rdMolStandardize.Normalizer()
  chooser = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)

  ## remove salts
  print('\033[34m## Desalting moleucles...\033[0m\n')
  df['mol'] = df.MOL.apply(remover.StripMol)
  ## choose largest fragment (most Hs)
  print('\033[34m## Choosing moleucles...\033[0m\n')
  df['mol2'] = df.mol.apply(chooser.choose)
  ## clean molecule (not really relevant?)
  print('\033[34m## Cleaning moleucles...\033[0m\n')
  df['mol3'] = df.mol2.apply(normzer.normalize)
  ## rewrite SMILES with newest mol3
  print('\033[34m## Converting moleucles...\033[0m\n')
  df['smiles'] = df.mol3.apply(Chem.MolToSmiles)

  if   args.format == 'sdf':
    rdpd.WriteSDF(df, args.outpref+'.'+args.format, molColName='mol3', idName=args.id, properties=['smiles'])
  elif args.format == 'smi':
    df.to_csv(args.outpref+'.'+args.format, index=False, sep=' ', columns=['smiles',args.id], header=True)



##########################################################################
## Read in SMILES or SDF input and add Hs to it
def RDkitRead( in_file, idnm, removeHs=False, add_Hs=False ):

  ## Read in SDF file; can choose to add hydrogens or not
  if re.search(r'.sdf', in_file):
    print(' \033[34m## Reading SDF ##\033[0m')
    df = rdpd.LoadSDF(  file_handle(in_file), removeHs=removeHs,
                        smilesName='smiles', molColName='MOL' )
    if add_Hs:
      df['MOL'] = df.mol.apply(Chem.AddHs)

  ## Read in SMILES file, check if there is a header "smiles"
  if re.search(r'.smi', in_file):
    print('  \033[34m## Reading SMI ##\033[0m')
    with file_handle(in_file) as fi:
      if re.search('smi', str(fi.readline()), re.IGNORECASE):
        print(' \033[36m# Smiles input has Header #\033[0m\n')
        df = pd.read_csv(in_file, sep='\s+').dropna()
        df.columns = ['smiles', idnm]
      else:
        print(' \033[35m# Smiles input has NO Header #\033[0m\n')
        df = pd.read_csv(in_file, header=None, sep='\s+', comment='#').dropna()
        df.columns = ['smiles', idnm]
    df['MOL'] = df.smiles.apply(Chem.MolFromSmiles)

  print('## Number of MOL read from \033[34m{0}: \033[31m{1}\033[0m\n'.format(in_file,len(df)))
  return df

#########################################################################
## Handle gzip and bzip2 file if the extension is right. otherwise, just open
## outuput: file handle
def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    handle = gzip.open(file_name, 'rb')
  elif re.search(r'.bz2$', file_name):
    handle = bz2.BZ2File(file_name, 'rb')
  else:
    if re.search(r'.smi', file_name):
      handle = open(file_name, 'r')
    else:
      handle = file_name
  return handle

#######################################################################
def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-in', dest='infile', required=True,
                  help='Input mol file (sdf|smi) * gzip/bzip2 accepted')
  p.add_argument('-fo', dest='format', required=True,
                  help='Output format (sdf|smi)')  
  p.add_argument('-out', dest='outpref', required=True,
                  help='Output mol file prefix')

  p.add_argument('-id', dest='id', required=False,
                  help='Optional: Identifier used in input file (Def: ID)')
  p.add_argument('-st', dest='start', required=False,
                  help='Optional: Starting number of mol (Def: 0) ')
  p.add_argument('-en', dest='end', required=False,
                  help='Optional: Ending number of mol (Def: -1) ')

  args=p.parse_args()
  return args

#######################################################################
if __name__ == '__main__':
  main( )


#######################################################################
#
#  Peter M.U. Ung @ MSSM/Yale
#
#  v1   20.03.26
#
#  Remove salts and only keep the largest fragment (by HA and H numbers)
#
#  MolStandardize available from 2019.09.01
#
#
######################################################################
