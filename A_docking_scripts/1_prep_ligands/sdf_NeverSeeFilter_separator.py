#!/usr/bin/env python3

import sys,os

msg = '''  > {0}
    -in   <>  [ Input  SDF/smi file                        ] ** gzip/bz2 okay **
    -pass <>  [ Output SDF/smi passed gRED NeverSee filter ] ** gzip/bz2 okay **
    -nsee <>  [ Output SDF/smi failed gRED NeverSee filter ] ** gzip/bz2 okay **\n'''.format(sys.argv[0])
if len(sys.argv) != 7: sys.exit(msg)

import re
import gzip,bz2
import pandas as pd
from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import PandasTools as rdpd
rdBase.DisableLog('rdApp.*')

from argparse import ArgumentParser


##########################################################################
def main():
  args = UserInput()

  df = RDkitRead( args.infile, args.id, removeHs=False, add_Hs=False )

  nsee_df = df[df['NeverSee_Groups'] == 'Y'] ; len(nsee_df)
  pass_df = df[df['NeverSee_Groups'] == 'N'] ; len(pass_df)

  print('\033[34m Passed NeverSee Filter: \033[32m{0}\033[0m'.format(len(pass_df)))
  print('\033[34m Failed NeverSee Filter: \033[31m{0}\033[0m'.format(len(nsee_df)))

  if re.search(r'.smi', args.nsee_file, re.IGNORECASE):
    nsee_df.smiles = nsee_df.MOL.apply(lambda m: Chem.MolToSmiles(Chem.RemoveHs(m)))
    nsee_df.to_csv(args.nsee_file, columns=['smiles', 'ID'], sep=' ', header=False, index=False)
  else:
    rdpd.WriteSDF(nsee_df, args.nsee_file, molColName='MOL', properties=list(nsee_df.columns))

  if re.search(r'.smi', args.pass_file, re.IGNORECASE):
    pass_df.smiles = pass_df.MOL.apply(lambda m: Chem.MolToSmiles(Chem.RemoveHs(m)))
    pass_df.to_csv(args.pass_file, columns=['smiles', 'ID'], sep=' ', header=False, index=False)
  else:
    rdpd.WriteSDF(pass_df, args.pass_file, molColName='MOL', properties=list(pass_df.columns))
  print('')


##########################################################################
## Read in SMILES or SDF input and add Hs to it
def RDkitRead( in_file, idnm, removeHs=False, add_Hs=False ):

  ## Read in SDF file; can choose to add hydrogens or not
  if re.search(r'.sdf', in_file):
    print(' \033[34m## Reading SDF ##\033[0m')
    df = rdpd.LoadSDF( file_handle(in_file), removeHs=removeHs, smilesName='smiles', molColName='MOL' )
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


##########################################################################
## Handle gzip and bzip2 file if the extension is right. otherwise, just open
## outuput: file handle
def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    handle = gzip.open(file_name, 'rb')
  elif re.search(r'.bz2$', file_name):
    handle = bz2.open(file_name, 'rb')
  else:
    if re.search(r'.smi', file_name):
      handle = open(file_name, 'r')
    else:
      handle = file_name
  return handle


##########################################################################
def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-in', dest='infile', required=True,
                  help='Input mol file (sdf|smi) * gzip/bzip2 accepted')
  p.add_argument('-pass', dest='pass_file', required=True,
                  help='Write out mols passed NeverSee Filter (sdf|smi)')
  p.add_argument('-nsee', dest='nsee_file', required=True,
                  help='Write out mols failed NeverSee Filter (sdf|smi)')

## Place holder only
  p.add_argument('-id', dest='id', required=False,
                  help='Optional: Identifier used in input file (Def: ID)')

  args=p.parse_args()
  return args


##########################################################################
if __name__ == '__main__':
  main( )

##########################################################################
#
#  Peter M.U. Ung @ gRED
#
#  v1.   20.07.16
#
#  Parse a SDF file tagged by NeverSee Filter, separate those with 'Y' flag
#  from the rest and write out the passed/flagged molecules in SDF/smi
#
#  Can read/write in .gz/.bz2 format too
#
#  rdkit  # 2020.03.01+
