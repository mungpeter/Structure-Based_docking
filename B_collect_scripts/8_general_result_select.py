#!/usr/bin/env python3

import sys
usage = '''\n\n  ## Usage: {0}
      -list < >  [ List of Mol ID to be extracted ]   
      -file <+>  [ SDF file(s): sdf ]   * gzip|xz|bzip2 accepted
      -pref < >  [ Output Prefix ]\n
  Optional:
      -tag  < >  [ SDF tag of MOL identifier (Def: 'Name') ]
      -sort < >  [ Sort based on SDF Tag, or "Name,Rank, Score" (Def: None) ]\n'''.format(sys.argv[0])
if len(sys.argv) < 4: sys.exit(usage)

import re,gc
import bz2,gzip,lzma
import pandas as pd

from rdkit import Chem
from rdkit.Chem import PandasTools as rdpd

from argparse import ArgumentParser

##########################################################################
def main():

  args = UserInput()

###############
  ## Read in the list of selected ligand ID 
  n_df     = pd.read_csv(args.mol_id, delimiter='\s+',header=None,comment='#',skip_blank_lines=True)
  keywords = n_df.loc[:, 0].to_list()
  print('\n > Number of items in <{}>: {}\n'.format( args.mol_id, len(keywords) ))
  print('first 5: ')
  print(keywords[:5])

  ## Extract the selected ligands from the supplied SDFs
  mol_sele = []
  for infile in args.infiles:
    df = RDkitRead(infile, removeHs=False)
    print(len(df))
    Items = df['ID'].apply(CheckID)
    df['Name']  = list(zip(*Items))[0]
    df['Rank']  = list(map(int, list(zip(*Items))[1]))
    df['Score'] = list(zip(*Items))[2]
    df['Soft']  = list(zip(*Items))[3]
    mol_sele.append( df[ df[args.id_tag].isin(keywords) ] )
    del df
    gc.collect()

#  all_df = pd.concat(mol_sele).set_index('ID').reindex(keywords).reset_index()
  all_df = pd.concat(mol_sele).reset_index(drop=True)
  print(all_df[:5])
  found_id  = all_df[args.id_tag].to_list()
  missed_id = [x for x in keywords if x not in set(found_id)]

  if missed_id is False:
    print('\033[31m  Info: \033[35m{0}\033[31m MOL cannot be found:\033[0m'.format(len(missed_id)))
    print(missed_id)

  ## Sort data, if needed
  if args.sort_tag:
    if re.search('Rank|Name', args.sort_tag, re.IGNORECASE):  asc = True
    else: asc = False
    all_df.sort_values(by=[args.sort_tag], ascending=asc, inplace=True)

  rdpd.WriteSDF(all_df, args.outpref+'.sdf.gz', molColName='mol', properties=list(all_df.columns))
  print('\033[31m  Info: \033[35m{0}\033[31m MOL output\033[0m'.format(len(all_df)))


##########################################################################
##
def CheckID( id ):
  if re.search(r'::', id):
    name, rank, score, soft = id.split('::')
    return [name, rank, score, soft]
  elif re.search(r' ', id):
    return [id.split()[0], 0, 0.0, '']
  else:
    return [id, 0, 0.0, '']

##########################################################################
## Read in SMILES or SDF input and add Hs to it
def RDkitRead( in_file, removeHs=True, add_Hs=False ):
  ## Read in SDF file; can choose to add hydrogens or not
  if re.search(r'.sdf', in_file):
    print(' # Reading SDF')
    df = rdpd.LoadSDF(  file_handle(in_file), removeHs=removeHs,
                        idName='ID', molColName='mol' )
    df['smiles'] = df.mol.apply(lambda m:Chem.MolToSmiles(Chem.RemoveHs(m)))
    if add_Hs:
      df['mol'] = df.mol.apply(Chem.AddHs)

  ## Read in SMILES file, check if there is a header "smiles"
  if re.search(r'.smi', in_file):
    print('# Reading SMI')
    with file_handle(in_file) as fi:
      if re.search('smi', str(fi.readline()), re.IGNORECASE):
        print('# Smiles input has Header #\n')
        df = pd.read_csv(in_file, sep='\s+', comment='#').dropna()
        df.columns = ['smiles','ID']
      else:
        print('# Smiles input has NO Header #\n')
        df = pd.read_csv(in_file, header=None, sep='\s+', comment='#').dropna()
        df.columns = ['smiles','ID']
    rdpd.AddMoleculeColumnToFrame(df, smilesCol='smiles', molCol='mol')
    df['smiles'] = df.mol.apply(Chem.MolToSmiles)

  ## Add Murcko scaffolds
  rdpd.AddMurckoToFrame(df, molCol='mol', MurckoCol='Hetero_Murcko',Generic=False)
  rdpd.AddMurckoToFrame(df, molCol='mol', MurckoCol='Generic_Murcko',Generic=True)

  print('## Number of MOL read from {}: {}\n'.format(in_file,len(df.smiles)))
  return df

#########################################################################
## Handle gzip and bzip2 file if the extension is right. otherwise, just open
## outuput: file handle
def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    handle = gzip.open(file_name, 'rb')
  if re.search(r'.xz$', file_name):
    handle = lzma.open(file_name, 'rb')
  elif re.search(r'.bz2$', file_name):
    handle = bz2.BZ2File(file_name, 'rb')
  else:
    if re.search(r'.smi', file_name):
      handle = open(file_name, 'r')
    else:
      handle = file_name
  return handle


##########################################################################
def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-list', dest='mol_id', required=True,
                  help='List of Mol_ID to be extracted')
  p.add_argument('-file', dest='infiles', required=True, nargs='+',
                  help='SDF File(s) * gzip,xz,bzip2 accepted')
  p.add_argument('-pref', dest='outpref', required=True,
                  help='Output prefix')

  p.add_argument('-tag', dest='id_tag', required=False, default='Name',
                  help='Optional: SDF Tag to be used to select MOL (Def: Name)')
  p.add_argument('-sort', dest='sort_tag', required=False, 
                  help='Optional: Sorting based on SDF Tag, or "Name, Rank, Score" (Def: None)')

  return p.parse_args()


##########################################################################
if __name__ == '__main__':
  main()


##########################################################################
#
#  Peter M.U. Ung @ gRED
#
#  v1   20.09.21
#  v2   21.11.29  add .xz compression capability
#  v3   21.12.24  add Murcko scaffold to molecules
  
#  *New Version based on RDkit PandasTools*
#
#  Extract and sort the selected molecules by ranking/score
#
##########################################################################
