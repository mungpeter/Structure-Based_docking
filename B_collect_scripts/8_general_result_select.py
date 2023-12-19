#!/usr/bin/env python3

import sys
usage = '''\n\n  ## Usage: {0}
      -id   < >  [ List of Mol ID to be extracted ]   
      -sdf  <+>  [ SDF file(s): sdf ]   * gzip|xz|bzip2 accepted
      -pref < >  [ Output Prefix ]\n
  Optional:
      -tag  < >  [ SDF tag of MOL identifier (Def: 'Name') ]
      -sort < >  [ Sort based on SDF Tag, or "Name,Rank, Score" (Def: None) ]\n
    FRED/HYBRID score tag: Chemgauss4
    GlideSP score tag:     r_i_glide_gscore\n'''.format(sys.argv[0])
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
  n_df = pd.read_csv(args.mol_id,delimiter='\s+',header=None,comment='#',skip_blank_lines=True)

  ## ZINC data uses numeric molname, enforce string
  mol_id = n_df.loc[:, 0].astype(str).to_list()

  ## Check how many IDs are unique/duplicate
  uniq_id = list(set(mol_id))
  dupl_id = duplicate_check(mol_id)

  print('\n\033[32m > Number of items in <{0}>:\033[0m {1}'.format( args.mol_id, len(mol_id) ))
  if dupl_id:
    print('\033[32m > Unique IDs:\033[0m {0}'.format(len(uniq_id)))
    print('\033[32m > Duplicates:\033[0m')
    print(dupl_id)
  print('\n\033[32m > First 5 items: \033[0m')
  print(uniq_id[:5])

  ## Extract the selected ligands from the supplied SDFs
  mol_sele = []
  for infile in args.infiles:
    df = RDkitRead(infile, removeHs=False)
#    df.mol.apply(lambda x: Chem.AddHs(x,addCoords=True)).apply(MolStandardize.rdMolStandardize.Normalize) # temp fix for mol without Hs

    print(len(df))
    Items = df['Title'].apply(CheckTitle)
    df['Name']  = list(zip(*Items))[0]
    df['Rank']  = list(map(int, list(zip(*Items))[1]))
    df['Score'] = list(zip(*Items))[2]
    df['Soft']  = list(zip(*Items))[3]
    mol_sele.append( df[ df[args.id_tag].isin(uniq_id) ] )
    del df
    gc.collect()

#  all_df = pd.concat(mol_sele).set_index('Title').reindex(uniq_id).reset_index()
  all_df = pd.concat(mol_sele).drop_duplicates(subset=args.id_tag).reset_index(drop=True)
  print(all_df[:5])
  found_id  = all_df[args.id_tag].to_list()
  missed_id = [x for x in uniq_id if x not in set(found_id)]

  if missed_id:
    print('\033[31m  Error: \033[31m{0}\033[31m MOL cannot be found, double check the IDs:\033[0m'.format(len(missed_id)))
    print(missed_id)
  else:
    print('\033[32m  Info: All MOL in hitlist is accounted for in SDF file\033[0m')

  ## Sort data, if needed
  if args.sort_tag:
    if re.search('Rank|Name', args.sort_tag, re.IGNORECASE):  asc = True
    else: asc = False
    all_df.sort_values(by=[args.sort_tag], ascending=asc, inplace=True)

  rdpd.WriteSDF(all_df, args.outpref+'.sdf.gz', molColName='mol', properties=list(all_df.columns))

  print('\033[32m  Info: \033[0m{0}\033[32m MOL output\033[0m'.format(len(all_df)))


##########################################################################
## Check if the list has duplicated mol_id
def duplicate_check( inp ):
  seen = set()
  uniq = []
  dupl = []
  for x in inp:
    if x not in seen:
      uniq.append(x)
      seen.add(x)
    else:
      dupl.append(x)
  return dupl

##########################################################################
## Check if 'Title' of molblock is modified to contain addition data
def CheckTitle( title ):
  if re.search(r'::', title):
    name, rank, score, software = title.split('::')
    return [name, rank, score, software]
  elif re.search(r' ', title):
    return [title.split()[0], 0, 0.0, '']
  else:
    return [title, 0, 0.0, '']

##########################################################################
## Read in SMILES or SDF input and add Hs to it
def RDkitRead( in_file, removeHs=False, add_Hs=False ):
  ## Read in SDF file; can choose to add hydrogens or not
  if re.search(r'.sdf', in_file):
    print('\n\033[31m # Reading SDF #\033[0m')
    df = rdpd.LoadSDF(  file_handle(in_file), removeHs=removeHs,
                        idName='Title', molColName='mol' )
#    df['smiles'] = df.mol.apply(lambda m:Chem.MolToSmiles(Chem.RemoveHs(m)))
    if add_Hs:
      df['mol'] = df.mol.apply(Chem.AddHs)

    ## ZINC library molname can be numeric, issue with pandas/numpy handling
    df['Title'] = df['Title'].astype(str)

  ## Read in SMILES file, check if there is a header "smiles"
  if re.search(r'.smi', in_file):
    print('# Reading SMI')
    with file_handle(in_file) as fi:
      if re.search('smi', str(fi.readline()), re.IGNORECASE):
        print('# Smiles input has Title #\n')
        df = pd.read_csv(in_file, sep='\s+', comment='#').dropna()
        df.columns = ['smiles','Title']
      else:
        print('# Smiles input has NO Title #\n')
        df = pd.read_csv(in_file, header=None, sep='\s+', comment='#').dropna()
        df.columns = ['smiles','Title']
    rdpd.AddMoleculeColumnToFrame(df, smilesCol='smiles', molCol='mol')
    df['smiles'] = df.mol.apply(Chem.MolToSmiles)

  ## Add Murcko scaffolds
  rdpd.AddMurckoToFrame(df, molCol='mol', MurckoCol='Hetero_Murcko',Generic=False)
  rdpd.AddMurckoToFrame(df, molCol='mol', MurckoCol='Generic_Murcko',Generic=True)

  print('## Number of MOL read from {}: {}\n'.format(in_file,len(df)))
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

  p.add_argument('-id', dest='mol_id', required=True,
                  help='List of Mol_ID to be extracted')
  p.add_argument('-sdf', dest='infiles', required=True, nargs='+',
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
#  v4   22.09.06  handle ZINC numeric molname; rename Title name from 'ID' to 'Title'
#  v5   23.04.12  check for missing molecule or duplicated input
#
#  *New Version based on RDkit PandasTools*
#
#  Extract and sort the selected molecules by ranking/score
#
##########################################################################
