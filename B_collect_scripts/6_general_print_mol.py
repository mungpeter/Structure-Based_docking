#!/usr/bin/env python3

import sys

msg = '''
    > {0} [Chemical Structure File: sdf|smi + gz|bz2]
'''.format(sys.argv[0])
if len(sys.argv) != 2: sys.exit(msg)

import glob,re,gzip,bz2

from rdkit import Chem
from rdkit_grid_print import grid_print

################################################################################

def main( filename ):

  mol_file = glob.glob(filename)[0]
  print('\n > File read: {}\n'.format(mol_file))

  if re.search(r'.sdf', mol_file):
    handle = file_handle(mol_file)

    Mol = [x for x in Chem.ForwardSDMolSupplier(handle, removeHs=True) 
            if x is not None]
    grid_print(mol_file.split('.sdf')[0], Mol, 'sdf')


  if re.search(r'.smi', mol_file):
    if re.search(r'.bz2$|.gz$', mol_file):
      print('\n  ## INFO: RDKit cannot take SMILES in zipped format, only ASCII\n')
    else:
      with open(mol_file, 'r') as fi:
        first_line = fi.readline()

      if re.search(r'smiles', first_line, re.IGNORECASE):
        Mol = [x for x in Chem.SmilesMolSupplier(mol_file, titleLine=True,
                 delimiter=' |\t|,') if x is not None]
      else:
        Mol = [x for x in Chem.SmilesMolSupplier(mol_file, titleLine=False,
                 delimiter=' |\t|,') if x is not None]

      grid_print(mol_file.split('.smi')[0], Mol, 'smi')



#################################################################################
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
      handle = open(file_name, 'rb')

  return handle


###############################################################################
if __name__ == '__main__':
  main( sys.argv[1] )

