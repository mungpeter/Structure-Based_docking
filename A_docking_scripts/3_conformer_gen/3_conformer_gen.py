#!/usr/bin/env python3

##########################################################################
#
#   Peter M.U. Ung  @ MSSM
#
#   v1.0    15.05.15
#   v2.0    19.03.21 - simplify sdf output
#   v3.0    19.12.24 - use the new ETKDGv2 conformer generation method
#
#   Purpose:    Use RDKit's 3D conformer generator to make ligand conformers
#               similar to OpenEye's OMEGA
#               Howeever, still very time consuming and no where near
#               the speed of OEMGA, and enforce all-protonation with no
#               pKa control
#
#   based on:   conformers_pre2018_09.py   Class file to generate conformers
#               https://github.com/skearnes/rdkit-utils/blob/master/rdkit_utils/conformers.py
#
#               ETKDG method (2015)  enable skipping force field minimzation
#               https://pubs.acs.org/doi/abs/10.1021/acs.jcim.5b00654
#               
##########################################################################
import sys
msg = '''
    Usage: {0}
            -in        Ligand File: sdf,smi   (gzip/bzip2 accepted)
            -out       Output Prefix          (output SDF.gz and log file)\n
          --Options
            -addh      Add Hydrogens when use SDF file (Def: False)
            -cpu       CPU for MPI: int       (Def: all available)
            -maxconfs  Fixed number of Max. conformers to generate
                        (default: adaptive)
            -ff        Minimization with UFF Forcefield (Def: false)
            -rmsd      RMSD cutoff for conformers to group together
                        (Def: 0.50)\n
    e.g.:  {0}
            -in lig.sdf.bz2 -out output -cpu 8 -rmsd 0.65 -addh -ff\n
'''.format(sys.argv[0])
if len(sys.argv) > 13 or len(sys.argv) <= 4: sys.exit(msg)

import re,os
import gzip,bz2
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import AllChem
from rdkit.Chem import MolStandardize
from rdkit.Chem.rdMolAlign import AlignMol

from rdkit.Chem import PandasTools as rdpd
print('**', 'RDKit:', rdBase.rdkitVersion)

from tqdm import tqdm
from pathos import multiprocessing
from argparse import ArgumentParser

##########################################################################
def main():
  args = UserInput()
  if not args.cpu: cpu = multiprocessing.cpu_count()
  else:            cpu = int(args.cpu)

  if not args.maxconfs: maxconfs = 0
  else:                 maxconfs = int(args.maxconfs)

  if not args.rmsd:     rmsd = 0.50
  else:                 rmsd = float(args.rmsd)

  mpi = multiprocessing.Pool(processes=cpu)
  cCF = ETKDG_ConfGenerator( numConfs=maxconfs, rmsd=rmsd, run_uff=args.ff )

  m_df = RDkitRead( args.in_file, keep_Hs=True, add_Hs=args.add_Hs )

  ## Generate conformer 3D coordinate vectors for each input molecule
  print('  ## Generateing conformer 3D coordinates ##')
  mol_in = list(zip(m_df.mol.to_numpy(), m_df.ID.to_numpy()))
  Confs = [x for x in tqdm( mpi.imap(cCF, mol_in), total=len(m_df) )]
#  Confs = [cCF(m) for m in mol_in]
  mpi.close()
  mpi.join()

  ## Convert conformer 3D coordinate vectors to actual coordinates
  ## !! cannot do rdkit tagging in mpi mode, have to do it serially
#  Mols = [x for x in tqdm( mpi.map(ConvertConformerResult, Confs), 
#                            total=len(Confs))]
  Mols = [ConvertConformerResult(conf) for conf in Confs]

  print('# Print results of mols: {}\n'.format(len(Mols)))
  fo = Chem.SDWriter(args.outpref+'.sdf')
  with open(args.outpref+'.conf.log', 'w') as log:
    for inp in Mols:
      mol, name, conf_num = inp
      log.write('{}\t{}\n'.format(name, conf_num))

      for conf in mol:
        fo.write(conf)
  fo.close()

  os.system('bzip2 -f {}.sdf'.format(args.outpref))


##########################################################################
##########################################################################
class ETKDG_ConfGenerator(object):
  def __init__( self, numConfs='', run_uff='', rmsd='' ):
    self.numConfs = numConfs    # max. allowed conformation
    self.run_uff  = run_uff     # on/off to run UFF minimization
    self.rmsd     = rmsd        # rmsd threshold as different cluster

  def __call__( self, inp ):
    return self.conf_generator( inp )

#########################
  # adaptive max conformer to generate based on rotatable bond count
  def adaptive_conf( self, mol ):
    '''
        default conformer number based on no. rotatable bond
          rb >= 13       max_conformers = 300
          10 < rb <= 13  max_conformers = 250
          8  < rb <= 10  max_conformers = 200
          5  < rb <= 8   max_conformers = 150
          3  < rb <= 5   max_conformers = 100
          rb <= 3        max_conformers = 50
        override these settings if numConfs != 0
    '''
    # enforce user input max conformer number if -maxconfs is > 0
    if self.numConfs != 0:
      max_conf = self.numConfs
    else:
      rb_num = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
      if  rb_num <= 3:
        max_conf = 50
      elif rb_num > 3  and rb_num <= 5:
        max_conf = 100
      elif rb_num > 5  and rb_num <= 8:
        max_conf = 150
      elif rb_num > 8  and rb_num <= 10:
        max_conf = 175
      elif rb_num > 10 and rb_num <= 13:
        max_conf = 200
      else:
        max_conf = 250

    return max_conf

############################
  ## Calculate energies of conformers stored in mol
  def get_conformer_energies( self, mol, ids ):
    energies = []
    for _id in ids:
      ff = AllChem.UFFGetMoleculeForceField(mol, confId=_id)
      energies.append( ff.CalcEnergy() )
    energies = np.asarray(energies, dtype=np.float32)
    return energies

############################
  ## Prune conformers from a molecule using an RMSD threshold, starting
  ## with the lowest energy conformer. Used when UFF minimization is ON
  def prune_conformers( self, mol, ids, max_conf ):
    if self.rmsd <= 0. or mol.GetNumConformers() <= 1:
      return mol

    energies = self.get_conformer_energies(mol, ids)
    rmsd = get_conformer_rmsd(mol)

    sort = np.argsort(energies)  # sort by increasing energy
    keep = []  # always keep lowest-energy conformer
    discard = []
    for i in sort:
      # always keep lowest-energy conformer
      if len(keep) == 0:
        keep.append(i)
        continue

      # discard conformers after max_conformers is reached
      if len(keep) >= max_conf:
        discard.append(i)
        continue

      # get RMSD to selected conformers
      this_rmsd = rmsd[i][np.asarray(keep, dtype=int)]

      # discard conformers within the RMSD threshold
      if np.all(this_rmsd >= self.rmsd):
        keep.append(i)
      else:
        discard.append(i)

  # create a new molecule to hold the chosen conformers
  # this ensures proper conformer IDs and energy-based ordering
    new = Chem.Mol(mol)
    new.RemoveAllConformers()
    conf_ids = [conf.GetId() for conf in mol.GetConformers()]
    for i in keep:
      conf = mol.GetConformer(conf_ids[i])
      new.AddConformer(conf, assignId=True)

    return new


############################
############################
  def conf_generator( self, inp ):
    mol, name = inp

    conf_parm = AllChem.ETKDGv2()
    conf_parm.pruneRmsThresh = self.rmsd
    conf_parm.randomSeed     = -1

    ## get rotatable bond-dependent adaptive conformation number
    max_conf = self.adaptive_conf(mol)

    ## Generate 3D conformers, map atom 3D vectors in 'ids' to 'mol'
    ## Hydrogens are supposed to be added beforehand
    ids = AllChem.EmbedMultipleConfs(mol, max_conf, conf_parm)

    ## align all conformers to 1st frame
    rmslist = []
    AllChem.AlignMolConformers(mol, RMSlist=rmslist)

    ## Minimize conformers with UFF, 2x slower than without
    ## with minimization, parameters can be used to cluster conformers
    if self.run_uff:
      for _id in ids: 
        AllChem.UFFOptimizeMolecule(mol, confId = _id)
      mol = self.prune_conformers(mol, ids, max_conf)

    return [mol, name]
    # the item 'ids' generated in EmbedMultipleConfs can be converted to
    # a list list(ids), but it is just a list of numbers corresponding to
    # the conformers stored in 'mol', so it is unnecessary, really


##########################################################################
## Calculate conformer-conformer pairwise RMSD.
def get_conformer_rmsd( mol ):
  rmsd = np.zeros((mol.GetNumConformers(), mol.GetNumConformers()),
                    dtype=np.float32)
  for i, ref_conf in enumerate(mol.GetConformers()):
    for j, fit_conf in enumerate(mol.GetConformers()):
      if i >= j:
        continue
      rmsd[i, j] = AlignMol(mol, mol, ref_conf.GetId(), fit_conf.GetId(),
                            maxIters=200)
      rmsd[j, i] = rmsd[i, j]

  return rmsd


##########################################################################
## Map the 3D coordinates in vectors back to the reference 2D molecule to
## generate the 3D structure. Append information to each molecule/conformer,
## write to a log file how many conformer actually generated
def ConvertConformerResult( inp ):
  cf, name = inp
  Mols = []

  for id in range(0, cf.GetNumConformers()):
    mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(cf,confId=id),removeHs=False)
    mol.SetProp('_Name', name+'_'+str(id+1))
    mol.SetProp('Group', name)
    mol.SetProp('Conformation', str(id+1))
    Mols.append(mol)

  return Mols, name, len(Mols)


##########################################################################
## Read in SMILES or SDF input and add Hs to it
def RDkitRead( in_file, keep_Hs=True, add_Hs=False ):
  ## Read in SDF file; can choose to add hydrogens or not
  if re.search(r'.sdf', in_file):
    df = rdpd.LoadSDF(  file_handle(in_file), removeHs=keep_Hs,
                        smilesName='smiles', molColName='mol' )
    if add_Hs:
      df['mol'] = df.mol.apply(Chem.AddHs)

  ## Read in SMILES file, check if there is a header "smiles"
  if re.search(r'.smi', in_file):
    with file_handle(in_file) as fi:
      if re.search('smi', str(fi.readline()), re.IGNORECASE):
        print('\n # Smiles input has Header #\n')
        df = pd.read_csv(in_file, sep='\s+', comment='#').dropna()
        df.columns = ['smiles','ID']
      else:
        print('\n # Smiles input has NO Header #\n')
        df = pd.read_csv(in_file, header=None, sep='\s+', comment='#').dropna()
        df.columns = ['smiles','ID']
    df['mol'] = df.smiles.apply(CleanSmilesForMol)
    
  print('## Number of MOL read from {}: {}\n'.format(in_file,len(df.smiles)))
  return df

##########################################################################
## Standardize and add Hs to SMILES input only
def CleanSmilesForMol( smi ):
  std_smi = MolStandardize.standardize_smiles(Chem.CanonSmiles(smi))
  tau_smi = MolStandardize.canonicalize_tautomer_smiles(std_smi)
  tau_mol = Chem.MolFromSmiles(tau_smi)
  mol_h   = Chem.AddHs(tau_mol)
  return mol_h

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


##########################################################################
def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-in', dest='in_file', required=True,
                  help='Input molecular file: sdf,smi (zipped okay)')
  p.add_argument('-out', dest='outpref', required=True,
                  help='Output prefix')

  p.add_argument('-addh', dest='add_Hs', action='store_true',
                  help='Enforce adding Hydrogens when using SDF input (def: False)')
  p.add_argument('-maxconfs', dest='maxconfs', required=False,
                  help='Maximum conformer to generate (def: 0 = adaptive)')
  p.add_argument('-rmsd', dest='rmsd', required=False,
                  help='RMSD threshold to group conformers as 1 cluster (def: 0.50)')
  p.add_argument('-ff', dest='ff', action='store_true',
                  help='Boolean flag to use UFF minimization (def: False)')
  p.add_argument('-cpu', dest='cpu', required=False,
                  help='MPI CPU number (def: 0 = all)')
          
  args=p.parse_args()
  return args


##########################################################################
if __name__ == '__main__':
  main()
