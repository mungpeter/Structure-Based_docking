#!/usr/bin/env python3

##########################################################################
#
#	Peter M.U. Ung @ MSSM
#
#	v1.0	14.05.10
#       v1.1    14.06.04 -- change Score files input format
#       v3.0    14.09.10 -- correct 'average' sdf output bug
#       v5.0    14.09.18 -- reqirement to reach top X % in Y number of
#                           model to be considered in the consensus result
#
#	Purpose: perform a consensus of results from a set of dockings to 
#	         a single receptor. A molecule has to be in the top <X>% of
#                <Y> number of model to be considered in the consensus result.
#
#		 For each molecule, the highest ranking score and pose is 
#		 saved. 
#		 The final consensus is the ranking of the single best pose.
#
##########################################################################

import sys,os

msg = '''\n    > {0}
        [List of Score Files]   ** SDF and Score files must have same prefix
        [Score Column in Score File]  ** Standard numbering
        [Output Prefix] 
        [Threshold Top Ranking to Count: int %]
        [Threshold Number of model for Consensus Result: int]
        [Running mode]\n
      ## [Running mode]: One of the following options:
      # Analyze ALL models
        [-a   single best rank]
      # Minimal number of model to pass the threhold number and % ranking
        [-b   average of rank | -c single best rank]\n
    e.g. <script> score_file.list 2 consensus_result 50  2 -c\n
    Purpose: Generate consensus of docking based on 
          --Averge of rankings | Single Best Pose--\n'''.format(sys.argv[0])
if len(sys.argv) != 7: sys.exit(msg)

import re,glob
import gzip,bz2
import decimal
import operator
import numpy as np
import pandas as pd

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem

##########################################################################
def main(score_list, score_column, prefix, top_percent, consens_thresh, mode):
  Score_Files = []

  df = pd.read_csv(score_list, comment='#', header=None, delimiter='\s+')
  Score_Files = df.dropna().loc[:,0].to_numpy()
#  with open(score_list, 'r') as f:
#    Score_Files = [line.rstrip() for line in f if line.strip() 
#                    and not re.search(r'^#', line)]

  score_best = False
  conss_best = False
  if   mode == '-a': score_best = True
  elif mode == '-b': score_best = False
  elif mode == '-c': conss_best = True
  else: sys.exit('Mode of ranking not specified: '+mode)

  # Create a dict of scores from all imported scores
  Scores = {}
  for scores in Score_Files:
    print('  Reading: {0}'.format(scores))
    ReadScores(Scores, score_column, scores)
  print('\n  Number of Score File: {0}\n'.format(len(Score_Files)))

  # Convert the dict of scores into a ranked list of consensus scores.
  # Write the consensus scores to file
  Consensus = ConsensusScores(Scores, score_best, conss_best, prefix,
                              top_percent, consens_thresh)
  
  # Extract the docking pose of molecules according to it top consensus score.
  # Write the docking poses to file
  ExtractMol(Consensus, score_best, prefix)


##########################################################################
## Create a dict of scores from all imported scores
def ReadScores(Dict, column, score_file):
  Mol  = []
  prev_dict = len(Dict)

  Tmp = pd.read_csv(score_file, sep='\s+',comment='#',header=None).values.tolist()
#  Tmp = np.genfromtxt(score_file, comments='#', skip_header=True,
#                dtype={'formats': ('S20', np.float16),'names': ('Title', 'Score')})
  Mol = [[l[0], l[column-1]] for l in Tmp if l is not None]
  print('  Found entry: {0}\n'.format(len(Mol)))

  for idx, Items in enumerate(sorted(Mol, key=lambda tup: tup[1])):
    # If the Mol exist in Dict, add a new list of score and filename to it,
    # if not, create a new entry for it
    # Dict[name] = [score, score_file, rank]
    if Items[0] in Dict:
      Dict[Items[0]].append([Items[1], score_file, idx])
    else:
      Dict[Items[0]] = [[Items[1], score_file, idx]]

  print('  Accumulated number of molecule found: {0}'.format(len(Dict)))
  print('  New entry: {0}\n'.format(len(Dict) - prev_dict) )


##########################################################################
## Convert the dict of scores into a ranked list of consensus scores. 
## Write the consensus scores to file
def ConsensusScores(Dict, score_best, conss_best, prefix,
                    top_percent, consens_thresh):

  Consensus = []
  for mol_id in Dict:
    # Conditions to calculate "Best Score/Average Rank" 
    # for ALL models, or only those with X models in top Y% of the rank
    if score_best is True:
      Mol    = sorted(Dict[mol_id], key=lambda tup: tup[0])
    else:
      accept, Passed = ConsensusThreshold( Dict[mol_id], len(Dict),
                                            top_percent, consens_thresh )
      if accept is True:
        Mol = sorted(Passed, key=lambda tup: tup[0])
      else:
        continue

    # Use the best scoring pose for consensus, 
    # else, use the sum of rankings for consensus
    Ranks = [x[2] for x in Mol]
    try:
      log = np.log10(np.sum(Ranks))
    except RuntimeWarning:
      log = 'n/a'
#      continue
    if score_best is True or conss_best is True:
            #(best_score, mol_id, mol_in_file, log(Rank_Sum))
      Consensus.append([Mol[0][0], mol_id, Mol[0][1], log, len(Mol)])
    else:
      Consensus.append([log, mol_id, Mol[0][1], Mol[0][0], len(Mol)])
  Consensus.sort(key=lambda tup: tup[0])    # ranked by top ligands

  out = open(prefix+'.txt', 'w')

  if score_best is True or conss_best is True:
    out.write('#{0}\t{1}\t{2}\t{3}\t{4}\n'.format('mol_id','score','file',
                                                  'B-Rank','Passed'))
  else:
    out.write('#{0}\t{1}\t{2}\t{3}\t{4}\n'.format('mol_id','rank','file',
                                                  'B-Score','Passed'))

  print('''\n  Number of Molecule satisfying Consensus Threshold:
  <{0} models @ {1} %> out of {2} total:\n  {3}\n'''.format(consens_thresh,
                                    top_percent, len(Dict), len(Consensus)))

  for Item in Consensus:
    print(Item)
    try:
      float(Item[3])
    except ValueError:
      Item[3] = 0.0
    out.write( '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(
                Item[1],(Item[0]),Item[2],Item[3],Item[4]) )
  
  # Print frequency of model with top scoring ligands
  for Item in TopModelFrequency(Consensus, top_percent):
    out.write('# {0}\t{1}\n'.format(Item[0], Item[1]))
  out.close()

  return Consensus


##########################################################################
# Check if the molecule has X number of models that is in top Y% of the list
def ConsensusThreshold(Mol, lig_num, top_percent, consens_thresh):

  # Get the ranking cutoff for the top% required
  lig_cutoff = round(decimal.Decimal(lig_num*top_percent/100), 0)

  # Get how many docking result pass the ranking threshold 
  pass_threshold = 0
  Passed = []
  for Items in Mol:
    rank = Items[2]
    if rank <= lig_cutoff:
      pass_threshold += 1
      Passed.append(Items)

  # Proceed if the molecule has enough docking that passes the threshold
  if pass_threshold >= consens_thresh:
    return True, Passed
  else:
    return False, 0


##########################################################################
# Write out the frequency of model with top scoring ligands
def TopModelFrequency(List, top_percent):
  File = {}
  top = int(round(len(List)*(top_percent/100)))
  for Item in List[0:top]:
    if Item[2] in File:
      File[Item[2]] += 1
    else:
      File[Item[2]] = 1
  return sorted(File, key=lambda tup: tup[1], reverse=True)


##########################################################################
## Extract the docking pose of molecules according to it top consensus score.
## Write the docking poses to file
def ExtractMol(List, score_best, prefix):
  ConsScore = {}
  Saved_Mol = []

  # cluster molecules based on which SDF file they belong to
  for Mol in sorted(List, key=lambda tup: tup[2]):
    if Mol[2] in ConsScore:
      ConsScore[Mol[2]].append(Mol)
    else:
      ConsScore[Mol[2]] = [Mol]

  # from each SDF file, extract the docked pose
  for file_id in tqdm(ConsScore, total=len(ConsScore)):
    file_prefix = file_id.split('txt')[0]
    SDF         = glob.glob(file_prefix+'sdf*')
    if len(SDF) == 0: sys.exit('{0} or related SD file not found.'.format(file_prefix+'sdf*'))
    else: sdf_file = SDF[0]

    handle = file_handle(sdf_file)
    Temp = [x for x in Chem.ForwardSDMolSupplier(handle ,removeHs=False)
            if x is not None]
    SDMol = {}
    for mol in Temp:
      name = mol.GetProp('_Name')
#      name = mol.GetProp('_Name').split()[0]   # if name is separated
      if re.search(r':', name): # when the SD file is processed from docking
        SDMol[name.split(':')[0]] = mol
      else:
        SDMol[name] = mol

    for Mol in ConsScore[file_id]:
      try:
        test = SDMol[Mol[1]]
      except KeyError:
        print('{0} is not registered in database. Skip.'.format(Mol[1]))
        continue
#      if score_best is True: 
#        Saved_Mol.append([Mol[0], Mol[1], SDMol[Mol[1]]])
#      else:
      Saved_Mol.append([Mol[0], Mol[1], SDMol[Mol[1]], Mol[3]])
    del Temp
    del SDMol

#############

  # Sort all mol based on score and write out
  saved_sdf = Chem.SDWriter(prefix+'.sdf')
  for M in sorted(Saved_Mol, key=lambda tup: tup[0]):
    saved_sdf.write(M[2])
  saved_sdf.flush()
  saved_sdf.close()


#######################################################################
## to handle the raw SDF format, python3's rdkit has a documented bug and
## hasn't been fixed since 2016. https://github.com/rdkit/rdkit/issues/1065
## To avoid it, the input file cannot be an object handle of a regular file,
## i.e. handle = open('xxx.sdf','r') will fail but handle = 'xxx.sdf' is fine.
## It only happens to regular file but not to gzip.open() or bz2.BZ2File() in
## python3 rdkit but not in python2 rdkit...
## Fix it by replace handle = open('xxx.sdf','r') with handle = 'xxx.sdf'

## Handle gzip and bzip2 file if the extension is right. otherwise, just open
## outuput: file handle
def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    handle = gzip.open(file_name, 'rb')
  elif re.search(r'.bz2$', file_name):
    handle = bz2.BZ2File(file_name, 'rb')
  else:
    handle = open(file_name, 'rb')
  return handle

##########################################################################

if __name__ == '__main__':
  main( sys.argv[1], int(sys.argv[2]), sys.argv[3],
        float(sys.argv[4]), int(sys.argv[5]), sys.argv[6])
