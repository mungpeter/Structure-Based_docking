#!/usr/bin/env python3

import sys
from rdkit import Chem
from rdkit.Chem import PandasTools as rdpd
from argparse import ArgumentParser

msg = '''  > {0}
    -in  < >   [ Input SDF File  ] * gzip ok *
    -out < >   [ Output SDF File ] * gzip ok *
    -id  < >   [ SDF Tag with ID ]\n
  Optional:
    -score < > [ Gold Score (def: plp) ]
    -nosort    [ No sorting of result by -score ]\n'''.format(sys.argv[0])
if len(sys.argv) <= 4: sys.exit(msg)

##########################################################################
def main():
  args = UserInput()
  if not args.score:
    score = 'Gold.PLP.Fitness'
  else:
    if re.search('plp', args.score, re.IGNORECASE):
      score = 'Gold.PLP.Fitness'

  df = rdpd.LoadSDF(args.infile, removeHs=False, molColName='ROMol').fillna('')

  df[score] = df[score].apply(gold_scale)

  ## sort by default
  if not args.nosort:
    df.sort_values(by=[score], ascending=True, inplace=True)
    df['Rank'] = df.reset_index(drop=True).index
  else:
    df['Rank'] = -1
  
  for idx, row in df.iterrows():
    df['ROMol'][idx].SetProp('_Name', '{0}::{1}::{2:.1f}::{3}'.format(
                             row[args.mol_id],row['Rank']+1,row[score],'gold'))

  csv_suf = args.outfile.split('.sdf')[0]
  rdpd.WriteSDF(df, args.outfile, properties=list(df.columns))
  df.to_csv(csv_suf+'.txt', header=False, index=False, sep='\t',
            columns=[args.mol_id, score], float_format='%.3f')  


###########################################################################
## Scale the score by -0.1x
def gold_scale( num ):
  return -float(num)/10

###########################################################################
def UserInput():
  p = ArgumentParser()
  p.add_argument('-in', dest='infile', required=True,
                  help='Input SDF file')
  p.add_argument('-out', dest='outfile', required=True,
                  help='Output SDF file')
  p.add_argument('-id', dest='mol_id', required=True,
                  help='SDF Tag with ID of mols')
  p.add_argument('-score', dest='score', required=False,
                  help='SDF Tag with GOLD Score (def: Gold.PLP.Fitness')
  p.add_argument('-nosort', dest='nosort', action='store_true',
                  help='No sorting of docking result by -score')
  return p.parse_args()


###########################################################################
if __name__ == '__main__':
  main()

##########################################################################
#
#  Peter M.U. Ung @ gRED
#
#  v1.0    20.08.10

#  Specifically used for gRED MOE-fastROCS search result -- 
#  1. replace MolID in SDF MOL (Title) by removing tag appended by GOLD with
#  separator '|'
#  2. change 'Gold.PLP.Fitness' score in opposite direction and scaled by 0.1x
#  to match other docking scores (fred, glide, etc)
#
