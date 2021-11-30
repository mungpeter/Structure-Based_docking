#!/usr/bin/env python3

import sys
from rdkit import Chem
from rdkit.Chem import PandasTools as rdpd
from argparse import ArgumentParser

msg = '''\n  > {0}
    -in    < >  [ Input hitlist SDF File  ] * gzip ok *
    -pref  < >  [ Output prefix ]\n
  Optional:
    -top   < >  [ Output the top <num> of mol (def: -1) ]
    -score < >  [ SDF Tag with Score (def: Chemgauss4) ]
    -dock  < >  [ Docking method fred|hybrid (def: fred) ]\n'''.format(sys.argv[0])
if len(sys.argv) < 3: sys.exit(msg)

##########################################################################
def main():
  args = UserInput()
  score = args.score
  dock = args.dock
  top = int(args.top)

  df = rdpd.LoadSDF(args.infile, removeHs=False, molColName='ROMol',
                    idName='ID')[:top].fillna('')
  print('\033[34m> select mol: \033[32m{0}\033[0m'.format(len(df)))
  df[score]  = df[score].apply(float)
  df['Rank'] = df.index

  for idx, row in df.iterrows():
    df['ROMol'][idx].SetProp('_Name', '{0}::{1}::{2:.2f}::{3}'.format(
                                row.ID, row['Rank']+1, row[score], dock))

  sdf_out = '{0}.{1}_docked.sdf'.format(args.outpref, dock)
  csv_out = '{0}.{1}_docked.txt.bz2'.format(args.outpref, dock)

  rdpd.WriteSDF(df, sdf_out, properties=list(df.columns))
  df.to_csv(csv_out, header=False, index=False, sep='\t',
            columns=['ID', score], float_format='%.3f')


###########################################################################
def UserInput():
  p = ArgumentParser()
  p.add_argument('-in', dest='infile', required=True,
                  help='Input hitlist SDF file *gzip okay*')
  p.add_argument('-pref', dest='outpref', required=True,
                  help='Output prefix')

  p.add_argument('-top', dest='top', required=False, default=-1,
                  help='Output top <num> mol (def: -1)')
  p.add_argument('-score', dest='score', required=False, default='Chemgauss4',
                  help='SDF Tag with Score (def: Chemgauss4)')
  p.add_argument('-dock', dest='dock', required=False, default='fred',
                  help='Docking method fred|hybrid (def: fred)')

  return p.parse_args()


###########################################################################
if __name__ == '__main__':
  main()

##########################################################################
#
#  Peter M.U. Ung @ gRED
#
#  v1.0    20.08.16
#
#  Convert OE_Orion-generated hitlist.sdf into common format
#  OE_Orion output a sorted hitlist.sdf file (typically 10,000 results) 
#  with <Chemgauss4> tag but without <mol id> tag, so need to extract 
#  mol_id from the title.
#
