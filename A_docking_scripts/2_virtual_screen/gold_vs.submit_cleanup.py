#!/usr/bin/env python3

import sys,os
import re

from rdkit import Chem
from rdkit.Chem import PandasTools as rdpd

import pandas as pd
from tqdm import tqdm
from pathos import multiprocessing
from argparse import ArgumentParser

msg = '''\n  > {0}
    -conf < >   [ GOLD setting file ]
    -rec  < >   [ Receptor structure file (.mol2) ]
    -cav  < >   [ Binding site cavity atoms file (.cavity.atoms) ] 
    -lig  < >   [ Ligand file (.sdf ; gzip/bz2 okay) ]
    -pref < >   [ Result prefix ]\n
  Optional:
    -genconf    [ Generate GOLD setting file: templ.conf.setting ]
    -nosort     [ No sorting of result ]
    -constr < > [ Constraint file for docking ]
    -top  < >   [ Save only the top <n> ligands; forced sorting ]
    -rlig < >   [ Reference ligand file (.mol2) ]\n'''.format(sys.argv[0])
if len(sys.argv) == 1: sys.exit(msg)

cwd = os.getcwd()

##########################################################################
def main():
  args = UserInput()
  if args.genconf:
    GenerateConfTemplInput()
    sys.exit()
  if args.savetop:
    try:
      savetop = int(args.savetop)
      args.nosort = True            # force sorting
    except TypeError:
      sys.exit('\033[31m  ERROR: -top must be an integer: \033[0m'+args.savetop)

########################
  ## Read input configure file
  settings = ReadConfSettings(args.conffile)
  settings['receptor']  = args.receptor
  settings['cavity']    = args.cavity
  settings['rslt_pref'] = args.rslt_pref

  if args.lig_ref: 
    settings['lig_ref'] = args.lig_ref
  if args.constr:
    settings['constr_file']  = args.constr

  ## handle sdf file in gzip/bzip2
  if   re.search('.gz$', args.ligand):
    ligand = args.ligand.split('/')[-1].split('.gz')[0]
    os.system('gunzip -c {0} > ./{1}'.format(args.ligand, ligand))
    settings['ligand'] = '{0}/{1}'.format(cwd, ligand)
  elif re.search('.bz2$', args.ligand):
    ligand = args.ligand.split('/')[-1].split('.bz2')[0]
    os.system('bunzip2 -c {0} > ./{1}'.format(args.ligand, ligand))
    settings['ligand'] = '{0}/{1}'.format(cwd, ligand)
  else:
    settings['ligand'] = args.ligand

  ## Write a list of gold.conf files
  Confs = GenerateConfFiles(settings)
  print('\033[34m## Generated subjobs: \033[33m{0}\033[0m'.format(len(Confs)))

  ## Run GOLD in parallel until all finished
  if int(args.cpu) > 0 and int(args.cpu) <= multiprocessing.cpu_count():
    core = int(args.cpu)
  else:
    core = multiprocessing.cpu_count()
  mpi = multiprocessing.Pool(core)
  tmp = [x for x in tqdm(mpi.imap(RunGOLD, Confs), total=len(Confs))]
  mpi.close()
  mpi.join()

############ Post-processing #############

  tmpdsf = settings['tmpdsf']
  findsf = settings['findsf']
  finssf = settings['finssf']
  mol_id = settings['mol_id']
  if settings['gold_funct'] == 'plp':
    score = 'Gold.PLP.Fitness'

  ## Modify each subjob docking result, summarize them all into 1 dataframe
  pref_list = [c.split('.conf')[0] for c in Confs]
  dock_list = []
  for pref in pref_list:
    os.chdir(pref)
    in_sdf  = '{0}.{1}'.format(pref, tmpdsf)
    out_sdf = '{0}.{1}'.format(pref, findsf)

    ## modify docked sdf file, collect them
    dock_list.append( RescaleRename( in_sdf, out_sdf, mol_id, score ) )
    os.system('bzip2 *sdf *lst')
    os.chdir(cwd)

  ## combine all subjob data, sort by ranking, output docked sdf and rank
  ## save only top ligands if needed
  xdf = pd.concat(dock_list)
  if not args.nosort:
    xdf.sort_values(by=[score], ascending=True, inplace=True)
    if args.savetop:
      xdf = xdf[:savetop]

  fin_sdf = '{0}.{1}'.format(settings['rslt_pref'], findsf)
  fin_scr = '{0}.{1}'.format(settings['rslt_pref'], finssf)
  rdpd.WriteSDF( xdf, fin_sdf, properties=list(xdf.columns) )
  xdf.to_csv( fin_scr, index=False, sep='\t', columns=[mol_id, score],
              header=False, float_format='%.3f' )
  os.system('bzip2 {0} {1}'.format(fin_sdf, fin_scr))


##########################################################################
##########################################################################

def GenerateConfFiles( settings ):

  st = settings

  step  = int( int(st['lig_mol_num']) / int(st['njobs']) )
  start = 1
  last  = step
  rslt_pref_orig = st['rslt_pref']

  cfiles = []
  ## for each njobs, write out a standalone gold.conf file
  for i in range(1, int(st['njobs'])+1, 1):

    ## incorporate job-specific variables
    st['lig_start'] = ' start_at_ligand  {0} '.format(start)
    st['lig_last']  = ' finish_at_ligand {0} '.format(last)
    st['rslt_pref'] = '{0}.{1}'.format(rslt_pref_orig, i)
    st['rslt_dock'] = '{0}.{1}.{2}'.format(rslt_pref_orig, i, st['tmpdsf'])

    ## special case: don't write out last ligand number in last njob
    if i == int(st['njobs']):
      st['lig_last'] = ' '

    ## fresh set of template data to generate subjob gold.conf
    keywords = SubstituteKeywords()

    cfile = []
    for line in GOLDconfTemplate():
      if not re.search('X(\w+)X', line):  # no modif if no variable 
        cfile.append(line)
        continue

      for x in keywords.keys():
        if re.search(keywords[x], line):
          if not st[x]:    # ignore line if variable is ''
            continue
          tmp = re.sub(keywords[x], st[x], line)
          if re.search('X(\w+)X', tmp):  # if multiple vars in line
            line = tmp
          else:
            cfile.append(tmp)
            break

    start = start + step
    last  = last  + step

    ## write out .conf file
    with open('{0}.conf'.format(st['rslt_pref']), 'w') as fo:
      for l in cfile:
        fo.write(l+'\n')
    cfiles.append('{0}.conf'.format(st['rslt_pref']))

  st['rslt_pref'] = rslt_pref_orig
  return cfiles


##########################################################################
def RunGOLD( conf ):
  os.system('gold_auto {0}'.format(conf))


## Rescale docking scores, rename MOL title to pre-GOLD title, output 
## modified SDF, keep dataframe in memory for later use
def RescaleRename( in_sdf, out_sdf, mol_id, score ):
  df   = rdpd.LoadSDF(in_sdf, removeHs=False, molColName='ROMol').fillna('')
  s_id = df[mol_id]

  def gold_scale( num ):       ## Scale the score by -0.1x
    return -(float(num)/10)

  for idx, row in df.iterrows():
    df['ROMol'][idx].SetProp('_Name', s_id[idx])
  df[score] = df[score].apply(gold_scale)

  return df


###########################################################################
def GOLDSettings():
  return {
    'scale': [0.67, '# Scaling of GA population (std = 1.0)\n'],
    'grid_origin': ['0,0,0', ' # Grid sphere center coordinates'],
    'grid_radius': [10, '# Grid sphere radius'],
    'lig_ref': ['""', '# Reference ligand; "" if none\n'],
    'gold_funct': ['plp','# CHEMPLP: plp (default); \n'],
    'ga_run': [10, '# GA run per ligand (def: 10)'],
    'mol_id': ['ID', '# SD Tag of mol ID\n'],
    'njobs':  [20, '# Number of subjobs'],
    'lig_mol_num': [500, '# mol no. in Library (dividable by subjobs)\n'],
    'constr_file': ['""', '# Text file with constraint terms; "" if none\n'],
    'tmpdsf': ['dock.sdf', '# Suffix of temporary docked SDF result'],
    'findsf': ['gold_docked.sdf', '# Suffix of final docked SDF result'],
    'finssf': ['gold_docked.txt', '# Suffix of final docked SDF ranking'],
  }


def SubstituteKeywords():
  return {
    'scale': 'XSCALEX',  'cavity': 'XCAVITYX',
    'grid_origin': 'XGRID_ORIGINX',  'grid_radius': 'XGRID_RADIUSX',
    'lig_ref': 'XLIG_REFX', 'gold_funct': 'XGOLD_FUNCTX',
    'receptor': 'XRECEPTORX', 'ligand': 'XLIGX', 'ga_run': 'XGA_RUNX',
    'lig_start': 'XLIG_STARTX', 'lig_last': 'XLIG_LASTX', 
    'rslt_pref': 'XRSLT_PREFX', 'rslt_dock': 'XRSLT_DOCKX',
    'constr_file': 'XCONSTRAINTSX'
  }


def GenerateConfTemplInput():
  settings = GOLDSettings()
  with open('templ.conf.setting', 'w') as f:
    for i in settings.keys():
      f.write('{0}\t{1}\t{2}\n'.format(i, settings[i][0], settings[i][1]))
  sys.exit('\033[34m## Generated - \033[32mtempl.conf.setting\033[0m')


def ReadConfSettings( infile ):
  temps = GOLDSettings()
  settings = {}
  for tmp in temps.keys():
    settings[tmp] = temps[tmp][0]   # drop the comment item

  s_df = pd.read_csv(infile, comment='#', sep='\s+', header=None).fillna('')
  for i in s_df.to_numpy():
    settings[i[0]] = re.sub(',',' ',i[1])
  return settings


###########################################################################
def UserInput():
  p = ArgumentParser()
  p.add_argument('-conf', dest='conffile', required=False,
                  help='GOLD setting file')
  p.add_argument('-cpu', dest='cpu', required=False,
                  help='CPU number; "0" to use all cores')

  p.add_argument('-rec', dest='receptor', required=False,
                  help='Receptor structure file (.mol2)')
  p.add_argument('-cav', dest='cavity', required=False,
                  help='Binding site cavity atoms file (.cavity.atoms)')
  p.add_argument('-lig', dest='ligand', required=False,
                  help='Ligand file (.sdf ; gzip okay)')
  p.add_argument('-pref', dest='rslt_pref', required=False,
                  help='Result prefix')

  p.add_argument('-genconf', dest='genconf', action='store_true',
                  help='Generate GOLD setting file: templ.conf.setting')
  p.add_argument('-nosort', dest='nosort', action='store_true',
                  help='No sorting of result')
  p.add_argument('-constr', dest='constr', required=False,
                  help='Constraint file (optional)')
  p.add_argument('-top', dest='savetop', required=False,
                  help='Save only the top <n> ligands; forced sorting (int ; optional)')
  p.add_argument('-rlig', dest='lig_ref', required=False,
                  help='Reference Ligand file (.mol2 ; optional)')

  return p.parse_args()


###########################################################################
## standard gold.conf file with tags (X(\w+)X) for substitution
def GOLDconfTemplate():
  conf = '''
#  Peter M.U. Ung @ gRED
# 
#  v1    20.08.10
#
#  Template GOLD config file for MPI runs
#

  GOLD CONFIGURATION FILE

  VARIABLE SETTINGS
autoscale = XSCALEX
origin = XGRID_ORIGINX
radius = XGRID_RADIUSX
cavity_file = XCAVITYX
ligand_reference_file = XLIG_REFX
gold_fitfunc_path = XGOLD_FUNCTX

protein_datafile = XRECEPTORX
ligand_data_file XLIGX XGA_RUNX XLIG_STARTX XLIG_LASTX

directory = XRSLT_PREFX
concatenated_output = XRSLT_DOCKX

  CONSTRAINTS
XCONSTRAINTSX
### Examples ###
###force_constraints = 1
###constraint sphere                  51.50 -26.63 35.55    2.250 1.0000 hydrophobic_atoms
###aconstraint pharmacophore D_A GAUSS 47.97 -27.65 29.29 1 3.5000 5.0000 1.0000
###constraint pharmacophore D_A GAUSS 49.77 -21.18 37.29 1 3.6000 5.0000 1.0000


#####################################################################
  POPULATION
popsiz          = auto
select_pressure = auto
n_islands       = auto
maxops          = auto
niche_siz       = auto

  GENETIC OPERATORS
pt_crosswt      = auto
allele_mutatewt = auto
migratewt       = auto

################################
  FLOOD FILL
do_cavity         = 1
floodfill_atom_no = 0
floodfill_center  = file

################################
  DATA FILES
param_file             = DEFAULT
set_ligand_atom_types  = 1
set_protein_atom_types = 0
tordist_file           = DEFAULT
make_subdirs           = 0
save_lone_pairs        = 0
fit_points_file        = fit_pts.mol2
read_fitpts            = 0

  FITNESS FUNCTION SETTINGS
initial_virtual_pt_match_max = 3
relative_ligand_energy       = 1
score_param_file             = DEFAULT

################################
  FLAGS
internal_ligand_h_bonds = 0
flip_free_corners       = 0
match_ring_templates    = 0
flip_amide_bonds        = 0
flip_planar_n           = 1 flip_ring_NRR flip_ring_NHR
flip_pyramidal_n        = 1
rotate_carboxylic_oh    = flip
use_tordist             = 1
postprocess_bonds       = 1
rotatable_bond_override_file = DEFAULT
solvate_all             = 1

  TERMINATION
early_termination = 0
n_top_solutions   = 3
rms_tolerance     = 1.5

  COVALENT BONDING
covalent = 0

#################################
  SAVE OPTIONS
save_score_in_file    = 1
save_protein_torsions = 0
output_file_format    = MACCS
bestranking_list_filename = bestrank.lst        # docked scores

clean_up_option delete_all_solutions
clean_up_option save_top_n_solutions 1
clean_up_option delete_redundant_log_files
clean_up_option delete_all_initialised_ligands
clean_up_option delete_empty_directories
clean_up_option delete_rank_file
clean_up_option delete_all_log_files
'''
  return list(filter(None,conf.split('\n')))


##########################################################################
if __name__ == '__main__':
  main()

##########################################################################
#
#  Peter M.U. Ung @ gRED
#
#  v1.0    20.08.10
#
#  Parallelize GOLD run for 1 ligand library and summarize the results into
#  format compariable to PU(me)-standardized FRED/Glide output.
#
#  The molecule title is converted back to the original name listed in user-
#  defined SD Tag, e.g. <ID>.
#  GOLD ChemPLP score is standardized by a factor of -0.1x to be in-line
#  with other scoring functions, e.g. FRED/Glide/Vina
#
#  rdkit      # 2020.03.1+
#  pandas     # 1.0.3+
#  GOLD       # 2020.1
#
###########################################################################
