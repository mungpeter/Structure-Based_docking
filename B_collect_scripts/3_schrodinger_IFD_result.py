#!/usr/bin/env python3

##########################################################################
#
#   Peter M.U. Ung @ MSSM
#
#   v1.0    16.08.03
#   v1.1    16.10.25    - bugfix
#   v2.0    17.03.02    read in molecule smi and rename molecule 
#   v3.0    17.12.20    use multiprocessing to convert maegz to pdb
#   v4.0    18.01.08    parse the result file for docking scoring
#                       simplify the folder management of the results
#   v5.0    19.12.25    remove dependence on separate .py scripts
#                       renamed from /glide/0_glide_ifd_mae2pdb.py
#
#   Consolidate Schrodinger's IFD top results and put them into a folder
#       * automatically extract the result structure files from csv 
#       * consolidate the result structures into pymol 
#
##########################################################################

import sys
msg = '''\n    > {0}
\t\t[ Ligand file: sdf ] *give full path to sdf file
\t\t[ directory of IFD result ]
\t\t[ output filename prefix ]\n
      e.g.> x.py atp.sdf InducedFit_1 type_IB.dyrk1a.hit1\n
      ***   start outside of the _ifd folder
            '''.format(sys.argv[0])
if len(sys.argv) != 4: sys.exit(msg)

import os,re
import shutil
import pandas as pd

from rdkit import Chem
from rdkit.Chem import PandasTools as rdpd

from tqdm import tqdm
from pathos import multiprocessing

ligand  = sys.argv[1]
inpref  = sys.argv[2]    # IFD result directory
outpref = sys.argv[3]    # prefix of results

##########################################################################
def file_handle(file_name):
  if re.search(r'.gz$', file_name):
    handle = gzip.open(file_name, 'rb')
  elif re.search(r'.bz2$', file_name):
    handle = bz2.BZ2File(file_name, 'rb')
  else:
    handle = file_name
  return handle

##########################################################################
## Going into IFD home directory
os.chdir(inpref)

## Extract molecule names and listing order from the original input ligands
Names = rdpd.LoadSDF( file_handle(ligand), molColName='mol').ID.tolist()

## Going into IFD scoring result directory
os.chdir(inpref+'_workdir/scoring_dir')

## create result folder
if not os.path.exists('1_data'):
  os.makedirs('1_data')

## Parse the result .csv to extract 'filenames', redirect results from
## glide_docking_dir_2 to scoring_dir
## ** sed "s/glide_docking_dir_2/scoring_dir/" report.csv > report_corr.csv **
os.system('sed "s/glide_docking_dir_2/scoring_dir/" report.csv > report_corr.csv')
os.system('cp report_corr.csv 1_data')

##########################################################################
## Parse the result .csv to extract 'filenames' for each ligand
# Ligands are identified by their order of appearance, as number in integer
df = pd.read_csv('report_corr.csv', delimiter=',')
Items = df.iloc[:, [1,5,2,3,4]].to_numpy().tolist()

#with open('report_corr.csv', 'r') as fi:
#  for line in fi:
#    if re.search(r'Entry', line): continue
#    Entries = line.split(',')
#    Items.append( [ int(Entries[1]),   Entries[5], 
#                    float(Entries[2]), float(Entries[3]), float(Entries[4]) ] )

## Extract the result structure in 'filenames'
LigGroups = {}
for item in Items:
  struct_name = item[1].rstrip()
  if not LigGroups.get(item[0]):
    LigGroups[item[0]] = [ [struct_name, item[2], item[3], item[4]] ]
  else:
    LigGroups[item[0]].append( [struct_name, item[2], item[3], item[4]] )
print('# Number of Ligand Group found: '+str(len(LigGroups)))

#for key in xrange(1,len(Structures)+1,1)):
keylist = sorted(LigGroups.keys())     # sort the 'key' by Ligand number
print(keylist)

## print out data into readable format
NewRslt = {}
Columns = ['Ligand','IFDScore','GlideScore','Prime_Energy']
for key in keylist:
  NewRslt[Names[int(key)-1]] = LigGroups[key]

#L = [(k, *t) for k, v in NewRslt.items() for t in v]            # python 3
L = [(k, t[1], t[2], t[3]) for k, v in NewRslt.items() for t in v]    # python 2
df = pd.DataFrame(L, columns=Columns)
df.to_csv('1_data/'+outpref+'.result.csv', header=True, index=False)
df.to_excel('1_data/'+outpref+'.result.xlsx', header=True, index=False)


##########################################################################
## convert maegz files to pdb
class Maegz2pdb( object ):

  def __init__( self, pref=None, group=None, names=None ):
    self.group = group
    self.names = names
    self.pref  = pref
  def __call__( self, keylist ):
    return self.run_convert(keylist)

########################
  def build_multi_pdb( self, tmp_pdb, key ):
    # Write to a bz2 multi-pdb file, with MODEL {n} as designation
    with bz2.open('{0}.{1}.pdb.bz2'.format(self.pref, key), 'wb') as fo:

      for idx, pdb in enumerate(tmp_pdb):
        fo.write('MODEL {0}\n'.format(idx+1))

        with open(pdb, 'r') as fi:
          for l in fi:
            if re.search('^MODEL|^ENDMDL', l): continue
            fo.write(l)
        fo.write('ENDMDL\n')
          
########################
  def run_convert( self, key ):
    print('key:  '+str(key))
    print('Name: '+self.names[int(key)-1])

    tmp_pdb = []
    for idx, mae in enumerate(self.group[key]):
      name = mae[0].split('.mae')[0]
      print(name)
      tmp_pdb.append('_temp.{0}.{1}.pdb'.format(key, idx))

      if not os.path.isfile('{0}.{1}.pdb'.format(self.pref, key)):
        os.system('$SCHRODINGER/utilities/structconvert -imae {0}.maegz -opdb _temp.{1}.{2}.pdb'.format(name, key, idx))

    self.build_multi_pdb(tmp_pdb, key)


##########################################################################
## Find the result files in the result csv and use them to build pymol session

## Multiprocess for converting maegz results to pdb format
m2p = Maegz2pdb(pref=inpref, group=LigGroups, names=Names)

mpi = multiprocessing.Pool()
tmp = [x for x in tqdm(mpi.imap(m2p, keylist), total=len(keylist))]
mpi.close()
mpi.join()

## correspond the order of result in .csv list to input molcules
pml = open('_temp.pml', 'w')
for key in keylist:
#  m2p(key)   # serial process, use when mulitprocess fails

  pml.write('load {0}.{1}.pdb.bz2, {2}.{1}\n'.format(inpref, key,
                                                     Names[int(key)-1]))
  pml.write('distance hb.{0}, {0}.{1} and poly, {0}.{1} and org, mode=2\n'.format(Names[int(key)-1], key))


##########################################################################
## write out settings to pml file for pymol session generation
pml.write('show cartoon\n')
pml.write('hide lines\nhide labels\nset valence\n')
pml.write('show sticks, org\nshow lines, org\n')
pml.write('show lines, byres poly within 6 of org\n')
pml.write('hide (h. and (e. c extend 1))\n')  #hide nonpolar hydrogens
pml.write('set light_count, 1\nset ray_trace_mode, 1\n')
pml.write('set ray_trace_gain, .008\nset ray_trace_color, black\n')
pml.write('set opaque_background, 0\nutil.cbaw poly\n')
pml.write('util.cbas org\n')
pml.write('center org\n')
pml.write('zoom org\n')
pml.write('set pse_export_version, 1.70\n')
pml.write('save {0}.pse\n'.format(inpref))
pml.close()
os.system('pymol -c _temp.pml')

##########################################################################
## manage result files
os.system('bzip2 {0}.pse'.format(inpref))
os.system('tar -cf {0}.pdb.tar {0}.*pdb.bz2'.format(inpref))
shutil.copy('{0}.pse.bz2'.format(inpref), '1_data')
shutil.copy('{0}.pdb.tar'.format(inpref), '1_data')

os.system('rm {0}*pdb* {0}.pse.bz2 _temp* '.format(inpref))

shutil.copytree('1_data', '../../1_data')
os.chdir('../..')

print('  ## Zipping up the IFD work directory ##\n')
os.system('tar -jcf {0}_workdir.tbz {0}_workdir'.format(inpref))


##########################################################################
