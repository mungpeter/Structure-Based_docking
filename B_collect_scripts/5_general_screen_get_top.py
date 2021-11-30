#!/usr/bin/env python3

#######################################################################
##
##	Peter M.U. Ung @ MSSM
##	
##  v1.0  - 13.11.13
##  v2.0  - 13.11.20 - add FRED tag to filename for future process
##	               FRED score is 1 decimal place
##  v3.0  - 13.11.25 - change to output number based on SDF, not score
##  v4.0  - 13.12.25 - change the read-write organization to use less 
##	               memory
##  v4.5  - 13.12.27 - fixed bug in sdf file reading added functions,
##	               read GZip and BZip2 files
##  v5.0  - 14.05.27 - change Histogram range, allow optional input
##  v6.0  - 16.12.21 - read sdf if name has "::"
##  v6.1  - 17.11.13 - chomp on SDF molname to avoid backspaces
##  v6.2  - 18.02.28 - set default upper/lower for FRED and Glide
##  v7.0  - 18.08.28 - enable SMARTS match to filter out substructures
##  v8.0  - 18.08.29 - rewrite
##  v8.1  - 18.10.30 - bugfix, SMARTS filters for selection and exclusion
##  v9.0  - 19.05.08 - use Seaborn to improve visual of histogram and add
##                     mpi to reading process, but need to watch out mem
##  v10.  - 19.10.15 - fixed a bug with int(args.all_top)*coll > len(d_df)
##  v11   - 20.09.11   use PandasTools, optional histogram generation
##  v12   - 20.12.20   remove sklearn processes, bugfix
##  v13   - 21.11.29   add xz compression capability
##
##	Take *_score.txt generated by OpenEye FRED docking to rank molecules.
##	Then read in corresponding SDFs to select ranked molecules for output.
##	Print out the top-ranking sdf molecules and generate a histogram.	
##    -select|-exclude option enables filtering of molecules with matching
##    substructure
##
##	Required:	fred_screen_preprocess.py
##			*.fred_score.txt
##			*.fred_docked.sdf(.bz2|.gz|.xz)
##
#######################################################################

import sys
MSG = """\n  ## Usage: x.py 
      -score   <+>  [ score files: txt ]    # gzip/xz/bzip2 okay
      -sdf     <+>  [ sdf files: sdf ]      # gzip/xz/bzip2 okay
      -top     < >  [ Number of Top MOL in output: int ]
      -dock    < >  [ docking software: fred | sch | etc ]
      -outpref < >  [ Prefix of Output sdf, png, and txt files ]\n
  Optional:
      -sdfname < >  [ SDF tag to match mol Name  to score file (def: original title) ]
      -sdfscor < >  [ SDF tag to match mol Score to score file (def: original title) ]
      -histo        [ Generate Histogram of score distribution (def: False) ]
      -png          [ Generate grid-image of SDF (def: False) ]
      -hmax    < >  [ < default fred:-14.0 | sch:-10.0 >: float ]
      -hmin    < >  [ < default fred: -2.0 | sch:-3.0  >: float ]\n
      -exclude < >  [ <SMARTS filter> (smt-clean) ]   removal filter
      -select  < >  [ <SMARTS filter> (smt-selec) ]   selection filter
                    [ use when SMARTS filtering is enabled ]\n
  e.g.: x.py -score "*_score.txt" -sdf "*.sdf" 
             -top 1000 -dock sch -outpref ksr-allost 
             -histo -hmax '-16.0' -hmin '-2.0'     # need double quote to work
             -exclude 'C(=O)[O-]|S(=O)(=O)[O-]|P(=O)(O)[O-]'  # acidic moieties\n"""
if len(sys.argv) == 1: sys.exit(MSG)

import glob,re,gc
import gzip,bz2,lzma
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools as rdpd

from pathos import multiprocessing
from argparse import ArgumentParser
from rdkit_grid_print import grid_print

############################################################################
def doit( ):

    args = UserInput()
    if args.dock == 'fred':
      if args.hmax is None:
        upper = -16.0
      else:
        upper = float(args.hmax)
      if args.hmin is None:
        lower = -2.0
      else:
        lower = float(args.hmin)
    else:
      if args.hmax is None:
        upper = -10.0
      else:
        upper = float(args.hmax)
      if args.hmin is None:
        lower = 0.0
      else:
        lower = float(args.hmin)

    if args.arg_exc is not None:
      arg_exc = args.arg_exc
    else:
      arg_exc = None
    if args.arg_sel is not None:
      arg_sel = args.arg_sel
    else:
      arg_sel = None

######################
    # Read in .fred_score.txt
    File_Names = list(glob.glob(args.all_txt))
    print("\033[31mScore Files:\033[0m")
    print(File_Names)

    ## format the output name based on number of top output
    if int(args.all_top) >= 1000:
      top_name = '{0}.{1}_top{2}k'.format( args.prefix, args.dock, int(int(args.all_top)/1000) )
    else:
      top_name = '{0}.{1}_top{2}'.format( args.prefix, args.dock, args.all_top )

    mpi  = multiprocessing.Pool()
    Data = [x for x in tqdm(mpi.imap(ExtractScore,File_Names),total=len(File_Names))]
#    Data = [ExtractScoreInfo(fn) for fn in File_Names]
    mpi.close()
    mpi.join()

    df = pd.concat(Data)
    d_df = df[df.columns[0:2]]
    d_df.columns = [args.sdfname, args.sdfscore]
    d_df.sort_values(by=[args.sdfscore], inplace=True)
    print('\033[31m## Entries read:\033[0m {0}\n'.format(len(df)))
    print(d_df[:5])

    ## Make histogram of ditribution of FRED scores
    if args.histo:
      Histogram( d_df.Score, int(args.all_top), top_name, args.dock, upper, lower )
      print("\n  ## Finished plotting overall Top-Ranks ##\n {0} / {1}\n\n".
            format(upper, lower))

##################
    # Read in SDF file name
    SDF_Names = glob.glob(args.all_sdf)

    ## Read in top SDF files and build ranked SDF file
    mpi = multiprocessing.Pool(processes=int(multiprocessing.cpu_count()/3))
    sdfdata = CollectSDFData(mol_df=d_df[:int(args.all_top)], name=args.sdfname)
#    Temp = [x for x in tqdm(mpi.imap(sdfdata, SDF_Names), total=len(SDF_Names))]
    Temp = [sdfdata(sdf) for sdf in SDF_Names]
    mpi.close()
    mpi.join()

    Top_sdf = pd.concat(Temp).sort_values(by=args.sdfscore).drop_duplicates(subset=args.sdfname,keep='first').reset_index(drop=True).reset_index()
    Top_sdf['NewID'] = ['{0}::{1}::{2:.2f}::{3}'.format(row[args.sdfname],i+1,row[args.sdfscore],args.dock) for i,row in Top_sdf.iterrows() ]
   
    if arg_exc is None and arg_sel is None:
      switch = None
      WriteSDFDataSelect( Top_sdf, args.arg_sel, top_name, args.sdfname, args.sdfscore, switch, args.grid )
    elif arg_exc is None and arg_sel is not None:
      switch = 'include'
      WriteSDFDataSelect( Top_sdf, args.arg_sel, top_name, args.sdfname, args.sdfscore, switch, args.grid )
    elif arg_exc is not None and arg_sel is None:
      switch = 'exlcude'
      WriteSDFDataSelect( Top_sdf, args.arg_exc, top_name, args.sdfname, args.sdfscore, switch, args.grid )
    else:
      sys.exit(' ## ERROR: Only 1 SMARTS filtering should be used ##\n')


#######################################################################
def ExtractScore( fname ):
  print('file_name: '+str(fname))
  df = pd.read_csv(fname, sep='\s+', comment='#')
  print('# Ligand Collected: {0}'.format(len(df)))
  return df


##########################################################################
## Build a database of molecules from SDF files
class CollectSDFData(object):
  def __init__( self, mol_df, name='ID' ):
    self.mol_df = mol_df
    self.name   = name

  def __call__( self, sdf_file ):
    return self._read_sdf(sdf_file)

  def _read_sdf( self, sdf_file ):
    ## Build a library of molecules found in the Top-Selction List
    df = rdpd.LoadSDF(sdf_file, molColName='ROMol', idName='ID', removeHs=False)
    print('  # SDF mol read in from > {0} = {1}'.format(sdf_file, len(df)))
    sel_df = df[ df[self.name].isin(self.mol_df[self.name]) ]
    top_df = pd.merge(sel_df, self.mol_df, on=self.name)
    del df
    gc.collect()    # active collection of memory to avoid crash
    return top_df


#######################################################################
# Remove molecules matching SMARTS strings into .smarts.* files until 
# reaching the targeted number of top-selected molecules
def WriteSDFDataSelect( Top_sdf, arg_pat, top_name, name, score, switch, grid ):

  if not switch:
    rdpd.WriteSDF(Top_sdf, '{0}.sdf.gz'.format(top_name), idName='NewID', properties=list(Top_sdf.columns))
    Top_sdf.to_csv('{0}.txt.gz'.format(top_name), columns=[name, 'index', score], sep='\t')
    if grid: grid_print(top_name, Top_sdf, 'sdf')
    return None
  else:

    match = []
    for smarts in [ p for p in arg_pat.split('|') ]:
      match.append(Top_sdf['ROMol'].apply(lambda x: x.HasSubstructMatch(Chem.MolFromSmarts(smarts))))

    Select = []
    for mol in zip(*match):
      if True in set(mol):
        Exlcude.append(True)
      else:
        Exclude.append(False)
    Top_sdf['match'] = Select

    if switch == 'include':
      sel_sdf = Top_sdf[Top_sdf['match'] == True ]
    elif switch == 'exclude':
      sel_sdf = Top_sdf[Top_sdf['match'] == False]

    rdpd.WriteSDF(sel_sdf, '{0}.{1}.sdf.gz'.format(top_name,switch), idName='NewID', properties=list(sel_sdf.columns))
    sel_sdf.to_csv('{0}.{1}.txt.gz'.format(top_name,switch), columns=[name, 'index', score], sep='\t')
    if grid: grid_print(top_name, sel_df, 'sdf')


#######################################################################
## if "ligand name" has '::' due to previous rdkit processing, remove the
## added data and just the "name" again
def OriginalNameSDF(sdfs):
  NewData = []
  for mol in sdfs:
    name = mol.GetProp('_Name')
    mol.SetProp('_Name', name.split('::')[0])
    NewData.append(mol)
  return NewData


#######################################################################
## Plot Histogram of Score distribution
def Histogram(Histo, top, top_name, dock, UPPER, LOWER ):

    ## if input top number is larger than database size:
    if len(Histo) < top:
      top = len(Histo)

    bin_size  = 0.2
    text_high = 0.275
    text_hori = 0.3

    plt.figure(num=1, figsize=(8,6))

    sns.set_context(  context='notebook', font_scale=1.4)
#          rc={'font.sans-serif': 'Arial', 'font.size':14 })
    
    ## plot histogram - kernel density estimation
    fig = sns.kdeplot(Histo, shade=True, bw='scott')
    sns.despine()

    fig.set_xlabel('Score', fontsize=20 )
    fig.set_ylabel('Fraction of Docked Molecules', fontsize=20 )
    fig.set_title( top_name+": "+str(len(Histo)), fontsize=20 )
    fig.set_xlim([UPPER, LOWER])
    fig.legend().remove()

    ## Draw a vertical line to indicate the Top hits
    fig.axvline( x=Histo[top-1], ymin=0, ymax=1000, 
                color='r', linewidth=3 )
    top_num = 'Top {0}: {1:.2f}'.format(top, Histo[top-1])
    fig.text( Histo[top-1]-text_hori, text_high, 
              top_num, rotation=90, color='black', fontsize=18 )

    ## Draw a vertical line to indicate the Median Score
    fig.axvline( x=np.median(Histo), ymin=0, ymax=1000, 
                color='b', linewidth=3 )
    median = 'Median:{0:.2f}'.format(np.median(Histo))
    fig.text( np.median(Histo)-text_hori, text_high, 
              median, rotation=90, color='black', fontsize=16 )

    ## Draw 2 vertical lines to indicate the standard deviation
    fig.axvline( x=(np.median(Histo)+np.std(Histo)), 
                ymin=0, ymax=1000, color='k', linewidth=1.5 )
    fig.axvline( x=(np.median(Histo)-np.std(Histo)), 
                ymin=0, ymax=1000, color='k', linewidth=1.5 )
    stdev = 'StDev: {0:.2f}'.format(np.std(Histo))
    fig.text( np.median(Histo)+np.std(Histo)-text_hori, text_high, 
              stdev, rotation=90, color='k', fontsize=16 )

#    plt.show()
    fig.figure.savefig( top_name+'.histo.png', dpi=300 )
    plt.close()


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
  elif re.search(r'.xz$', file_name):
    handle = lzma.open(file_name, 'rb')
  elif re.search(r'.bz2$', file_name):
    handle = bz2.BZ2File(file_name, 'rb')
  else:
    handle = open(file_name, 'rb')
  return handle

###########################################################################
#### Default boundary constant for Histogram and changes ####
grid= False

def UserInput():
  p = ArgumentParser(description='Command Line Arguments')

  p.add_argument('-score', dest='all_txt', required=True,
                  help='Score Files: txt')
  p.add_argument('-sdf', dest='all_sdf', required=True,
                  help='sdf files: sdf')
  p.add_argument('-top', dest='all_top', required=True,
                  help='Numner of Top mol in output: int')
  p.add_argument('-dock', dest='dock', required=True,
                  help='docking software: fred | sch | etc')
  p.add_argument('-outpref', dest='prefix', required=True,
                  help='Prefix of Output sdf, png, and txt files')

  p.add_argument('-sdfname', dest='sdfname', required=False, default='ID',
                  help='Optional: SDF Tag to match SDF name to score file (def: original title)')
  p.add_argument('-sdfscor', dest='sdfscore', required=False, default='Score',
                  help='Optional: SDF Tag to match SDF name to score file (def: original title)')

  p.add_argument('-histo', dest='histo', required=False, action='store_true',
                  help='Optional: Generate Histogram of score distribution (def: False)')
  p.add_argument('-hmax', dest='hmax', required=False,
                  help="optional: -hmax=< default fred:-14.0 | sch:-10.0 >: ''-float''")
  p.add_argument('-hmin', dest='hmin', required=False,
                  help="optional: -hmin=< default fred: -2.0 | sch:-3.0  >: ''-float''")
  p.add_argument('-png', dest='grid', required=False, action='store_true', help='?? ignore')

  p.add_argument('-exclude', dest='arg_exc', required=False,
                  help='[Optional: -exclude=<SMARTS filter> (smt-clean)] removal filter')
  p.add_argument('-select', dest='arg_sel', required=False,
                  help='[optional: -select=<SMARTS filter>  (smt-selec)] selection filter\n[             use when SMARTS filtering is enabled ]')

  args = p.parse_args()
  return args

############################################################################
if __name__ == '__main__':
    doit(  )
