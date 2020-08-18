#!/usr/bin/env python3

import sys
MSG = '''\n  ## Usage: x.py 
             [score files: txt]
             [Number of Top MOL in output: int]
             [docking software: fred | sch | etc]
             [Prefix of Output png, and txt files]\n
             [optional: -hmax=< default fred:-15.0 | sch:-10.0 >: float]
             [optional: -hmin=< default fred: -2.0 | sch:-3.0  >: float]\n
         ##  TXT and SDF files can also be in GZip/BZip2 format\n
         e.g.: x.py "*_score.txt"
               1000 sch ksr-allost -hmax=-16.0 -hmin=-2.0\n'''
if len(sys.argv) < 5 or len(sys.argv) > 7: sys.exit(MSG)

import glob, re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm

#### Default boundary constant for Histogram and changes ####
for argv in sys.argv:
  if re.search('sch',sys.argv[3]):
    upper, lower = -10., 0.
  else:
    upper, lower = -15., -2.  
  if re.search(r'-hmax=', argv): upper = float(argv.split('=')[1])
  if re.search(r'-hmin=', argv): lower = float(argv.split('=')[1])

##########################################################################
def doit( all_txt, all_top, dock, prefix ):

    # Read in .fred_score.txt
    File_Names = glob.glob(all_txt)
    print("Score File: ")
    print(File_Names)

    ## format the output name based on number of top output
    if all_top >= 1000:
      top_name = '{0}.{1}_top{2}k'.format( prefix, dock, all_top/1000 )
    else:
      top_name = '{0}.{1}_top{2}'.format( prefix, dock, all_top )

    all_df = ExtractScoreInfo( File_Names )
    ## Make histogram of ditribution of FRED scores
    Histogram( all_df.Score, all_top, top_name, dock, upper, lower )
    print("\n  ## Finished plotting overall Top-Ranks ##\n {0} / {1}\n\n".
            format(upper, lower))


#######################################################################
def ExtractScoreInfo( File_Names ):
  ## From the .fred_score.txt, extract the scores for ranking
#  data = np.genfromtxt(File_Names[0], comments=['#','Title', 'Name'],
#          dtype={'formats': ('S20', float),'names': ('Title', 'Score')})

#  d_df = pd.DataFrame(data, columns=['Title','Score']).sort_values(by=['Score'])
  d_df = pd.read_csv(File_Names[0], header=None, comment='#', sep='\s+')
  d_df.columns = ['Title', 'Score']
  
#  print(d_df[:20])
  x_df = d_df.loc[:5, ['Title', 'Score'] ]
  y_df = x_df[['Score','Title']]
  print(y_df)
  print(y_df.values)
  #print(x_df.set_index('Title').to_dict('list')['Score'])
  #print(d_df.loc[ d_df['Title'] == 'NCGC00187954-01' ] )
  print('# Total Ligand Collected: {0}'.format(len(d_df)))
  return d_df

#######################################################################
## Plot Histogram of Score distribution
def Histogram(all_df, top, top_name, dock, UPPER, LOWER ):

    bin_size  = 0.2
    text_high = 0.275
    text_hori = 0.4
    
#    fig, ax = plt.subplots()
    plt.figure(num=1, figsize=(8,6))

#    sns.set_style(style='white')
    sns.set_context(  context='notebook', font_scale=1.4)
#          rc={'font.sans-serif': 'Arial', 'font.size':14 })
    
    ## plot histogram - kernel density estimation
    fig = sns.kdeplot(all_df, shade=True, bw='scott')
    sns.despine()

    fig.set_xlabel('Score', fontsize=20)
    fig.set_ylabel('Fraction of Docked Molecules', fontsize=20 )
    fig.set_title( top_name+": "+str(len(all_df)), fontsize=20 )
    fig.set_xlim([UPPER, LOWER])
    fig.legend().remove()

    ## Draw a vertical line to indicate the Top hits
    print(all_df[top-1])
    fig.axvline( x=all_df[top-1], ymin=0, ymax=1000, 
                  color='r', linewidth=3 )
    top_num = 'Top {0}: {1:.2f}'.format(top, all_df[top-1])
    fig.text( all_df[top-1]-text_hori, text_high, 
              top_num, rotation=90, color='black', fontsize=18 )

    ## Draw a vertical line to indicate the Median Score
    fig.axvline( x=np.median(all_df), ymin=0, ymax=1000, 
                  color='b', linewidth=3 )
    median = 'Median:{0:.2f}'.format(np.median(all_df))
    fig.text( np.median(all_df)-text_hori, text_high, 
              median, rotation=90, color='black', fontsize=16 )

    ## Draw 2 vertical lines to indicate the standard deviation
    fig.axvline( x=(np.median(all_df)+np.std(all_df)), 
                  ymin=0, ymax=1000, color='k', linewidth=1.5 )
    fig.axvline( x=(np.median(all_df)-np.std(all_df)), 
                  ymin=0, ymax=1000, color='k', linewidth=1.5 )
    stdev = 'StDev: {0:.2f}'.format(np.std(all_df))
    fig.text( np.median(all_df)+np.std(all_df)-text_hori+.1, text_high, 
              stdev, rotation=90, color='k', fontsize=14 )

    plt.show()
    fig.figure.savefig( top_name+'.histo.png', dpi=300,  )
    #fig.close()


############################################################################
if __name__ == '__main__':
    doit( sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4] )
