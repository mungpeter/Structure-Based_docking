# Structure-Based_docking
Scripts to handle structure-based molecular docking, single-receptor virtual screening, consensus results of multiple receptor docking, and etc.

There are 2 major folders with primary running scripts (A_ and B_), 1 folder with miscellaneous scripts (Z_), and 1 folder with example scripts to run those primary scripts (Examples) :
```
  /A_docking_scripts   --   Setting up and running of docking
  /B_collect_scripts   --   Collecting the docking results
  /Z_misc_scripts      --   collection of less used scripts to manage docking results
  /Examples
      |------ /1_prep_ligands             - prepare ligand library
      |------ /2_virtual_screen           - running docking
      |------ /3_single_VS_rslt           - collect docking result for 1 structure
      |------ /4_consensus_VS_rslt        - consensus of multiple docking results
      |------ /5_extract_VS_data          - extract a subset of ligands from library
      |------ /6_print_select_2d          - print out docking result in 2D
      |------ /7_search_analogs_substruct - substructure search with SMILES patterns
      |------ /8_search_analogs_FP        - search SMILES patterns with fingerprints
```
#######################################################################################
# Docking

#######################################################################################
# Collecting results
```
> B_collect_scripts/5_general_screen_get_top.py
      -score    [ score file(s) ]  # accept multiple files "*.txt,x.txt"
      -sdf      [ sdf files(s)  ]  # same prefix as files in -score
      -top      [ number of top mol to output: int ]
      -dock     [ docking software: fred | sch | etc ]
      -outpref  [ output prefix for sdf, png, txt files ]

  Optional:
      -hmax     [ histogram max score: Def: < fred: -14 | sch: -10 > ]
      -hmin     [ histogram min score: Def: < fred:  -2 | sch:  -3 > ]
      -coll     [ collect X times top MOL in memory (Def: 2) ]
      -exclude  [ <SMARTS pattern> (smt-clean) ] removal filter
      -select   [ <SMARTS pattern> (smt-selec) ] selection filter

  e.g.>  *.py  -score "*docked.txt.bz2"   -sdf "*docked.sdf.bz2" 
               -top 1000    -dock sch     -outpref single_VS.summary
               -hmax "-16.0" -hmin "-2.0"     # need double quote to worl
               -exclude 'C(=O)[O-]|S(=O)(=O)[O-]|P(=O)(O)[O-]'  # no acidic moieties
               
 result > single_VS.summary.sch_top1k.sdf, single_VS.summary.sch_top1k.txt, 
          single_VS.summary.sch_top1k.histo.png
```
- Collect the top scoring results from one or multiple docking files and generate a histogram to show the distribution of the docking scores (.png). The score files (.txt) and the docking pose files (.sdf) can be in zipped format but must have the same prefix, i.e.:
> test.sch_docked.txt.bz2
> test.sch_docked.sdf.bz2


```
> B_collect_scripts/7_general_docking_cluster.py
             [ Docking input(s): sdf ]     - accept multiple files "*.sdf,x.sdf"
             [ Tanimoto cutoff: float ]  - 0.4 works for lead-/drug-like molecules
             [ output prefix ] 
             [ receptor structure: pdb ]
     Optional:
             [ -top   <no. of mol>: int ]
             [ -dl    DayLight Fingerprint      
               -ec    ECFP_4 Fingerprint (Default)
               -ms    MACCS 166-bit Key ]

             # also output a table of clustered molecules in PDF

 e.g.>  *.py vina.sdf,fred.sdf,ehits.sdf 
             0.4   output.clust-04   recpt.pdb.bz2   -dl

 result >  output.clust-04.sdf,    output.clust-04.pse,
           output.clust-04.pdf
```
- Cluster the supplied molecules using a fingerprint-based chemical similarity method with a user-defined Tanimoto cutoff value. A value of 0.40 works for most cases of lead-/drug-like molecules (MW 300-500+), while fragment-like may need a lower value, e.g. 0.37. The clustered results are placed into a PyMOL session file and the docking receptor is used. The clustered molecules are also output into a PDF file.

```
> B_collect_scripts/9_consensus_best_pose.py
        [ List of Score Files ]   ** SDF and Score files must have same prefix
        [ Score Column in Score File ]  ** Standard numbering
        [ Output Prefix ] 
        [ Threshold Top Ranking to Count: int % ]
        [ Threshold Number of model for Consensus Result: int ]
        [ Running mode ]
        
      ## [ Running mode ]: One of the following options, -c is recommended
      # Analyze ALL models
        [-a   single best rank]
      # Minimal number of model to pass the threhold number and % ranking
        [-b   average of rank | -c single best rank]
        
    e.g.>  *.py   score_file.list   2   consensus_result.cons_50-2   50    2    -c
 
 result >  consensus_result.cons_50-2.sdf,   consensus_result.cons_50-2.txt
 ```
- Perform a consensus of results from a set of dockings to a single receptor. A molecule has to be in the top <X>% of <Y> number of model to be considered in the consensus result. For each molecule, the highest ranking score and pose is saved. The final consensus is the ranking of the single best pose or the average of rank, tho single best rank is better overall. SDF and Score files used must have same prefix.

#######################################################################################
#######################################################################################
#######################################################################################

#######################################################################################
# Required software / packages
```
OpenEye           # 2018+
  fred
  omega2
  flipper
  fixpka
  make_receptor
```
```
Schrödinger       # 2019-01+
  ligprep
  glide
  glide_sort
  prepwizard
```
```
csh/tcsh          # shell
python            # 3.6.8+
  numpy           # 1.16.2+
  pandas          # 0.24.2+
  matplotlib      # 3.0.3+
  seaborn         # 0.9.0+
  rdkit           # 2019.03.01
  tabulate        # 0.8.3+
  pathos          # 0.2.3+
  tqdm            # 4.31.1+
  argparse        # 1.1
  html            # 
  gzip            # 
  bz2             # 
```
