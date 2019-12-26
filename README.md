# Structure-Based_docking

Scripts to handle structure-based molecular docking, single-receptor virtual screening, consensus results of multiple receptor docking, and etc.

There are 2 major folders with _primary running_ scripts (**A_docking_scripts** and **B_collect_scripts**), 1 folder with _miscellaneous_ scripts (**Z_misc_scripts**), and 1 folder with _example_ scripts to run those primary scripts (**Examples**) :

```
  /A_docking_scripts   --   Setting up and running of docking
      |------ /1_prep_ligands             - scripts to prepare ligand for docking
      |------ /2_virtual_screen           - scripts to run docking
      |------ /3_conformer_gen            - generate 3D conformers
                    |------ /CSD_confgen  - CSD's generator

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
# Steps to Prepare Ligand Library for Docking - 19.07.28

**Prepare initial ligand library**

- 0. Download library subsets in SMILES format.
- 1. check total number of ligands, collect all SMILES from different subsets
   into one single file, the shuffle the SMILES to randomize the order
```
    > cat */*/*smi | grep -v 'smiles' > all.smi
    > sort -uR all.smi > shuffled.smi
    > cat shuffled.smi | wc -l
```
- 2. to the shuffled SMILES file, split it into a number of library subsets with each subset  has ~100,000 ligands (a '.' is needed before 'smi' and after 'new_name'). This _subset number_ will need to be precalculated by knowing how many ligands are there (from previous step, 'cat shuffled.smi | wc -l')
```
    > split --numeric-suffixes=1         \
            --additional-suffix=".smi"   \
            --number=30                  \
            shuffled.smi "new_name."
```

> one issue with these shuffled SMILES files -- the first line is somehow always missing some characters, rendering the first line (molecule) of all smiles files to be bad and any conversion will fail. But that only affect 75 molecules out of ~ 6M, ignore for now.

- for **lead-like** (300 < MW <= 450) and **fragment-like** (MW <= 300) libraries, keep each subset file at 100,000 mol max.
- for **drug-like** (MW > 450) library, keep each subset file at  80,000 mol max.

#################################

**Prepare ligand library for Glide docking**

- from Schoridnger/2016-03+, the _'-r 1'_ flag is decrepated. Before, it controls the ring conformation "add input ring conformation if available". 18.08.22 added property filter to remove useless compds (known **reactives** and **PAINS**). while removing known reactives (acyl chloride, aldehyde, etc) is useful for fragment-like library, there are very few reactives in lead-like and drug-like libraries. 

- This setup will use [OPLS3e force field](https://doi.org/10.1021/acs.jctc.8b01026) to generate the 3D structure for the input ligand. Use the universal SDF.gz format. It generates one(1) most preferred charged state/tautomer for each ligand at pH=7.2 +/- 0.3, becuase most screening are against targets in the [cytosol (pH~7.2)](https://doi.org/10.1042/bj2500001) or blood/extraceullar domain (pH~7.4). For targets in [specialized organelles](https://doi.org/10.1152/physiologyonline.2002.17.1.1) (peroxisome, mitocondria, lysosome, etc) that have lower pH, a new set of ligands with charged state/tautomer appropiate for that pH should be generated. It also generates the stereoismers, tho it is recommended to use input with pre-defined chiral centers available.

- Athough removing PAIN compounds is a good idea, the SMARTS patterns used in the substructure recognition are greedily "dirty"[(1)](https://doi.org/10.1021/acs.jcim.6b00465) [(2)](https://doi.org/10.1021/acschembio.7b00903) [(3)](https://blogs.sciencemag.org/pipeline/archives/2018/08/06/pains-filters-in-the-real-world) and can flag down okay compounds as false positives. This has happened to the lead compound in [Dyrk1A type-II inhibitor project](https://doi.org/10.1016/j.ejmech.2018.08.007), which turns out it is a false positive. Though I now do not recommend doing a PAINS pre-filtering, it is better to understand the key features [(1)](https://doi.org/10.1021/jm901137j) [(2)](https://www.nature.com/news/chemistry-chemical-con-artists-foil-drug-discovery-1.15991) and pick out potential PAINS by eye in the hit-pick stage, then suggest watching these guys out during testing in orthogonal assays.

```
time ${SCHRODINGER}/ligprep -WAIT -LOCAL \
        -i 2   -epik -We,-ph,7.2,-pht,0.3,-ms,1 \
        -s 1   -t 1 \
        -bff   16 \
#        -f $homedir/all-1.sch.propt_filter.cflt \
        -ismi  $smi_dir/SDIN.smi \
#       -isd   /sc/hydra/projects/schlea02a/8_lib/zinc_lead13/glide/SDIN.sch.sdf.gz \
        -osd   SDIN.sch.sdf.gz \
        -HOST  localhost:$schcpu \
        -NJOBS $schnjob
echo $!

## Use Schrodinger Canvas to filter certain pains ligands
#  $SCHRODINGER/utilities/canvasSearch -WAIT -LOCAL \
#        -isd   SDIN.sdf.gz \
#        -osd   SDIN.sch.sdf.gz \
#        -osd2  SDIN.fail.sdf.gz -filter  \
#        -file  $homedir/all-1.sch.pains_filter.cflt \
#        -JOB   JBNAME \
#        -HOST  localhost:$schcpu
```

- To run this in batch mode on a HPC cluster, use the shell script **run_ligprep_lsf.csh** to modify the template LSF script and submit the jobs:
```
> A_docking_scripts
        |---- /1_prep_ligands
                    |---- run_ligprep_lsf.csh
                    |---- ligprep.tmpl.lsf
```

#################################

**Prepare ligand library for OpenEye Docking**

- OpenEye preparation requires multiple steps. First with _'fixpka'_ to generate one(1) most relevant tautomoer/charged state for the ligand as pH=7.4, then _'flipper'_ to generate multiple stereoisomers (tho it is recommended to use starting ligands with their chiral center predefined), then use _'omega2'_ to generate conformer library. For conformer generation, recommend to get _-maxconfs 300_, since for larger molecules (like sorafenib) the default 200 maximum conformer is not enough to cover the conformational space of ligands with more rotatable bonds.
- Although OMEGA uses the older MMFF94 force field with electrostatics modification by the original developers of force fields (Bayly, Darden are with [OE](https://www.eyesopen.com/staff)!), the generated small molecule conformations are actually similar to crystallographic poses and appropiate for use in docking. For atoms with unique hybridization state (e.g. "the almost aromatic" ring NC1=NN=CCS1 in lead compound in [Dyrk1A type-II inhibitor](https://doi.org/10.1016/j.ejmech.2018.08.007)), conformation will look bad in comparison to those generated by the all-new OPLS3e or QM calculation. But given the extreme efficiency and quality for vast majority of the ligands, results from OMEGA is still very much usable.

```
fixpka  -in input.smi    -out input.oe.smi
flipper -in input.oe.smi -out input.oe.ism
omega2 \
  -mpi_np   $cpu      \
  -in       input.oe.ism    \
  -out      input.oeb.gz    \
  -prefix   prefix    \
  -maxconfs 300       \
  -progress none      \
  -strict   false
```

- To run this in batch mode on a HPC cluster, use the shell script **run_omega_lsf.csh** to modify the template LSF script and submit the jobs:
```
> A_docking_scripts
        |---- /1_prep_ligands
                    |---- run_omega_lsf.csh
                    |---- omega.tmpl.lsf
```

#######################################################################################
# Docking
```
> ${OPENEYE}/bin/make_receptor

> ${OPENEYE}/bin/fred -mpi_np $cpu
```
- OpenEye docking has 2 steps, (1) generating the _rigid_ receptor grid file, and (2) _rigid_ docking of ligand conformers to receptor grid. For the most parts the default settings by OpenEye will work for 95% of the cases becuase they _really really_ optimized everything. As for the scoring function ChemGauss4, it is a combination of a knowledge-based Gaussian potential for shape complementary plus a modified ChemScore/PFP/GoldScore(?) empirical-based chemical feature scoring.

- **Grid generation** - this is done through the command-line GUI *make_receptor*. Default setting will work for almost all cases. May want to incresae the size of the docking pocket. Can setup directional H-bond (donor or acceptor) or SMARTS-defined spherical constraints if these interactions are known and critical.

- **Docking** - this is done through the command-line program *fred*. The default setup will work for 99% of the case. Write out the docking results to **sdf.gz** format, not the *oeb.gz* format. Be sure to change the **-hitlist_size** to **0** (save all), otherwise only _500_ (default) will be saved and the rest are lost.

- To run this in batch mode on a HPC cluster, use the shell script **fred_submit_lsf.csh** to modify the template LSF script and submit the jobs:
```
> A_docking_scripts
        |---- /2_virtual_screen
                    |---- fred_submit_lsf.csh
                    |---- fred_docking.repeat.lsf
```

#################################

```
> ${SCHRODINGER}/glide  glide-grid_generation.inp

> ${SCHRODIMGER}/glide  glide-dock.inp
```
- Schrodinger's Glide docking has 2 steps, (1) generating the _rigid_ receptor grid file, and (2) _flexible_ ligand docking to receptor grid. Most of the inital setup is done through Schrodinger's Maestro GUI, but certain settings can only be done via altering the input files written out by the GUI, and then run the modified input files in command-line mode. As for the scoring function GlideScore, it is an empirical-based scoring function, plus the options to consider H-bond-like _Aromatic H_ and _Halogen_ bonding [interactions](https://doi.org/10.1021/jm100112j).

- **Grid generation** - the default vdw scaling of 1.0x is alright, but I have found a softened vdw scaling (0.75x - 0.80x) to be optimal for most cases (kinases, transporers, etc.). But really, need to try a couple vdw scalings (0.60x - 1.20x, in 0.05 increment) on knonwn ligands to test what is the best vdw scaling factor for the particular receptor. If there is no knonwn ligand available, then a 0.80x vdw scaling is a safe bet. Cannot edit this in the GUI, have to edit it in the input text file written out by the GUI. **Always** turn on **Aromatic H**and **Halogen bonds** options. If you know certain interaction **must** happen, e.g. kinase inhibitors H-bonding the backbone amide of hinge residue, type-II kinase inhibitors occupying the hydrophobic DFG-pocket, then you can setup contraint location (spherical, or directional for H-bond) and preferred ligand types in SMARTS patterns.

- **Docking** - the default _Standard Precision (SP)_ setting works for most cases. **Do not** use _High-throughput virtual screen (HTVS)_ setting to "save time" as the manual recommended. Tried it and this setting failed to capture many known ligands even in the top 50% ranking. **Extreme Precision (XP)** is useless and slow for most cases as it uses very harsh and rigid vdw scaling, which fails in most cases too. Only ever use it if you are working on a congeneric series of a scaffold that is developed based on a crystal structure. **Always** turn on **Aromatic H** and **Halogen bonds** options. Can use constraints previously defined in Grid generation by requiring _X_ of _n_ constraints must be fulfilled. Also, save the docking results with **no** receptor pose, and in **SDF** format, don't use **MAE** format.

- To run this in batch mode on a HPC cluster, use the shell script **glide_submit_lsf.csh** to modify the template LSF script and submit the jobs. Schrödinger does not write out a readable score file (.rept) like it used to after ~ ver.2016-09, and the pose files are in .gz format. The data will need to be cleaned up to get the more compact .bz2 format and to extract the scores from the pose files. Run the **glide_clean_folder.csh** with a list of the docking folder name (concatinate the _*.clean_ files generated during the initial submission step) to get the **/1_data** folder with all the cleaned data.
```
> A_docking_scripts
        |---- /2_virtual_screen
                    |---- glide_submit_lsf.csh           # control glide submission
                    |---- glide-dock.template.lsf        # docking LSF template
                    |---- glide-dock.HTVS_SP.template.in # docking input template
                    |---- glide-grid_template.in         # grid input template
                    |
                    |---- glide_clean_folder.csh   # summarize batch docking folders
                    |---- glide_clean_score.py     # make readable Glide score file
```

#################################

```
> ${SCHRODINGER}/ifd
```
- **Induced-Fit Docking (IFD)** from Schrodinger performs _flexible_ ligand docking to _flexible_ receptor. It combines softened docking with receptor minimization, and optional loop modeling. Most of the setup is done through _Maestro IFD GUI_. 

- The initial rigid docking uses a much softened _vdw radii_ (recommend 0.5-0.6x vdw scaling), followed by receptor "side chain" minimization, followed by regular rigid Glide docking (recommend 0.8x vdw scaling) and rescoring for the final result. User can define which residue's side chain can be minimized (default: those within 5A distance from initial ligand docking pose), or user-defined residues. User can also designate which residue side chain is "removed" in the initial docking to avoid steric hindrance or overwhelming bias in the intial docking (see the supplementary of article [_(Lazarus Ung JACS 2019 ULK4)_](https://doi.org/10.1021/jacs.9b10458), where IFD for prodrug R788 required such maneuver). 
- **Also, a hidden option to use "Flexible Loops" is accessible by modifying the input text file written out by IFD**. This is described in the manual but is not available in the Meastro IFD GUI. It is useful for cases where binding site involves a flexible loop, e.g. P-loop of protein kinase.

#######################################################################################
# Collecting Results of Standard Rigid-Receptor Docking

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
- Collect the top scoring results from _one_ or _multiple_ docking files and generate a histogram to show the distribution of the docking scores (.png). The score files (.txt) and the docking pose files (.sdf) can be in zipped format but must have the same prefix, i.e.:
> test.sch_docked.txt.bz2
> test.sch_docked.sdf.bz2

#################################

```
> B_collect_scripts/7_general_docking_cluster.py
             [ Docking input(s): sdf ]   - accept multiple files "*.sdf,x.sdf"
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
- Cluster the supplied molecules using a _fingerprint-based chemical similarity_ method with a user-defined Tanimoto cutoff value. A value of **0.40** works for most cases of lead-/drug-like molecules (MW 300-500+), while fragment-like may need a lower value, e.g. 0.37. The clustered results are placed into a PyMOL session file and the docking receptor is used. The clustered molecules are also output into a PDF file.

#################################

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

#################################

# Collecting Results of Schrodinger Induced-Fit Docking (IFD)
```
> B_collect_scripts/3_glide_IFD_mae2pdb.py
       [ ligand file: sdf ]         # full path to file
       [ directory of IFD result ]  # full path to directory
       [ output filename prefix ]

   e.g.> *.py   atp.sdf   InducedFit_1   test_IFD

result > test_IFD.pse
      ***   start outside of the parent IFD folder, 
      ***   result files will be in the IFD folder

```
- Consolidate Schrodinger's IFD top results scattered within the IFD working directory and put them into a PyMOL session file. The IFD scores are within a result CSV file, which lists the IFD docking pose file associating with the ligand used. This script extracts those information to file the corresponding .MAE files, converts them back into .PDB files, and collects them into the PyMOL session file.


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
  ifd
  prepwizard
  structconvert
```
```
csh/tcsh          # shell
python            # 3.6.8+
  numpy           # 1.16.2+
  pandas          # 0.24.2+
  matplotlib      # 3.0.3+
  seaborn         # 0.9.0+
  rdkit           # 2019.09.02
  tabulate        # 0.8.3+
  pathos          # 0.2.3+
  tqdm            # 4.31.1+
  argparse        # 1.1
  html            # 
  gzip            # 
  bz2             # 
```
