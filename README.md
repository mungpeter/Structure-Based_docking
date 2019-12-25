# Structure-Based_docking
Scripts to handle structure-based molecular docking, single-receptor virtual screening, consensus results of multiple receptor docking, and etc.

There are 2 major parts:
```
  /A_docking_scripts -- Setting up and running of docking
  /B_collect_scripts -- Collecting the docking results
```
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

############################################################################################
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
