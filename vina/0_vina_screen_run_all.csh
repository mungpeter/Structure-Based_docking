#!/bin/csh

##########################################
##
##      Peter M.U. Ung @MSSM
##
##      v1.0    - 13.11.03
##	v2.0	- 13.11.21 - retired 4_vina_screen_pymol.pl
##
##      Purpose: To run through Vina screening result   
##
##########################################


## Directory to scirpt
set script   = ~/Dropbox/9_scripts/1_Docking/vina

## Home Directory with the Vina docking results set in folders
set home_dir = /home/pmung/xxx_data/1_kinase/3_vina/4_ksr-ksr/2_lead

## Preprocess Directory, where to store the extracted score *.vina_score.txt
set proc_dir = $home_dir/result

## Reference Directory with reference list of pdbqt->ZINC
set ref_dir  = ~/8_libraries/1_zinc/zinc_lead_now_130204/zinc-ref

## Score Directory with vina_score.txt, the docking/scoring result list
set top_dir  = $home_dir/result

###################################################################
##
set pre_proc  = 0	# Preprocess to extract the Vina scores
set get_topV  = 1	# Rank and output the top-scoring poses
set extract   = 1	# Extract the top-ranking molecule
##set gen_pymol = 1       # Generate the pymol session with top-ranking mol
set cluster   = 1	# Cluster the top-ranking mol and print out PDF

###################################################################
## Individual ligand result folder for pre-processing
set lig_fold  = "21_p0.*"

## Rank and output the top-scoring poses
#set vina_rslt = .vina_score.txt
set sub_num   = 100
set all_num   = 2500
set prefix    = ksr-ksr

## Extract the top-ranking molecules
set top_list  = $prefix.vina_top$all_num.txt
set format    = sdf	# For clustering process, must use SDF
set extr_file = $prefix.vina_top$all_num.$format

## Receptor PDB for PyMOL
set receptor  = $home_dir/dimer.ftmap.2.pdb

## Clustering of the top-ranking molecules
set in_sdf    = $extr_file		# SDF for clustering
set cutoff    = 0.4			# Tanimoto coefficient cutoff
set clust_sdf = $prefix.vina_top$all_num.clust-$cutoff	# SDF with clusters
###################################################################
cd $home_dir
if (! -e $proc_dir) then
  mkdir $proc_dir
endif
cd $proc_dir
rm TEMP.*.png

if ($pre_proc == 1) then
  cd $home_dir
  echo "1_vina_screen_preprocess.py"
  $script/1_vina_screen_preprocess.py $lig_fold
  cp "*.vina_score.txt" $proc_dir  
endif

if ($get_topV == 1) then
  cd $top_dir
  echo "2_vina_screen_get_top.py"
  $script/2_vina_screen_get_top.py f="*_score.txt*" $sub_num $all_num $prefix
  ## output lists of top-rank molecule and the histogram of score distribution
  ## [set name].vina_top$sub_num.txt	[set name].vina_top$sub_num.png
  ##     $prfix.vina_top$all_num.txt	       all.vina_top$all_num.png
endif

if ($extract == 1) then
  cd $home_dir
  echo "3_vina_screen_extract_mol.pl"
  echo $top_list
  $script/3_vina_screen_extract_mol.pl $ref_dir $top_dir/$top_list $format $extr_file $receptor
  mv $extr_file $extr_file.pse $top_dir
endif

##if ($gen_pymol == 1) then
##  cd $home_dir
##  echo "4_vina_screen_pymol.pl"
##  $script/4_vina_screen_pymol.pl $top_dir/$top_list $ref_dir $receptor
##endif

if ($cluster == 1) then
  cd $top_dir
  echo "5_vina_screen_cluster.py"
  $script/5_vina_screen_cluster.py $in_sdf $cutoff $all_num $clust_sdf $receptor
endif
