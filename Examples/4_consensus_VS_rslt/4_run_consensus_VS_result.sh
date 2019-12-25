#!/bin/csh

## Steps:
## 1. for each receptor, collect all subsets of data into a top set
## 2. run through all receptors' result
## 3. generate consensus result
## 3. general cluster of the top set with ECFP_4 tanimoto cutoff @ 0.4

## Score file and Docking Pose file have the same prefix, differ only
## by their suffix, *.txt and *.sdf, respectively.
##
## indy_4f35.sch_top500.sdf.bz2 and indy_5ul9.sch_top500.sdf.bz2 are
## top-500 docking results from different xtal structures.
## score_file.list - a list of the name of the score files


# step 3, consensus result
../../9_consensus_best_pose.py        \
  score_file.list                     \
  2                                   \
  recpt.zfg19.sch_top500.cons_10-2    \
  10                                  \
  2                                   \
  -c

# step 4, cluster results
../../7_general_docking_cluster.py           \
  recpt.zfg19.sch_top500.cons_10-2.sdf       \
  0.4                                        \
  recpt.zfg19.sch_top500.cons_10-2.clust-04  \
  ../recpt.pdb.bz2                           \
  -ec

bzip2 *pse *sdf *txt

## indy_top500.cons_10-2.sdf
#  consensus docking result of 2 different receptors, of which those
#  ranked in top 10% in 2 of 2 receptors are selected. The best 
#  scoring pose of the selected ligand is used
#
## indy_top500.cons_10-2.txt
#  consensus scoring result of 2 different receptors, corresponds to
#  the consensus docking results, the best score of the selected
#  ligand docking pose is used
#
## indy_consensus_top500.clust-04.pml
#  pymol setup file for summary in MyMOL session
#
## indy_consensus_top500.clust-04.pse  (not generated here)
#  PyMOL session file of docking summary, with clustering at ECFP_4 
#  taminmoto cutoff at 0.4

## 19.12.19

