#!/bin/csh

## Steps:
## - collect all subsets of data into a top set
## - general cluster of the top set with ECFP_4 tanimoto cutoff @ 0.4

../../B_collect_scripts/5_general_screen_get_top.py    \
  -score   "*.sch_docked.txt.bz2"    \
  -sdf     "*.sch_docked.sdf.bz2"    \
  -top     1000                      \
  -dock    sch                       \
  -outpref single_VS.summary         \
  -hmax    -12                       \
  -hmin    -2


../../B_collect_scripts/7_general_docking_cluster.py          \
  single_VS.summary.sch_top1k.sdf           \
  0.4                                       \
  single_VS.summary.sch_top1k.clust-04      \
  ../recpt.pdb.bz2                          \
  -ec

bzip2 *sdf *txt *pse

## single_VS.summary.sch_top1k.sdf
#  top 1000 molecules, ranked 3D structures, collected from all subsets
#
## single_VS.summary.sch_top1k.txt
#  top 1000 molecules, ranked scores, collected from all subsets
#
## single_VS.summary.sch_top1k.histo.png
#  score distribution of docked molecules from subsets, including top 1000
#  cutoff score, and median/stdev scores
#
## single_VS.summary.sch_top1k.clust-04.sdf
#  top 1000 molecules, arranged according to cluster orders
#
## single_VS.summary.sch_top1k.clust-04.pml
#  setup file to construct summary of docking result in clusters in PyMOL
#
## single_VS.summary.sch_top1k.clust-04.pse (not generated here)
#  PyMOL session file of docking summary, with clustering at ECFP_4 
#  taminmoto cutoff at 0.4
#
## single_VS.summary.sch_top1k.clust-04.pdf
#  summary of docking result in clusters, in PDF format
#
## recpt.pdb.z2
#  a receptor PDB used in the VS for structural reference
#

## 19.12.19
