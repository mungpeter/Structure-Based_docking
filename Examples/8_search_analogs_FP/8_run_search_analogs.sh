#!/bin/csh

## search the supplied molecules for analogs that match with
## Fingerprint similarity >= Tanimoto cutoff 

../../1_FP_search_analogs.py          \
  query_mol.smi                       \
  recpt.zfg19.sch_top500.pick.sdf.bz2 \
  recpt.zfg19.query_search.fp_all_03  \
  0.3                                 \
  all


## recpt.zfg19.query_search.fp_all_03.smi
#  unique analogs selected from the library filfulling "FP Tc total >= 0.3"

## recpt.zfg19.query_search.fp_all_03.pdf
#  summary of molecules in .sdf.bz2 found to match query_mol.smi, using
#  all FPs at 0.3 Tanimoto cutoff

## 19.12.19
