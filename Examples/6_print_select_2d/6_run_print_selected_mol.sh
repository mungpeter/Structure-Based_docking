#!/bin/csh

## When the desired molecules are known and their docking poses needed to 
## be extracted, generate a list of single column of molecules' ID in docking
## result
## For sdf format, zipped file can be accepted
## for smi format, only ASCII format is accepted

../../B_collect_scripts/6_general_print_mol.py              \
  recpt.zfg19.sch_top500.pick.sdf.bz2

## recpt.zfg19.sch_top500.pick.pdf
#  specific 'picked' docking results in PDF format

## 19.12.19
