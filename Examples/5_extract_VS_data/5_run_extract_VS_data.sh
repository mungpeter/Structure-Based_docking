#!/bin/csh

## When the desired molecules are known and their docking poses needed to 
## be extracted, generate a list of single column of molecules' ID in docking
## result

## pick out specific molecules from original file(s) based on mol ID
../../8_general_result_select.py        \
  recpt.zfg19.sch_top500.pick.txt.bz2   \
  recpt.zfg19.sch_top500.sdf.bz2        \
  Score


bzip2 *.sdf *.txt


## recpt.zfg19.sch_top500.pick.sdf
#  specific 'picked' docking poses collected from original file
#
## recpt.zfg19.sch_top500.pick.txt
#  specific 'picked' scoring result collected from original file
#
## recpt.zfg19.sch_top500.pick.pdf
#  specific 'picked' docking results in PDF format

## 19.12.19
