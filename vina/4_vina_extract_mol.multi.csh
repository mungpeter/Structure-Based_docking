#!/bin/csh

# extract and convert multiple sets of vina pdbqt to sdf for consensus 

if ($#argv != 1) then
  echo 
  echo "  ## Usage: x.csh [list of vina selection files] "
  exit
endif

foreach list (`cat vina.list`)
  set name = `basename $list .txt`
  ~/9_scripts/1_Docking/vina/4_vina_extract_mol.py $list $name.docked

end
