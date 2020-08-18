#!/bin/csh

##############################################################################
##
## 	Peter M.U. Ung @ MSSM
##
## 	v1.0 -- 13.10.11
##
## 	Run AD Vina as a batch; blasting thru the whole directory for ligand
##
##############################################################################

set fold = /home/pmung/8_libraries/zinc_frag_now_130214/autodock
set home = /home/pmung/1_kinase/2_ksr/1_dimer/vina
set lib  = 22_p0.0


if (! -e $lib) then
  mkdir $lib
endif

@ y = `ll $fold/$lib/*.pdbqt | wc -l`
@ i = 0
echo "  ## Total Ligand = $y ##"
foreach lig ($fold/$lib/*.pdbqt)
  set name = `basename $lig .pdbqt`
  echo Processing Ligand $name

 time \
  vina --config config.txt        \
       --ligand $lig              \
       --out    $lib/$name.pdbqt  \
       --log    $lib/$name.log

  set x = `grep -m 1 "VINA" $lib/$name.pdbqt`
  echo "$lib/$name.pdbqt::$x" >> $lib/$lib.vina_score.txt

  @ i++
  @ y--
  echo "  ## Finished $i | Remaining $y ##"
end
