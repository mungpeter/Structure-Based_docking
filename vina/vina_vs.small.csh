#!/bin/csh

##############################################################################
##
##      Peter M.U. Ung @ MSSM
##
##      v1.0 -- 13.10.11
##
##      Run AD Vina as a batch; blasting thru the whole directory for ligand
##
##############################################################################

if ($#argv != 3) then
  echo "    ##  Usage: x.csh [lib_name] [start no.] [ending no.]  ##"
  exit
endif


set fold = /home/pmung/8_libraries/1_zinc/zinc_lead_now_130204/autodock/
set home = /home/pmung/xxx_data/1_kinase/3_vina
set lib  = $argv[1]	# 22_p0.0
set conf = config.1.allosteric.txt

set start = $argv[2]	# starting no. of ligand
set last  = $argv[3]	# ending no. of ligand


if (! -e $lib) then
  mkdir $lib
endif

set i = $start
@ y = $last - $i
while ($i <= $last)
  echo Processing Ligand $lib.$i

 time \
  vina --config $conf                     \
       --ligand $fold/$lib/$lib.$i.pdbqt  \
       --out    $lib/$lib.$i.pdbqt        \
       --log    $lib/$lib.$i.log

  set x = `grep -m 1 "VINA" $lib/$lib.$i.pdbqt`
  echo "$lib/$lib.$i.pdbqt::$x" >> $lib/$lib.$start.txt

  else
    echo "  #### Ligand $lib.$i.pdbqt is not found ####"
  endif

  @ i++
  @ y--
  echo "  ## Finished $i | Remaining $y ##"
end
