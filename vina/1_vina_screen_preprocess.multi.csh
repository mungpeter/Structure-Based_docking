#!/bin/csh

# for consensus of multiple sets of pdbqt results
# combine the temporary vina_docked score files and output them
# to the parent directory

# fold.list = a list of protein directory
# lig.list  = a list of ligand directory, under the protein directory

if ($#argv != 1) then
  echo
  echo "  ## Usage: x.csh [prefix of ligand dataset]"
  echo "              * need fold.list and lig.list "
  echo
  exit
endif

set prefix = $argv[1]	# Prefix of type of ligand dataset

foreach folder (`cat fold.list`)
  cd $folder
  echo $folder
  @ i = 0
  foreach lig (`cat ../lig.list`)
    cd $lig
    echo $lig
    cat *.temp | uniq | sort > $lig.vina_docked.txt
    cp  $lig.vina_docked.txt $lig.vina_score.txt
    cp $lig.vina_score.txt ../$folder.$lig.vina_score.txt
    cd ..
    cat $folder.*.vina_score.txt > x.x
    mv x.x ../$folder.$argv[1].vina_score.txt
  end
  cd ..
end
