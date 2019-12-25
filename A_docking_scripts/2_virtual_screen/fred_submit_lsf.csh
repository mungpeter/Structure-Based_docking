#!/bin/csh

if ($#argv != 2) then
  echo ''
  echo '  > x.csh [prot folder name (oeb.list in folder)] '
  echo '          [ligand list]'
  echo ''
  echo ''
  exit
endif

set folder = $argv[1]	# Folder name 
set liglst = $argv[2]	# ligand list

@ i = 1
foreach prot (`cat $folder/oeb.list`)
 foreach lig (`cat $liglst`)
  set hit = 0

  sed "s/PROT/$prot/" fred_docking.repeat.lsf | \
  sed "s/LIG/$lig/" | \
  sed "s/FOLDER/$folder/" | \
  sed "s/NUM/$i/" | \
  sed "s/HIT/$hit/" \
      > fred_docking.repeat.lsf.$i

  bsub < fred_docking.repeat.lsf.$i

  @ i++
 end
end
