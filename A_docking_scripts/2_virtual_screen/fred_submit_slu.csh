#!/bin/tcsh

if ($#argv != 3) then
  echo ''
  echo '  > x.csh [prot folder name (oeb.list in folder)] '
  echo '          [ligand list]'
  echo '          [template slu file]'
  echo ''
  exit
endif

set folder = $argv[1]	# Folder name 
set liglst = $argv[2]	# ligand list
set slufl  = $argv[3]

@ i = 0
foreach prot (`cat $folder/oeb.list`)
 foreach lig (`cat $liglst`)
  set hit = 125

  sed "s/PROT/$prot/" $slufl | \
  sed "s/LIG/$lig/" | \
  sed "s/FOLDER/$folder/" | \
  sed "s/NUM/$i/" | \
  sed "s/HIT/$hit/" \
      > $slufl.$i

  sbatch $slufl.$i

  @ i++
 end
end
