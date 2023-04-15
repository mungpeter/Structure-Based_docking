#!/bin/tcsh

if ($#argv < 4 ) then
  echo ''
  echo '  > x.csh'
  echo '      [ prot folder name ]'
  echo '      [ receptor oeb list, in prot folder ]'
  echo '      [ ligand list ]'
  echo '      [ template slu file ]'
  echo '      [ fred|hybrid (def: fred) ]'
  echo '      [ saved hit per subset (def: 125) ]'
  echo ''
  exit
endif

set folder = $argv[1]	# Folder name
set recpt  = $argv[2]	# oeb list in Folder 
set liglst = $argv[3]	# ligand list
set slufl  = $argv[4]	# template slu file
set dock   = 'fred'	# docking method
set hit    = 125	# hit saved per subset
if ($#argv >= 5) then
  set dock = $argv[5]
  if ($#argv == 6) then
  set hit  = $argv[6]
  endif
endif

@ i = 0
foreach prot (`cat $folder/$recpt`)
  foreach lig (`cat $liglst`)

    cat $slufl              |\
    sed "s/PROT/$prot/"     |\
    sed "s/LIG/$lig/"       |\
    sed "s/DOCK/$dock/"     |\
    sed "s/FOLDER/$folder/" |\
    sed "s/NUM/$i/"         |\
    sed "s/HIT/$hit/"        \
      > $slufl.$i

  sbatch $slufl.$i

  @ i++
 end
end
