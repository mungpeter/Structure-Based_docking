#!/bin/csh

if ($#argv < 2 || $#argv> 3) then
  echo "";echo "    Usage: ${0}"
  echo "        [Protein PDB list] [Box info] [Optional: Constraint file]";echo ''
  exit 1
endif

set list = $argv[1]
set box  = $argv[2]
if ($#argv == 3) then
  set constr = $argv[3]
endif

foreach prot (`cat $list`)

  ~/Dropbox/9_scripts/1_Docking/fred/0_fred_receptor_setup.csh  \
    $prot $box $constr

end
