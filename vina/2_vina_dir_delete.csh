#!/bin/csh

if ($#argv != 1) then
  echo ''
  echo "  ## Usage: x.csh [directory prefix] ##"
  echo "      e.g.  x.csh 21_p0"
  echo "          : Delete all .pdbqt .log .txt .temp"
  echo ''
  exit
endif

set name = $argv[1]

set folder = "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16"
set num    = "1 2 3 4 5 6 7 8 9"

foreach fold ($folder)
  cd $name.$fold

  foreach i ($num)
    rm $name.$fold.$i*.pdbqt
    rm $name.$fold.$i*.log
    rm $name.*.temp
    rm $name.*.txt
  end
  cd ..
  rmdir $name.$fold
end
