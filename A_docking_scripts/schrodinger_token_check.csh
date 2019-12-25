#!/bin/csh

set schcpu = 15
echo $schcpu
@ limit = (($schcpu * 5) + 1)
echo $limit
@ avail = (2000 - `$SCHRODINGER/utilities/licutil -used | grep 'SUITE_' | awk '{print $2}'`)
echo $avail
if ($avail < $limit ) then
  echo '\n'
  echo " ## GLIDE_MAIN Token is fewer than $limit, stopped the job ##"
  echo '\n'
  exit
endif

