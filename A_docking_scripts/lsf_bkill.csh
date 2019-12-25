#!/bin/csh

if ( $#argv != 2 ) then
  echo ''
  echo '  > x.csh [start LSF id] [end LSF id]'
  echo ''
  exit
endif

set start = $argv[1]
set end  = $argv[2]
@ i = $start
@ j = 0
while ($i <= $end)
  bkill $i
  @ i++
  @ j++
end

echo 'Killed jobs: '$j
