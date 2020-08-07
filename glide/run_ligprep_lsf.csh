#!/bin/csh

if ($#argv != 1) then
  echo ''
  echo '    Usage: x.csh'
  echo '             [List of Ligand Sets]'
  echo ''
  echo '    e.g.:  x.csh lig.list 4 "\/work\/dir\/with\/escape\/for\/slash"'
  echo ''
  exit
endif
if (! -e 'ligprep.tmpl.inp') then
  echo ''
  echo '    Error: missing ligprep.tmpl.in'
  echo ''
  exit
endif
if (! -e 'ligprep.tmpl.lsf') then
  echo ''
  echo '    Error: missing ligprep.tmpl.lsf'
  echo ''
  exit
endif
##########################################################################
set list = $argv[1]

foreach file (`cat $list`)
  set name = `basename $file .sdf.gz`
  echo $name
  set inp  = $name.inp
  set lsf  = $name.lsf

  sed "s/XFILENAME/$file/g" ligprep.tmpl.inp | \
  sed "s/XLIGNAME/$name/g" \
    > $inp

  sed "s/FILE/$file/g" ligprep.tmpl.lsf | \
  sed "s/LIGNAME/$name/g" | \
  sed "s/JBNAME/l_$name/g" | \
  sed "s/INP/$inp/g" \
    > $lsf
  
  qsub $name.lsf

end
