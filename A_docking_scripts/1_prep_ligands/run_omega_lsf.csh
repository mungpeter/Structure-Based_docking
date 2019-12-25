#!/bin/csh

if ($#argv != 2) then
  echo ''
  echo '    Usage: x.csh  [smi or sdf]'
  echo '             [List of Ligand Sets]'
  echo ''
  echo '    e.g.:  x.csh smi lig.list'
  echo ''
  exit
endif
#if (! -e 'ligprep.tmpl.inp') then
#  echo ''
#  echo '    Error: missing ligprep.tmpl.in'
#  echo ''
#  exit
#endif
if (! -e 'omega.tmpl.lsf') then
  echo ''
  echo '    Error: missing omega.tmpl.lsf'
  echo ''
  exit
endif
##########################################################################
set type = $argv[1]
set list = $argv[2]

foreach file (`cat $list`)
  if ($type == 'smi') then
    set name = `basename $file .smi`
  else if ($type == 'sdf') then
    set name = `basename $file .sdf.gz`
  else
    echo '  # Failed - not the correct file type'
    exit
  endif

  echo $name
#  set inp  = $name.inp
  set lsf  = $name.lsf

  sed "s/SDIN/$name/g" omega.tmpl.lsf | \
  sed "s/SDFILE/$file/g"  | \
  sed "s/JBNAME/$name/g"    \
    > $lsf
  
  qsub $name.lsf

end
