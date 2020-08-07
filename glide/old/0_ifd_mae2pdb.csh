#!/bin/tcsh

#	Peter M.U. Ung @ MSSM
#	
#	v1.0	16.08.03
#
#	convert Schrodinger IFD-generated protein-lignad .maegz into .pdb
#	and append them into a multi-pdb file
#	** Obsolite -- cannot parse the result csv file. use the .py version
#

if ($#argv != 2) then
  echo ''
  echo '    Usage: x.csh [list of .maegz] [Prefix of multi-pdb file]'
  echo '    ** Obsolite -- Use the newer, .csv-parsing .py version '
  echo ''
  exit
endif


foreach mae (`cat $argv[1]`)

  set name = `basename $mae .maegz`

  $SCHRODINGER/utilities/structconvert \
    -imae $mae -opdb _temp.pdb
  grep -v 'MODEL' _temp.pdb | grep -v 'ENDMDL' > $name.pdb

  echo $name.pdb >> _temp.list 
end

/home/pmung/Dropbox/9_scripts/3_program/structures/0_build_multi_pdb.py \
  _temp.list $argv[2]

rm _temp.list
