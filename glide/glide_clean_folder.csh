#!/bin/csh

#	Peter M.U. Ung @ MSSM
#	v1.0 	- 15.06.20
#
#	Rename the Glide result filenames and clean up the folders
# 	Create 1_data folder that contains all the Glide docking results

set proj = /home/pmung/Dropbox/9_scripts/1_Docking/glide
if (! -e $proj/glide_clean_score.py) then
  echo ''
  echo '  Error: Cannot find glide_clean_score.py'
  echo ''
  exit
endif

if ($#argv != 1) then
  echo ''
  echo '   Usage: x.csh [list of directory name (Glide result filename)]'
  echo ''
  exit
endif

echo `pwd`
set dir = $argv[1]

mkdir 1_data

foreach folder (`cat $dir`)
  cd $folder
  set name = `basename $folder .sch`

  mv $folder.rept $name.sch_docked.raw
  mv $folder\_lib.sdfgz $name.sch_docked.sdf.gz
  grep 'FOR LIGAND' $folder\_subjobs.log > $name.sch_log.txt
  grep 'SKIP LIG' $folder\_subjobs.log >> $name.sch_log.txt
  grep 'Total elapsed time' $folder\_subjobs.log >> $name.sch_log.txt

  $proj/glide_clean_score.py \
    $name.sch_docked.raw $name.sch_docked.txt

  gunzip $name.sch_docked.sdf.gz
  bzip2 $name.sch_* &

  rm $folder-*_in.sdf.gz *_subjobs.* temp.schrodinger.hosts

  cp $name.sch_* ../1_data
  cd ..

end

cd 1_data
mkdir 1_log
mv *log.txt.bz2 *raw* 1_log
