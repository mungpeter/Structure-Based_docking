#!/bin/csh

#       Peter M.U. Ung @ MSSM
#       v1.0    - 15.06.20
#
#       Rename the Glide result filenames and clean up the folders
#       Create 1_data folder that contains all the Glide docking results

set source = /home/pmung/Dropbox/9_scripts/1_Docking/glide

if (! -e $source/glide_clean_score.py) then
  echo ''
  echo '  Error: Cannot find glide_clean_score.py'
  echo ''
  exit
endif

set dir = $argv[1]	# List of rept file

mkdir 1_data

foreach rept (`cat $dir`)
#  cd $folder
  set name = `basename $rept .rept`

  mv $name.rept $name.sch_docked.raw
  mv $name\_lib.sdfgz $name.sch_docked.sdf.gz
  grep 'FOR LIGAND' $name\_subjobs.log >  $name.sch_log.txt
  grep 'SKIP LIG'   $name\_subjobs.log >> $name.sch_log.txt
  grep 'Total elapsed time' $name\_subjobs.log >> $name.sch_log.txt

  $source/glide_clean_score.py $name.sch_docked.raw $name.sch_docked.txt

  gunzip $name.sch_docked.sdf.gz
  bzip2 $name.sch_*

  rm $name-*_in.sdf.gz *_subjobs.* temp.schrodinger.hosts

  cp $name.sch_* ./1_data
#  cd ..

end

cd 1_data
mkdir 1_log
mv *log.txt.bz2 1_log

