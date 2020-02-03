#!/bin/csh

#	Peter M.U. Ung @ MSSM
#	v1.0 	- 15.06.20
#
#	Rename the Glide result filenames and clean up the folders
# 	Create 1_data folder that contains all the Glide docking results

set proj = /sc/hydra/projects/schlea02a/1_pmung/9_scripts/2_schrodinger/2_glide
#set proj = ~/1_kinase/1_super 

if ($#argv != 2) then
  echo ''
  echo '  Usage: x.csh [list of Glide docking directory]'
  echo '               [output folder name]'
  echo ''
  echo "  Source folder: $proj"
  echo ''
  exit
endif
if (! -e $proj/glide_clean_score.py) then
  echo ''
  echo '  Error: Cannot find glide_clean_score.py'
  echo ''
  exit
endif

echo current folder: `pwd`
set dir = $argv[1]
set out = $argv[2]

if (! -e $out) then
  mkdir $out
endif

foreach folder (`cat $dir`)
  echo $folder
  if (! -e $folder) then
    echo Cannot find folder: $folder
  else
    echo --Enter folder: $folder
    cd $folder
    set name = `basename $folder .sch`
    echo --Work on: $name

    mv $folder.rept $name.sch_docked.raw
    mv $folder\_lib.sdfgz $name.sch_docked.sdf.gz
    grep 'FOR LIGAND' $folder\_subjobs.log > $name.sch_log.txt
    grep 'SKIP LIG' $folder\_subjobs.log >> $name.sch_log.txt
    grep 'Total elapsed time' $folder\_subjobs.log >> $name.sch_log.txt

    $proj/glide_clean_score.py \
      $name.sch_docked.raw $name.sch_docked.txt

    gunzip $name.sch_docked.sdf.gz

    rm $folder-*_in.sdf.gz *_subjobs.* 
#    rm temp.schrodinger.hosts

    cp $name.sch_* ../$out
    cp $name.*in   ../$out
    cd ..

  endif

end

cd $out
mkdir 1_log
mv *log.txt *raw *in 1_log
bzip2 *.sdf &
bzip2 *.txt &
cd 1_log
bzip2 *raw *log.txt &


