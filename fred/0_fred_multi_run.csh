#!/bin/csh

if ($#argv != 3) then
  echo ""
  echo "  Usage: > x.csh"
  echo "             [list of Protein oeb.gz] *mod.oeb.gz"
  echo "             [ligand library]"
  echo "             [number of CPU]"
  echo ""
  exit
endif

set pdb  = $argv[1]
set lig  = $argv[2]
set cpu  = $argv[3]

foreach oeb (`cat $pdb`)
  set name = `basename $oeb .mod.oeb.gz`
  echo " ## Running on $name"

#  timeout  3h  \
  fred -mpi_np $cpu \
         -receptor $oeb \
         -dbase    $lig \
         -prefix   $name \
         -hitlist_size 0 \
#         -save_component_scores true \
         -docked_molecule_file $name.fred_docked.sdf \
         -score_file $name.fred_docked.txt

  bzip2 $name.fred_docked.sdf $name.fred_docked.txt &
  echo "  ## Done with $name"
end


#######################################33

# v1	14.03.27
# v2	19.02.18 fred does not use oempirun anymore
