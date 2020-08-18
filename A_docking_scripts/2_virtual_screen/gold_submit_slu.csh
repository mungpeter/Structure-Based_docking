#!/bin/csh

##########################################################################
#
#       Peter M.U. Ung @ Genentech
#       v1.0    20.08.11
#
#       Run GOLD in MPI
#
#       Require all setup files in the home directory
#       Receptor .mol2 file is treated as Grid file
#
##########################################################################

set schro_ver  = 'csd/2019.0'
set dock_templ = ''
set slu_templ  = ''
#set dock_templ = gold_vs.template.conf
#set slu_templ  = gold-dock.template.slu

if ($#argv != 7) then
  echo ''
  echo '    Usage: x.csh'
  echo '             [List of Protein] [List of Ligand Database]'
  echo '             [Constraint File] [Protein Directory]'
  echo '             [GOLD-dock template SLU] [GOLD-dock template input]'
  echo '             [CPU per job]'
  echo ''
  echo "    - $schro_ver"
  echo "    - $slu_templ"
  echo "    - $dock_templ"
  echo ''
  exit
endif
#if (! -e $dock_templ) then
#  echo ''
#  echo "    Error: Missing $dock_templ"
#  echo ''
#  exit
#endif
#if (! -e $slu_templ) then
#  echo ''
#  echo "   Error: Missing $slu_templ"
#  echo ''
#  exit
#endif

set mol2_list  = $argv[1]
set lig_list   = $argv[2]
set constr     = $argv[3]
set mol2_dir   = $argv[4]
set slu_templ  = $argv[5]
set dock_templ = $argv[6]
set cpu_num    = $argv[7]

##########################################################################

foreach mol2_file (`cat $mol2_list`)

  set mol2_name = `basename $mol2_file .mol2`


  foreach lig_db (`cat $lig_list`)

    set lig_name = `basename $lig_db .sdf.gz`

    ## modify the tempalte slrum submission file
    sed "s/JBNAME/$mol2_name.$lig_name/" $slu_templ | \
    sed "s/XCAVITYX/$mol2_name.cavity.atoms/"   | \
    sed "s/PDBDIRECTORY/$mol2_dir/"  | \
    sed "s/XMOL2NAMEX/$mol2_file/"   | \
    sed "s/XCONSTRX/$constr/"  | \
    sed "s/XLIGNAMEX/$lig_db/" | \
    sed "s/XCPUX/$cpu_num/"    | \
    sed "s/INP/$dock_templ/"     \
      > $mol2_name.$lig_name.slu

    echo $mol2_name.$lig_name > $mol2_name.$lig_name.clean

    sbatch $mol2_name.$lig_name.slu
  end

end

