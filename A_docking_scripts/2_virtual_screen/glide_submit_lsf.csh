#!/bin/csh

##########################################################################
#
#	Peter M.U. Ung @ MSSM
#	v1.0	- 15.06.16	
#
#	Run Schrodinger's Glide. Optional to choose to prepare protein(s)
#	with default settings.
#
#	Require all setup files in the home directory
#	Dock precision level: HTVS or SP
#	Grid file has the extension: .grid.zip
#
##########################################################################

set schro_ver  = 'schrodinger/2019-1'
set grid_templ = glide-grid_template.in
set dock_templ = ''
set lsf_templ  = ''
#set dock_templ = glide-dock.HTVS_SP.template.in
#set lsf_templ  = glide-dock.template.lsf
#set grid_templ = glide-grid_template.mkk7-cido.in
#set dock_templ = glide-dock.HTVS_SP.template.mkk7-cido.in
#set dock_templ = glide-dock.ulk4_d-in.template.in 

if ($#argv != 7) then
  echo ''
  echo '    Usage: x.csh'
  echo '             [List of Protein] [List of Ligand Database]'
  echo '             [PDB Directory]'
  echo '             [glide-dock template LSF] [glide-dock template input]'
  echo '             [Dock Precision (HTVS|SP)]'
  echo '             [Prepare Protein? (1|0)]'
  echo '             [Run Docking?     (1|0)]'
  echo ''
  echo "    - $schro_ver"
  echo "    - $lsf_templ"
  echo "    - $grid_templ"
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
if (! -e $grid_templ) then
  echo ''
  echo "    Error: Missing $grid_templ"
  echo ''
  exit
endif
#if (! -e $lsf_templ) then
#  echo ''
#  echo "   Error: Missing $lsf_templ"
#  echo ''
#  exit
#endif

set pdb_list  = $argv[1]
set lig_list  = $argv[2]
set pdb_dir   = $argv[3]
set lsf_templ = $argv[4]
set inp_templ = $argv[5]
set precision = $argv[6]
set prepare   = $argv[7]
set docking   = $argv[8]

##########################################################################

foreach pdb_file (`cat $pdb_list`)
  set pdb_name   = `basename $pdb_file .pdb`

  ## Prepare Protein/Grid for Glide
  if ($prepare == 1) then
    echo ''
    echo "##  Using $schro_ver ##"
    echo ''
    module load $schro_ver

    cd $pdb_dir
#    if (! -e 'schrodinger.hosts') then
#cat > schrodinger.hosts << EOF
# Set "tmpdir" to a temporary directory in your scratch folder.
#     e.g. /sc/hydra/scratch/your_user_ID/.tmp_schrodinger
#name:        localhost
# tmpdir:      /hpc/users/ungp01/.tmp_schrodinger
# Name: mothra_serial_interactive
# Host: interactive1
# Processors: 2500
# Queue: LSF
# Qargs:  -m mothra -W 30
# schrodinger: /hpc/packages/minerva-common/schrodinger/2016-1
#EOF
#    endif

    ## skip protein preparation if already prepared
    if (! -e $pdb_name.maegz) then
      ${SCHRODINGER}/utilities/prepwizard -WAIT -SAVE -NOJOBID -NOLOCAL \
        -disulfides -rehtreat -captermini -fillsidechains \
        -propka_pH '7.0' -fix -f '3' \
        -j $pdb_name.run \
        $pdb_file $pdb_name.maegz
    endif

    sed "s/GNAME/$pdb_name/g" ../$grid_templ | \
    sed "s/GPROTNAME/$pdb_name.maegz/g" \
      > glide-grid.$pdb_name.in

    ${SCHRODINGER}/glide -WAIT -SAVE -NOJOBID -NOLOCAL -OVERWRITE \
      glide-grid.$pdb_name.in
    cd ..
  endif

#######################################################
  if ($docking == 1) then
    foreach lig_db (`cat $lig_list`)
  
      if (`echo $lig_db | grep '.sdf.gz'` != '') then
        set lig_name = `basename $lig_db .sdf.gz`
      else
        set lig_name = `basename $lig_db .maegz`
      endif

      sed "s/GGRIDNAME/$pdb_name.grid.zip/g" $dock_templ | \
      sed "s/GDOCKPRECIS/$precision/g" | \
      sed "s/GLIGNAME/$lig_db/g" \
        > $pdb_name.$lig_name.in

      sed "s/JBNAME/$pdb_name.$lig_name/g" $lsf_templ | \
      sed "s+PDBDIRECTORY+$pdb_dir+g"  | \
      sed "s/INP/$pdb_name.$lig_name.in/g" \
        > $pdb_name.$lig_name.lsf
      echo $pdb_name.$lig_name > $pdb_name.$lig_name.clean

      qsub $pdb_name.$lig_name.lsf
    end
  endif

end
