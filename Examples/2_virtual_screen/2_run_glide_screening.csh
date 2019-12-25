#!/bin/csh

../../A_docking_scripts/2_virtual_screen/glide_submit_lsf.csh             \
  pdb.sch.list                   \
  zld19.sch.list                 \
  pdb_directory                  \
  ../../A_docking_scripts/2_virtual_screen/glide-dock.template.lsf        \
  ../../A_docking_scripts/2_virtual_screen/glide-dock.HTVS_SP.template.in \
  SP                             \
  1                              \
  1
