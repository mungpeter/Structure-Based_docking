#!/bin/csh

if ($#argv < 2 || $#argv> 3) then
  echo "";echo "    Usage: ${0}"
  echo "        [Protein PDB] [Box info] [Optional: Constraint file]";echo ''
  exit 1
endif
if ($#argv == 3) then
  set constr = $argv[3]
endif
set lig = $argv[2]

foreach pdb (`ls $argv[1]`)
  set out = `basename $pdb .pdb`
  echo "  ## Working on $out with $lig"
#  pdb2receptor \
#    -pdb $pdb \
#    -ligand_residue $lig \
#    -receptor $out.oeb.gz

  receptor_setup \
    -protein $pdb \
    -box $lig \
    -receptor $out.oeb.gz    

  if ($#argv == 3) then
    receptor_toolbox \
      -receptor $out.oeb.gz \
      -set_custom_constraints $constr 
  endif

end
