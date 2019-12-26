
## Use RDKit ETKDG conformation generation method
## problem is it is very slow compare to OMEGA or even ligprep
## and sometimes amide bond is in linear form, need UFF minimization
## to fix this issue, but that will slow it down even further.
## Around 1.5s per molecules in the lead-like range (5-8 rb)
## Definitely not optimal for large volume generation

../../A_docking_scripts/3_conformer_gen/3_rdkit_conformer_gen.py    \
  -in antibiotics_small.sdf.bz2   \
  -out test_small                 \
  -rmsd .62                       \
  -addh                           \
  -ff

## test.conf.log
#  number of conformer generated for each input molecule
#
## test.sdf.bz2    
#  sdf file with all conformers of all input molecules
#  took ~4 mins to generate 19 molecules
#
## test_small.sdf.bz2
#  smaller set, 5 molecules
#
#  19.12.24

