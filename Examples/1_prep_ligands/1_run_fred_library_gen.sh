

## steps to generate OpenEye ligand conformer library for FRED docking

fixpka                            \
  -in  antibiotics.smi            \
  -out antibiotics.oe.smi


flipper                           \
  -in  antibiotics.oe.smi         \
  -out antibiotics.oe.ism           


omega2                            \
  -mpi_np   4                     \
  -in       antibiotics.oe.ism    \
  -out      antibiotics.oeb.gz    \
  -prefix   antibiotics           \
  -maxconfs 300                   \
  -progress percent               \
  -strict   false

## antibiotics.smi 
#  original input, noted all chiral centers are fixed
#
## antibiotics.oe.smi
#  molecules are checked for pKa and ONE tautomer most likely to be true at 
#  pH=7.4 is saved
#
## antibiotics.oe.ism
#  stereoisomers of each molecules are generated
#
## antibiotics.oeb.gz
#  3d conformer library 
#
#  19.12.24
