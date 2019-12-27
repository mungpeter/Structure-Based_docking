

##

time ${SCHRODINGER}/ligprep -WAIT -LOCAL        \
        -i 2   -epik -We,-ph,7.2,-pht,0.3,-ms,1 \
        -s 1   -t 1                             \
        -bff   16                               \
        -ismi  antibiotics.smi        \
        -osd   antibiotics.sch.sdf.gz \
        -HOST  localhost:4            \
        -NJOBS 4
#        -f all-1.sch.propt_filter.cflt \


## Use filter to pick out molecules that fulfill the PAINS SMARTS patterns
## however, getting flagged as PAINS does not mean it must be PAINs cmpds,
## just that features of the molecules get picked up by the SMARTS patterns.
## There are cases compounds flagged as PAINS are actually safe and useful
## in hit-to-lead. So I stopped using this filtering and just inspect the
## final virtual screening of entire library visually

time  $SCHRODINGER/utilities/canvasSearch -WAIT -LOCAL \
        -isd   antibiotics.sch.sdf.gz         \
        -osd   antibiotics.sch_safe.sdf.gz    \
        -osd2  antibiotics.sch_pain.sdf.gz    \
        -filter  \
        -file  ../../A_docking_scripts/1_prep_ligands/all-1.sch.pains_filter.cflt \
        -JOB   pains_filter         \
        -HOST  localhost:4



## antibiotics.sch.sdf.gz
#  LigPrep prepared single conformer 3D structure of input molecule
#
#  Tautomer fixed at pH=7.2 +- 0.3 with only ONE tautomer output,
#  only 1 stereoisomer, using OPLS3e force field (v2018.0x and on)
#
## antibiotics.sch_safe.sdf.gz
#  compounds that do not get flagged by PAINS SMARTS patterns
#
## antibiotics.sch_pain.sdf.gz
#  compounds that get flagged by PAINS SMARTS patterns
#
#  19.12.24
