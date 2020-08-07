#!/bin/csh

set HOME = /hpc/users/ungp01/1_kinase/3_vina
set FOLD = /scratch/ungp01/zinc-lead-now/autodock
set RSLT = /scratch/ungp01/vs-ksr-mek2
set PBS  = vina_vs.run.pbs
set FAIL = vina_vs.redock.pbs
set LIB  = 21_p0.$argv[1]
set STEP = 1000

echo "  ## checking number of ligand ##"
set CONF  = config.2.ksr-mek.z-lead.txt
set TOTAL = 70000	# `ll $FOLD/$LIB | wc -l`
set sub   = 1		# submit PBS: ON = 1; OFF = 0

set QUEUE = "small_24hr"
set WALL  = "walltime=24:00:00"
set NODES = "nodes=1:ppn=4"

echo "  ## Starting pbs_gen.vina.pl ##"
./pbs_gen.vina.pl \
  $HOME $FOLD $RSLT $LIB    \
  $TOTAL $STEP $PBS $CONF   \
  $QUEUE $WALL $NODES $sub

## Just generate the redock file
./pbs_gen.vina.pl \
  $HOME $FOLD $RSLT $LIB    \
  $TOTAL $STEP $FAIL $CONF  \
  $QUEUE $WALL $NODES 0

exit
