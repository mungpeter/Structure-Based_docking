#!/bin/csh

##
## Peter MU Ung @ MSSM/Yale
##
## v1  20.01.31
##
## Split a large smiles files into smaller pieces, each with
## a maximum number of smiles
##

if (@$argv != 3) then
  echo ""
  echo "  x.csh"
  echo "    [ input smiles file ]"
  echo "    [ output prefix ]"
  echo "    [ maximum smiles in each piece ]"
  echo ""
  echo "e.g. x.csh lead.smi enm19.lead 100000"
  echo ""
  exit
endif

set inp_file = $argv[1]
set outpref  = $argv[2]
set lib_size = $argv[3]

set name = `basename $inp_file .smi`

grep -v smiles $inp_file | sort -uR > $name.temp.smi
echo "> number of smiles lines: " `cat $name.temp.smi | wc -l`

set max = `cat $name.temp.smi | wc -l`

@ i = $max
@ x = 1
while ($i > 0)

  echo " x = $x   --  i = $i"
  tail -$i $name.temp.smi | head -$lib_size > $outpref.$x.smi
 
  @ i = $i - $lib_size
  @ x = $x + 1
end
