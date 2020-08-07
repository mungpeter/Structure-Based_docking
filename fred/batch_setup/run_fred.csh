#!/bin/csh


set cpu  = 12

foreach oeb (`cat $argv[1]`)
  set name = `basename $oeb .1atp.oeb.gz`
  echo " ## Running on $name"

#  timeout  3h  \
  oempirun -np $cpu \
    fred -receptor $oeb \
         -dbase    /home/pmung/xxx_data/1_kinase/d_irak1/1_struc/1_drug/drug.oeb.gz \
         -prefix   $name \
#         -hitlist_size 0 \
#         -save_component_scores true \
         -docked_molecule_file $name.fred_docked.sdf \
         -score_file $name.fred_docked.txt

  echo "  ## Done with $name"
end
