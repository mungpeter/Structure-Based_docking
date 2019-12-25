

## search the supplied molecules for analogs that match
## the substructure(s)


../../1_search_analogs_substructure.py            \
  recpt.zfg19.sch_top500.pick.sdf.bz2             \
  "c1ccccc1N,C(F)(F)F"                            \
  recpt.zfg19.substr_phen,recpt.zfg19.substr_cf3  \
  sdf


../../0_search_analogs_substructure.py            \
  recpt.zfg19.sch_top500.pick.smi.bz2             \
  "c1ccccc1N,C(F)(F)F"                            \
  recpt.zfg19.substr_phen,recpt.zfg19.substr_cf3  \
  smi


## recpt.zfg19.substr_phen.sdf/.smi
#  supplied molecules that matched substructure 'c1ccccc1N'
#
## recpt.zfg19.substr_cf3.sdf/.smi
#  supplied molecules that matched substructure 'C(F)(F)F'
#
#  19.12.20
