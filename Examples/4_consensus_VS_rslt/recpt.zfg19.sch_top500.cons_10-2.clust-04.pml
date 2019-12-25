load ../recpt.pdb.bz2, recpt
show cartoon, poly
hide lines
color white, poly
color cyan, org
show sticks, org and not resn NMA+ACE
set_bond stick_radius, .15, recpt and org
create ref_lig, recpt and org and not resn NMA+ACE
show lines, byres poly within 5 of ref_lig
hide sticks, recpt and org
load recpt.zfg19.sch_top500.cons_10-2.sdf, recpt.zfg19.sch_top500.cons_10-2
load _TEMP.clust.0.sdf, clust.1
dist HB.1, poly, clust.1, mode=2
load _TEMP.clust.1.sdf, clust.2
dist HB.2, poly, clust.2, mode=2
load _TEMP.clust.2.sdf, clust.3
dist HB.3, poly, clust.3, mode=2
show sticks, org
set valence
hide (h. and (e. c extend 1))
util.cbas
center org
zoom org
set dash_gap, 0.25
hide labels
show mesh, byres poly within 4 of ref_lig
set mesh_width, 0.1
set light_count, 1
set ray_opaque_background, off
color white, poly
color cyan, recpt and org
color cyan, ref_lig
util.cnc
set ray_trace_mode, 1
set ray_trace_gain, 0.008
set ray_trace_color, black
set pse_export_version, 1.70
save recpt.zfg19.sch_top500.cons_10-2.clust-04.pse
quit
