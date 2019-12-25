load 5ul9_prep.pdb.bz2, 5ul9_prep
show cartoon, poly
hide lines
color white, poly
color cyan, org
show sticks, org and not resn NMA+ACE
set_bond stick_radius, .15, 5ul9_prep and org
create ref_lig, 5ul9_prep and org and not resn NMA+ACE
show lines, byres poly within 5 of ref_lig
hide sticks, 5ul9_prep and org
load single_VS.summary.sch_top1k.sdf, single_VS.summary.sch_top1k
load _TEMP.clust.0.sdf, clust.1
dist HB.1, poly, clust.1, mode=2
load _TEMP.clust.1.sdf, clust.2
dist HB.2, poly, clust.2, mode=2
load _TEMP.clust.2.sdf, clust.3
dist HB.3, poly, clust.3, mode=2
load _TEMP.clust.3.sdf, clust.4
dist HB.4, poly, clust.4, mode=2
load _TEMP.clust.4.sdf, clust.5
dist HB.5, poly, clust.5, mode=2
load _TEMP.clust.5.sdf, clust.6
dist HB.6, poly, clust.6, mode=2
load _TEMP.clust.6.sdf, clust.7
dist HB.7, poly, clust.7, mode=2
load _TEMP.clust.7.sdf, clust.8
dist HB.8, poly, clust.8, mode=2
load _TEMP.clust.8.sdf, clust.9
dist HB.9, poly, clust.9, mode=2
load _TEMP.clust.9.sdf, clust.10
dist HB.10, poly, clust.10, mode=2
load _TEMP.clust.10.sdf, clust.11
dist HB.11, poly, clust.11, mode=2
load _TEMP.clust.11.sdf, clust.12
dist HB.12, poly, clust.12, mode=2
load _TEMP.clust.12.sdf, clust.13
dist HB.13, poly, clust.13, mode=2
load _TEMP.clust.13.sdf, clust.14
dist HB.14, poly, clust.14, mode=2
load _TEMP.clust.14.sdf, clust.15
dist HB.15, poly, clust.15, mode=2
load _TEMP.clust.15.sdf, clust.16
dist HB.16, poly, clust.16, mode=2
load _TEMP.clust.16.sdf, clust.17
dist HB.17, poly, clust.17, mode=2
load _TEMP.clust.17.sdf, clust.18
dist HB.18, poly, clust.18, mode=2
load _TEMP.clust.18.sdf, clust.19
dist HB.19, poly, clust.19, mode=2
load _TEMP.clust.19.sdf, clust.20
dist HB.20, poly, clust.20, mode=2
load _TEMP.clust.20.sdf, clust.21
dist HB.21, poly, clust.21, mode=2
load _TEMP.clust.21.sdf, clust.22
dist HB.22, poly, clust.22, mode=2
load _TEMP.clust.22.sdf, clust.23
dist HB.23, poly, clust.23, mode=2
load _TEMP.clust.23.sdf, clust.24
dist HB.24, poly, clust.24, mode=2
load _TEMP.clust.24.sdf, clust.25
dist HB.25, poly, clust.25, mode=2
load _TEMP.clust.25.sdf, clust.26
dist HB.26, poly, clust.26, mode=2
load _TEMP.clust.26.sdf, clust.27
dist HB.27, poly, clust.27, mode=2
load _TEMP.clust.27.sdf, clust.28
dist HB.28, poly, clust.28, mode=2
load _TEMP.clust.28.sdf, clust.29
dist HB.29, poly, clust.29, mode=2
load _TEMP.clust.29.sdf, clust.30
dist HB.30, poly, clust.30, mode=2
load _TEMP.clust.30.sdf, clust.31
dist HB.31, poly, clust.31, mode=2
load _TEMP.clust.31.sdf, clust.32
dist HB.32, poly, clust.32, mode=2
load _TEMP.clust.32.sdf, clust.33
dist HB.33, poly, clust.33, mode=2
load _TEMP.clust.33.sdf, clust.34
dist HB.34, poly, clust.34, mode=2
load _TEMP.clust.34.sdf, clust.35
dist HB.35, poly, clust.35, mode=2
load _TEMP.clust.35.sdf, clust.36
dist HB.36, poly, clust.36, mode=2
load _TEMP.clust.36.sdf, clust.37
dist HB.37, poly, clust.37, mode=2
load _TEMP.clust.37.sdf, clust.38
dist HB.38, poly, clust.38, mode=2
load _TEMP.clust.38.sdf, clust.39
dist HB.39, poly, clust.39, mode=2
load _TEMP.clust.39.sdf, clust.40
dist HB.40, poly, clust.40, mode=2
load _TEMP.clust.40.sdf, clust.41
dist HB.41, poly, clust.41, mode=2
load _TEMP.clust.41.sdf, clust.42
dist HB.42, poly, clust.42, mode=2
load _TEMP.clust.42.sdf, clust.43
dist HB.43, poly, clust.43, mode=2
load _TEMP.clust.43.sdf, clust.44
dist HB.44, poly, clust.44, mode=2
load _TEMP.clust.44.sdf, clust.45
dist HB.45, poly, clust.45, mode=2
load _TEMP.clust.45.sdf, clust.46
dist HB.46, poly, clust.46, mode=2
load _TEMP.clust.46.sdf, clust.47
dist HB.47, poly, clust.47, mode=2
load _TEMP.clust.47.sdf, clust.48
dist HB.48, poly, clust.48, mode=2
load _TEMP.clust.48.sdf, clust.49
dist HB.49, poly, clust.49, mode=2
load _TEMP.clust.49.sdf, clust.50
dist HB.50, poly, clust.50, mode=2
load _TEMP.clust.50.sdf, clust.51
dist HB.51, poly, clust.51, mode=2
load _TEMP.clust.51.sdf, clust.52
dist HB.52, poly, clust.52, mode=2
load _TEMP.clust.52.sdf, clust.53
dist HB.53, poly, clust.53, mode=2
load _TEMP.clust.53.sdf, clust.54
dist HB.54, poly, clust.54, mode=2
load _TEMP.clust.54.sdf, clust.55
dist HB.55, poly, clust.55, mode=2
load _TEMP.clust.55.sdf, clust.56
dist HB.56, poly, clust.56, mode=2
load _TEMP.clust.56.sdf, clust.57
dist HB.57, poly, clust.57, mode=2
load _TEMP.clust.57.sdf, clust.58
dist HB.58, poly, clust.58, mode=2
load _TEMP.clust.58.sdf, clust.59
dist HB.59, poly, clust.59, mode=2
load _TEMP.clust.59.sdf, clust.60
dist HB.60, poly, clust.60, mode=2
load _TEMP.clust.60.sdf, clust.61
dist HB.61, poly, clust.61, mode=2
load _TEMP.clust.61.sdf, clust.62
dist HB.62, poly, clust.62, mode=2
load _TEMP.clust.62.sdf, clust.63
dist HB.63, poly, clust.63, mode=2
load _TEMP.clust.63.sdf, clust.64
dist HB.64, poly, clust.64, mode=2
load _TEMP.clust.64.sdf, clust.65
dist HB.65, poly, clust.65, mode=2
load _TEMP.clust.65.sdf, clust.66
dist HB.66, poly, clust.66, mode=2
load _TEMP.clust.66.sdf, clust.67
dist HB.67, poly, clust.67, mode=2
load _TEMP.clust.67.sdf, clust.68
dist HB.68, poly, clust.68, mode=2
load _TEMP.clust.68.sdf, clust.69
dist HB.69, poly, clust.69, mode=2
load _TEMP.clust.69.sdf, clust.70
dist HB.70, poly, clust.70, mode=2
load _TEMP.clust.70.sdf, clust.71
dist HB.71, poly, clust.71, mode=2
load _TEMP.clust.71.sdf, clust.72
dist HB.72, poly, clust.72, mode=2
load _TEMP.clust.72.sdf, clust.73
dist HB.73, poly, clust.73, mode=2
load _TEMP.clust.73.sdf, clust.74
dist HB.74, poly, clust.74, mode=2
load _TEMP.clust.74.sdf, clust.75
dist HB.75, poly, clust.75, mode=2
load _TEMP.clust.75.sdf, clust.76
dist HB.76, poly, clust.76, mode=2
load _TEMP.clust.76.sdf, clust.77
dist HB.77, poly, clust.77, mode=2
load _TEMP.clust.77.sdf, clust.78
dist HB.78, poly, clust.78, mode=2
load _TEMP.clust.78.sdf, clust.79
dist HB.79, poly, clust.79, mode=2
load _TEMP.clust.79.sdf, clust.80
dist HB.80, poly, clust.80, mode=2
load _TEMP.clust.80.sdf, clust.81
dist HB.81, poly, clust.81, mode=2
load _TEMP.clust.81.sdf, clust.82
dist HB.82, poly, clust.82, mode=2
load _TEMP.clust.82.sdf, clust.83
dist HB.83, poly, clust.83, mode=2
load _TEMP.clust.83.sdf, clust.84
dist HB.84, poly, clust.84, mode=2
load _TEMP.clust.84.sdf, clust.85
dist HB.85, poly, clust.85, mode=2
load _TEMP.clust.85.sdf, clust.86
dist HB.86, poly, clust.86, mode=2
load _TEMP.clust.86.sdf, clust.87
dist HB.87, poly, clust.87, mode=2
load _TEMP.clust.87.sdf, clust.88
dist HB.88, poly, clust.88, mode=2
load _TEMP.clust.88.sdf, clust.89
dist HB.89, poly, clust.89, mode=2
load _TEMP.clust.89.sdf, clust.90
dist HB.90, poly, clust.90, mode=2
load _TEMP.clust.90.sdf, clust.91
dist HB.91, poly, clust.91, mode=2
load _TEMP.clust.91.sdf, clust.92
dist HB.92, poly, clust.92, mode=2
load _TEMP.clust.92.sdf, clust.93
dist HB.93, poly, clust.93, mode=2
load _TEMP.clust.93.sdf, clust.94
dist HB.94, poly, clust.94, mode=2
load _TEMP.clust.94.sdf, clust.95
dist HB.95, poly, clust.95, mode=2
load _TEMP.clust.95.sdf, clust.96
dist HB.96, poly, clust.96, mode=2
load _TEMP.clust.96.sdf, clust.97
dist HB.97, poly, clust.97, mode=2
load _TEMP.clust.97.sdf, clust.98
dist HB.98, poly, clust.98, mode=2
load _TEMP.clust.98.sdf, clust.99
dist HB.99, poly, clust.99, mode=2
load _TEMP.clust.99.sdf, clust.100
dist HB.100, poly, clust.100, mode=2
load _TEMP.clust.100.sdf, clust.101
dist HB.101, poly, clust.101, mode=2
load _TEMP.clust.101.sdf, clust.102
dist HB.102, poly, clust.102, mode=2
load _TEMP.clust.102.sdf, clust.103
dist HB.103, poly, clust.103, mode=2
load _TEMP.clust.103.sdf, clust.104
dist HB.104, poly, clust.104, mode=2
load _TEMP.clust.104.sdf, clust.105
dist HB.105, poly, clust.105, mode=2
load _TEMP.clust.105.sdf, clust.106
dist HB.106, poly, clust.106, mode=2
load _TEMP.clust.106.sdf, clust.107
dist HB.107, poly, clust.107, mode=2
load _TEMP.clust.107.sdf, clust.108
dist HB.108, poly, clust.108, mode=2
load _TEMP.clust.108.sdf, clust.109
dist HB.109, poly, clust.109, mode=2
load _TEMP.clust.109.sdf, clust.110
dist HB.110, poly, clust.110, mode=2
load _TEMP.clust.110.sdf, clust.111
dist HB.111, poly, clust.111, mode=2
load _TEMP.clust.111.sdf, clust.112
dist HB.112, poly, clust.112, mode=2
load _TEMP.clust.112.sdf, clust.113
dist HB.113, poly, clust.113, mode=2
load _TEMP.clust.113.sdf, clust.114
dist HB.114, poly, clust.114, mode=2
load _TEMP.clust.114.sdf, clust.115
dist HB.115, poly, clust.115, mode=2
load _TEMP.clust.115.sdf, clust.116
dist HB.116, poly, clust.116, mode=2
load _TEMP.clust.116.sdf, clust.117
dist HB.117, poly, clust.117, mode=2
load _TEMP.clust.117.sdf, clust.118
dist HB.118, poly, clust.118, mode=2
load _TEMP.clust.118.sdf, clust.119
dist HB.119, poly, clust.119, mode=2
load _TEMP.clust.119.sdf, clust.120
dist HB.120, poly, clust.120, mode=2
load _TEMP.clust.120.sdf, clust.121
dist HB.121, poly, clust.121, mode=2
load _TEMP.clust.121.sdf, clust.122
dist HB.122, poly, clust.122, mode=2
load _TEMP.clust.122.sdf, clust.123
dist HB.123, poly, clust.123, mode=2
load _TEMP.clust.123.sdf, clust.124
dist HB.124, poly, clust.124, mode=2
load _TEMP.clust.124.sdf, clust.125
dist HB.125, poly, clust.125, mode=2
load _TEMP.clust.125.sdf, clust.126
dist HB.126, poly, clust.126, mode=2
load _TEMP.clust.126.sdf, clust.127
dist HB.127, poly, clust.127, mode=2
load _TEMP.clust.127.sdf, clust.128
dist HB.128, poly, clust.128, mode=2
load _TEMP.clust.128.sdf, clust.129
dist HB.129, poly, clust.129, mode=2
load _TEMP.clust.129.sdf, clust.130
dist HB.130, poly, clust.130, mode=2
load _TEMP.clust.130.sdf, clust.131
dist HB.131, poly, clust.131, mode=2
load _TEMP.clust.131.sdf, clust.132
dist HB.132, poly, clust.132, mode=2
load _TEMP.clust.132.sdf, clust.133
dist HB.133, poly, clust.133, mode=2
load _TEMP.clust.133.sdf, clust.134
dist HB.134, poly, clust.134, mode=2
load _TEMP.clust.134.sdf, clust.135
dist HB.135, poly, clust.135, mode=2
load _TEMP.clust.135.sdf, clust.136
dist HB.136, poly, clust.136, mode=2
load _TEMP.clust.136.sdf, clust.137
dist HB.137, poly, clust.137, mode=2
load _TEMP.clust.137.sdf, clust.138
dist HB.138, poly, clust.138, mode=2
load _TEMP.clust.138.sdf, clust.139
dist HB.139, poly, clust.139, mode=2
load _TEMP.clust.139.sdf, clust.140
dist HB.140, poly, clust.140, mode=2
load _TEMP.clust.140.sdf, clust.141
dist HB.141, poly, clust.141, mode=2
load _TEMP.clust.141.sdf, clust.142
dist HB.142, poly, clust.142, mode=2
load _TEMP.clust.142.sdf, clust.143
dist HB.143, poly, clust.143, mode=2
load _TEMP.clust.143.sdf, clust.144
dist HB.144, poly, clust.144, mode=2
load _TEMP.clust.144.sdf, clust.145
dist HB.145, poly, clust.145, mode=2
load _TEMP.clust.145.sdf, clust.146
dist HB.146, poly, clust.146, mode=2
load _TEMP.clust.146.sdf, clust.147
dist HB.147, poly, clust.147, mode=2
load _TEMP.clust.147.sdf, clust.148
dist HB.148, poly, clust.148, mode=2
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
color cyan, 5ul9_prep and org
color cyan, ref_lig
util.cnc
set ray_trace_mode, 1
set ray_trace_gain, 0.008
set ray_trace_color, black
set pse_export_version, 1.70
save single_VS.summary.sch_top1k.clust-04.pse
quit
