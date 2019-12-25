#!/usr/bin/env python3

import os,re
import html
import tabulate
#from PIL import Image
#import cairo

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import DrawingOptions
from rdkit.Chem.rdmolops import AssignStereochemistryFrom3D


#######################################################################
##
##  Peter M.U. Ung @ MSSM
##
##  v1  14.?
##  v2  18.10.31  improve image quality with svg intermediate
##  v3  19.04.04  deprecate HTML.Table() and use 'tabulate' for formatting
##  v4  19.05.16  improve fidelity of the drawing chiral molecules
##
## Print out a HTML page, in which every row has a maximum of 5 compound png.
## Every major cluster of compounds is grouped together.
## List the Name of the compound, then the Rank and Score.

## SDF Molecules will be formatted into
## [img_Name, Name, Rank, Score, Type]

## Type = sdf		- simple sdf imported by rdkit from source
## Type = smi		- must have header line
##	= formated	- sdf formatted already
##
########################################################################

def grid_print( output_name, Mol_List, Type, column=5 ):

  Data   = []

  ## If input type is simple list of sdf mol, just get the info from sdf
  if Type == 'sdf':
    Img  = []
    for m in Mol_List:
      line = m.GetProp('_Name')
      if re.search(r'::', line):
        I = line.split('::')
      else:
        try:
          rank = m.GetProp('Rank')
        except KeyError:
          rank = '0'
        try:
          score = m.GetProp('Score')
        except KeyError:
          score = '0.0'
        I = [line, rank, score, '']

      svg_name = '_TEMP.{0}.svg'.format(I[0])

      AssignStereochemistryFrom3D(m)
      rdMolDraw2D.PrepareMolForDrawing(m)
      m = Chem.RemoveHs(m)
      AllChem.Compute2DCoords(m)
      DrawingOptions.atomLabelFontSize=18

      Draw.MolToFile(m, svg_name, size=(225,225))
#      cairosvg.svg2png( url=svg_name, write_to=png_name, dpi=240 )
      img_link = '<img src="{0}">'.format(svg_name)

      Img.append([img_link, I[0], I[1],"%.1f" % float(I[2]), I[3]]) 
                # [img_Name, Name, Rank, Score, Type]
    Data.append(Img)
  ## If input type is SMILES, which has no additional info, ignore info
  if Type == 'smi':
    Img  = []
    for m in Mol_List:
      line = m.GetProp('_Name')
#      mol_name = line.split('::')[0]
      mol_name = line

      svg_name = '_TEMP.{0}.svg'.format(mol_name)

      rdMolDraw2D.PrepareMolForDrawing(m)
      m = Chem.RemoveHs(m)
      AllChem.Compute2DCoords(m)
      DrawingOptions.atomLabelFontSize=18

      Draw.MolToFile(m, svg_name, size=(225,225))
#      cairosvg.svg2png( url=svg_name, write_to=png_name, dpi=240 )
      img_link = '<img src="{0}">'.format(svg_name)

      Img.append([img_link, mol_name, '', '', '']) 
    Data.append(Img)
  
  ## if input is already formatted for grid generation, skip all
  if Type == 'formatted':
    Data = Mol_List


##########################################################################
## Writing to HTML
  PAGE = open(output_name+'.html', 'w')
  for idx, C in enumerate(Data):
#    table  = HTML.Table()
    table  = []
    c_temp = []   # Compound Row
    i_temp = []   # ID number Row
    s_temp = []   # Rank/Score Row
    for num, img in enumerate(C):
      if num == 0:
        c_temp.append(img[0])
        i_temp.append('<center>'+img[1]+'</center>')
        s_temp.append('<center>{0}: {1} | {2}</center>'.format(img[4],img[2],img[3]))
      else:
        if num % column != 0:
          c_temp.append(img[0])
          i_temp.append('<center>'+img[1]+'</center>')
          s_temp.append('<center>{0}: {1} | {2}</center>'.format(img[4],img[2],img[3]))
        else:
          table.append(c_temp)
          table.append(i_temp)
          table.append(s_temp)
          c_temp = [img[0]]
          i_temp = ['<center>'+img[1]+'</center>']
          s_temp = ['<center>{0}: {1} | {2}</center>'.format(img[4],img[2],img[3])]

        if num == len(C)-1:
          table.append(c_temp)
          table.append(i_temp)
          table.append(s_temp)
    
#    htmlcode = str(table)
    htmlcode = tabulate.tabulate(table, tablefmt='html')
    if idx == len(Data)-1:
      print(PAGE.write('<p><b>  ###  No Cluster  ###</b></p>'))
    print(PAGE.write("<p><b>  ## Cluster: {0} -- {1} Hits ##</b></p>".format(idx+1, len(C))))

    print(PAGE.write(htmlcode))
  print(PAGE.write("<p><b>  ## no. of cluster: {0} ##</b></p>".format(len(Data))))
  PAGE.close()

  os.system('wkhtmltopdf {0}.html {0}.pdf'.format(output_name))
#  os.system("google-chrome "+output_name+".html")
  os.system('rm _TEMP*.svg {0}.html'.format(output_name))

##########################################################################
## Generate 2D molecule images and insert the images into Cluster Table 
## for PDF file generation
def GenClustTable( Mol_List, output_name, column=5 ):
  Img_Data = []
  for idx, Mols in enumerate(Mol_List):
    Img = []

    for mol in Mols:
     # Get molecule info
      try: m1 = mol.GetProp('Name')
      except KeyError:
        name = mol.GetProp('_Name')
        if re.search(r'::', name):
          m1 = name.split('::')[0]
        else:
          m1 = name
      try: m2 = mol.GetProp('Rank')
      except KeyError:
        m2 = ''
      try: m3 = mol.GetProp('Score')
      except KeyError:
        m3 = ''
      try: m4 = mol.GetProp('Type')
      except KeyError:
        m4 = ''

     # Create tag and write out to sdf file
      mol.SetProp('Cluster', str(idx+1))
      mol.SetProp('SMILES' , Chem.MolToSmiles(mol))
      AssignStereochemistryFrom3D(mol)

      svg_name = '_TEMP.'+m1+'.svg'
      mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
      rdMolDraw2D.PrepareMolForDrawing(mol)
      mol = Chem.RemoveHs(mol)
      AllChem.Compute2DCoords(mol)
      DrawingOptions.atomLabelFontSize=18

      Draw.MolToFile(mol, svg_name, size=(225,225))
#      cairosvg.svg2png( url=svg_name, write_to=png_name, dpi=240 )
      img_link = '<img src="'+svg_name+'">'
          # Img = (image_link, Name, Rank, Score, Type)
      Img.append([img_link, m1, m2, m3, m4])

    Img_Data.append(Img)

## Print out a HTML page, in which every row has a maximum of 5 compound png.
## Every major cluster of compounds is grouped together.
## List the Name of the compound, then the Rank and Score.

  grid_print(output_name, Img_Data, 'formatted', column=5)

  os.system('rm ./_TEMP.*.svg ./{0}.html'.format(output_name))


#######################################################################
## Takes in a molecular list that has list of lists (clusters of molecules)
## and generate a HTML-based table
def GenPyMOLClust(Mol_List, output_name, ref_pdb, sdf=False):
  pymol_pml = open(output_name+".pml", 'w')

  if sdf: m_out = Chem.SDWriter(output_name+'.sdf')

  ref_name = ref_pdb.split('/')[-1].split('.pdb')[0]

  pymol_pml.write("load "+ref_pdb+"\nshow cartoon, poly\nhide lines\ncolor white, poly\ncolor cyan, org\nshow sticks, org and not resn NMA+ACE\n")

  pymol_pml.write("set_bond stick_radius, .15, "+ref_name+" and org\n")
  pymol_pml.write("create ref_lig, "+ref_name+" and org and not resn NMA+ACE\n")
  pymol_pml.write("show lines, byres poly within 5 of ref_lig\n")
  pymol_pml.write("hide sticks, "+ref_name+" and org\n")

  for idx, Mols in enumerate(Mol_List):
    pse_sdf = Chem.SDWriter('_TEMP.clust.{0}.sdf'.format(idx))
    for mol in Mols:
      pse_sdf.write(mol)
      if sdf: m_out.write(mol)

    pse_sdf.flush()
    pse_sdf.close()
    pymol_pml.write("load _TEMP.clust."+str(idx)+".sdf, clust."+str(idx+1)+"\n")
    pymol_pml.write("dist HB."+str(idx+1)+", poly, clust."+str(idx+1)+", mode=2\n")
  pymol_pml.write("show sticks, org\ncenter org\nzoom org\n")
  pymol_pml.write("hide (h. and (e. c extend 1))\n")
  pymol_pml.write("util.cbas\n")
  pymol_pml.write("set mesh_width, 0.1\n")
  pymol_pml.write("set light_count, 1\nset ray_opaque_background, off\n")
  pymol_pml.write("color white, poly\ncolor cyan, "+ref_name+" and org\n")
  pymol_pml.write("set valence\nset dash_gap, 0.25\nhide labels\n")
  pymol_pml.write("set ray_trace_mode, 1\nset ray_trace_gain, .008\n")
  pymol_pml.write("set ray_trace_color, black\n")
  pymol_pml.write('set pse_export_version, 1.70\n')
  pymol_pml.write("save "+output_name+".pse\nquit\n")
  pymol_pml.close()
  if sdf:
    m_out.flush()
    m_out.close()

  os.system("pymol -c -q -Q {0}.pml".format(output_name))
  os.system("rm ./_TEMP.clust.*.sdf")
  os.system("rm {0}.pml".format(output_name))

##########################################################################

