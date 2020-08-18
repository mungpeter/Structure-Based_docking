#
# AutoDockVina.py
#


from Tkinter import *
from tkFileDialog import *
from pymol import cmd, stored, util
#from pymol.cgo import *
import tkSimpleDialog
import tkMessageBox
import os
import sys
import math
import subprocess
import shutil
import fnmatch

import Tkinter,Pmw
import Pmw

g_ligand_protein = [0 for i in range(100000)]
g_resname = [[] for i in range(100000)]
g_net_charge = [0 for i in range(100000)]
g_library_name = ""
g_protein_name_slide = ""
g_ligand_name = [[] for i in range(10000)]




#######################################################################################
def __init__(self):
	cmd.set("retain_order", 1)
	cmd.set("pdb_use_ter_records", 1)
	read_settings(self)


	# Simply add the menu entry and callback
	self.menuBar.addmenu('AutoDockVina', 'Run docking with AutoDockVina',tearoff=TRUE)
	self.menuBar.addmenuitem('AutoDockVina', 'command', 'Export ligand library',
														label = 'Export ligand library',command = lambda s=self : export_ligands(s))

	self.menuBar.addmenuitem('AutoDockVina', 'separator')
       
	self.menuBar.addmenuitem('AutoDockVina', 'command', 'Prepare target',
														label = 'Prepare system and start AutoDock',command = lambda s=self : prepare_protein(s))
    
	self.menuBar.addmenuitem('AutoDockVina', 'separator')

	self.menuBar.addmenuitem('AutoDockVina', 'command', 'Import/Monitor Results',
														label = 'Import Results',command = lambda s=self : import_results(s))

#	self.menuBar.addmenuitem('AutoDockVina', 'command', 'Translate Results into Alignment for RAPTOR',
														#label = 'Translate Results into Alignment for RAPTOR',command = lambda s=self : results_to_raptor(s))
    
	self.menuBar.addmenuitem('AutoDockVina', 'separator')

	self.menuBar.addmenuitem('AutoDockVina', 'command', 'Modify Settings_Linux.txt',
														label = 'Modify Settings_Linux.txt file',command = lambda s=self : modify_settings(s))

#######################################################################################
class ligand_select_dialog(tkSimpleDialog.Dialog):
    def body(self, master):
        self.ok_flag = 0

        Label(master, text="The following objects contain residues with names matching amino acids (or certain cofactors).\nWhich of the the following objects are ligands?").grid(row=0, columnspan=1, sticky=N)

        self.cb_aa = [[] for i in range(10000)]
        self.v_aa = []
        cj = 0
	for na in cmd.get_names("objects"):
            cmd.select('pro', na)
            cmd.select('pro2', 'pro and resn ALA+ARG+ASH+ASN+ASP+CYS+CYX+CY1+GLH+GLN+GLU+GLY+HIS+HIE+HID+HIP+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL+HEM+HEO+NDP+NAP')
            L = cmd.count_atoms('pro2')
            if L > 0:
                var = IntVar()
		self.cb_aa[cj] = Checkbutton(master, text=na, variable=var)
		self.cb_aa[cj].grid(row=cj+1, column=0, sticky=W, padx=5)
                self.v_aa.append(var)
                self.cb_aa[cj].deselect()
                cj += 1
            cmd.delete("pro2")
            cmd.delete("pro")
        
    def validate(self):
        try:
            return 1

        except ValueError:
            tkMessageBox.showwarning(
                "Bad input",
                "Illegal values, please try again"
            )
            return 0

    def apply(self):
        self.ok_flag = 1

    def cancel(self, event=None):
        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

#######################################################################################
class namesDialog(tkSimpleDialog.Dialog):
	def body(self, master):
		self.var_check = []

		self.ok_flag = 0

		self.e1 = [[] for i in range(100000)]
                self.e2 = [[] for i in range(100000)]
		Label(master, text="Please specify/check name and net charge of ligands").grid(row=0, column=0, sticky=W)

                Label(master, text="Ligand name ").grid(row=1, column=0, sticky=W)
                Label(master, text="Net charge of ligand ").grid(row=1, column=1, sticky=W)

#                frame = Frame(master)
#                scrollbar = Scrollbar(frame, orient=VERTICAL)
# 	        self.container_list = Listbox(frame, yscrollcommand=scrollbar.set)
# 	        scrollbar.config(command=self.container_list.yview)
# 	        scrollbar.pack(side=RIGHT, fill=Y)


		cj = 0
                ci = 0
                cm = 0
                cc = 0
                
		for na in cmd.get_names("objects"):
                    g_resname.insert(ci, na)
                    if g_ligand_protein[ci] == 1:
                        if cm >= 25:
                            cc += 1
                            cm = 0
			self.e1[cj] = Entry(master, width=12)
			self.e1[cj].insert(0, g_resname[ci])
			self.e1[cj].grid(row=cm+11, column=3*cc+0)

                        stored.net_charge = 0
                        cmd.iterate(na, "stored.net_charge = stored.net_charge + partial_charge")
                        netch = int(math.floor(stored.net_charge + 0.05))
                        
                        self.e2[cj] = Entry(master, width=5)
                        self.e2[cj].insert(0, netch)
                        self.e2[cj].grid(row=cm+11, column=3*cc+1)

			cj += 1
                        cm += 1

                    ci += 1

		return self.e1[0] # initial focus

	def validate(self):
 		try:
                    cj = 0
                    ci = 0
                    for na in cmd.get_names("objects"):
                        if g_ligand_protein[ci] == 1:
                            g_net_charge[cj] = int(self.e2[cj].get())
                            cj += 1
                        ci += 1
			return 1

		except ValueError:
			tkMessageBox.showwarning(
				"Bad input",
				"Illegal values, please try again"
			)
			return 0

	def apply(self):
		ci = 0
                cj = 0
		for na in cmd.get_names("objects"):
                    if g_ligand_protein[ci] == 1:
			g_resname.insert(ci, self.e1[cj].get())
                        g_net_charge[ci] = int(self.e2[cj].get())
			cj += 1
                    ci += 1

		self.ok_flag = 1

	def cancel(self, event=None):
		# put focus back to the parent window
		self.parent.focus_set()
		self.destroy()
                
#######################################################################################
class clusterTrajDialog(tkSimpleDialog.Dialog):
	def body(self, master):
		self.ok_flag = 0

		Label(master, text="Please specify criteria for clustering different trajectory frames").grid(row=0, column=0, sticky=W)

                Label(master, text="Maximum RMSD between center of cluster and other members [A]").grid(row=1, column=0, sticky=W)
                self.e1 = Entry(master, width=12)
                self.e1.insert(0, "0.75")
                self.e1.grid(row=1, column=1)
                Label(master, text="Minimum number of frames in cluster ").grid(row=2, column=0, sticky=W)
                self.e2 = Entry(master, width=12)
                self.e2.insert(0, "1")
                self.e2.grid(row=2, column=1)
                Label(master, text="Aim for number of clusters [allowed range: 0.5x - 1.5x value] ").grid(row=3, column=0, sticky=W)
                self.e3 = Entry(master, width=12)
                self.e3.insert(0, "50")
                self.e3.grid(row=3, column=1)

		return self.e1 # initial focus

	def validate(self):
 		try:
			return 1

		except ValueError:
			tkMessageBox.showwarning(
				"Bad input",
				"Illegal values, please try again"
			)
			return 0

	def apply(self):
                self.rmsd_max = float(self.e1.get())
                self.min_in_cluster = int(self.e2.get())
                self.aim_cluster_num = int(self.e3.get())
		self.ok_flag = 1

	def cancel(self, event=None):
		# put focus back to the parent window
		self.parent.focus_set()
		self.destroy()
                
#######################################################################################
class clusterDialog(tkSimpleDialog.Dialog):
	def body(self, master):
		self.ok_flag = 0

		Label(master, text="Please specify criteria for clustering docking poses in different trajectory snapshots").grid(row=0, column=0, sticky=W)

		self.v_clust_yes = IntVar()
		self.cb_clust_yes = Checkbutton(master, text="Perform clustering on solutions", variable=self.v_clust_yes)
		self.cb_clust_yes.grid(row=2, column=0, sticky=W, padx=5)
		self.v_clust_yes.set(1)

                Label(master, text="Maximum RMSD between center of cluster and other members [A]").grid(row=11, column=0, sticky=W)
                self.e1 = Entry(master, width=12)
                self.e1.insert(0, "1.5")
                self.e1.grid(row=11, column=1)
                Label(master, text="Minimum number of poses in cluster ").grid(row=12, column=0, sticky=W)
                self.e2 = Entry(master, width=12)
                self.e2.insert(0, "3")
                self.e2.grid(row=12, column=1)

		return self.e1 # initial focus

	def validate(self):
 		try:
			return 1

		except ValueError:
			tkMessageBox.showwarning(
				"Bad input",
				"Illegal values, please try again"
			)
			return 0

	def apply(self):
                self.rmsd_max = float(self.e1.get())
                self.min_in_cluster = int(self.e2.get())
		self.ok_flag = 1

	def cancel(self, event=None):
		# put focus back to the parent window
		self.parent.focus_set()
		self.destroy()
                
#######################################################################################
class detailsDialog1:
        def __init__(self, parent, numlig):
                top = self.top = Toplevel(parent)
        
		self.ok_flag = 0

                cm = order2[numlig]
                
                Label(top, text="Detailed results for ligand %s" % (g_ligand_name[cm]), font=("Helvetica", "16")).grid(row=0, columnspan=3, sticky=N)
                Label(top, text="  ").grid(row=1, columnspan=6, sticky=N)
                
                Label(top, text="Conformation        ").grid(row=10, column=0, sticky=N)
                Label(top, text="Score ").grid(row=10, column=1, sticky=N)
                               
                ci = 0
		for na in range(res_num_poses):
                    if na < num_solutions[cm]:
                        Label(top, text="%s %d" % (sol_name[cm][ci], ci+1)).grid(row=ci+11, column=0, sticky=W)
                        Label(top, text="%.2f" % (sol_affinity[cm][ci])).grid(row=ci+11, column=1, sticky=E)
                        
                        ci += 1


                
                b = Button(top, text="Close", padx=5, command=self.ok)        
                b.grid(row=num_ligands+30, column=1, sticky=S)


        def ok(self):
                self.top.destroy()

#######################################################################################
class detailsDialog2:
        def __init__(self, parent, numlig):
                top = self.top = Toplevel(parent)
                self.numlig = numlig
        
		self.ok_flag = 0

                cm = order2[numlig]
                self.cm = cm
                
                Label(top, text="Detailed results for ligand %s" % (g_ligand_name[cm]), font=("Helvetica", "16")).grid(row=0, columnspan=3, sticky=N)
                Label(top, text="  ").grid(row=1, columnspan=6, sticky=N)
                
                Label(top, text="Cluster representative   ").grid(row=9, column=0, sticky=N)
                Label(top, text="Mean score (+/- stdv)").grid(row=9, column=1, sticky=N)
                               
                Label(top, text="all frames").grid(row=10, column=0, sticky=W)
                Label(top, text="%.2f (+/- %.2f)" % (mean_affinity[cm], stdv_affinity[cm])).grid(row=10, column=1, sticky=E)
                Button(top, text="Show histogram", padx=3, command=self.showhistogram1).grid(row=10, column=2, sticky=E)

                for na in range(num_clusters[cm]):
                    Label(top, text=sol_name[cm][cluster_number[cm][na]]).grid(row=na+11, column=0, sticky=W)
                    Label(top, text="%.2f (+/- %.2f)" % (cluster_mean_affinity[cm][na], cluster_stdv_affinity[cm][na])).grid(row=na+11, column=1, sticky=E)
                    action = lambda x=na: self.showhistogram2(x)
                    Button(top, text="Show histogram", padx=3, command=action).grid(row=na+11, column=2, sticky=E)


                
                b = Button(top, text="Close", padx=5, command=self.ok)        
                b.grid(row=num_clusters[cm]+30, column=1, sticky=S)


        def showhistogram1(self):
                # minimum bin with occupancy
                for j in range(400):
                    if bin_affinity[self.cm][j] > 0:
                        min_x = j-1
                        break
                if min_x < 0:
                    min_x = 0
                elif min_x > 399:
                    min_x = 399
                # maximum bin with occupancy
                for j in range(399, 0, -1):
                    if bin_affinity[self.cm][j] > 0:
                        max_x = j+1
                        break
                if max_x < 0:
                    max_x = 0
                elif max_x > 399:
                    max_x = 399
                # minimum value occupancy = 0
                f_min_y = 0
                # maximum value occupancy = 0
                f_max_y = 0
                for j in range(400):
                    if bin_affinity[self.cm][j] > f_max_y:
                        f_max_y = bin_affinity[self.cm][j]

                d_hist = histogramDialog(self.top, 1, self.cm, 0, min_x, max_x, f_min_y, f_max_y)

        def showhistogram2(self, numclust):
                # minimum bin with occupancy
                for j in range(400):
                    if bin_affinity[self.cm][j] > 0:
                        min_x = j-1
                        break
                if min_x < 0:
                    min_x = 0
                elif min_x > 399:
                    min_x = 399
                # maximum bin with occupancy
                for j in range(399, 0, -1):
                    if bin_affinity[self.cm][j] > 0:
                        max_x = j+1
                        break
                if max_x < 0:
                    max_x = 0
                elif max_x > 399:
                    max_x = 399
                # minimum value occupancy = 0
                f_min_y = 0
                # maximum value occupancy = 0
                f_max_y = 0
                for j in range(400):
                    if bin_affinity[self.cm][j] > f_max_y:
                        f_max_y = bin_affinity[self.cm][j]

                d_hist = histogramDialog(self.top, 2, self.cm, numclust, min_x, max_x, f_min_y, f_max_y)

        def ok(self):
                self.top.destroy()

#######################################################################################
class histogramDialog:
        def __init__(self, parent, type_hist, lig_i, numcluster, min_x, max_x, f_min_y, f_max_y):
                top = self.top = Toplevel(parent)
                
                wid = 400
                hei = 400
                xmin = 50
                xmax = wid-xmin
                ymin = 50
                ymax = hei-ymin
                
                f_min_x = min_x*0.25 - 95.0 
                f_max_x = max_x*0.25 - 95.0 
                min_x_int = int(math.floor(f_min_x))
                max_x_int = int(math.ceil(f_max_x))
                step_x = 1
                if max_x_int - min_x_int > 500:
                    step_x = 100
                elif max_x_int - min_x_int > 200:
                    step_x = 50
                elif max_x_int - min_x_int > 100:
                    step_x = 20
                elif max_x_int - min_x_int > 50:
                    step_x = 10
                elif max_x_int - min_x_int > 20:
                    step_x = 5
                elif max_x_int - min_x_int > 10:
                    step_x = 2

                xlabels = []
                for i in range(min_x_int, max_x_int+1, step_x):
                    xlabels.append(i)

                min_y_int = int(math.floor(f_min_y))
                max_y_int = int(math.ceil(f_max_y))
                step_y = 1
                if max_y_int - min_y_int > 500:
                    step_y = 100
                elif max_y_int - min_y_int > 200:
                    step_y = 50
                elif max_y_int - min_y_int > 100:
                    step_y = 20
                elif max_y_int - min_y_int > 50:
                    step_y = 10
                elif max_y_int - min_y_int > 20:
                    step_y = 5
                elif max_y_int - min_y_int > 10:
                    step_y = 2

                ylabels = []
                for i in range(min_y_int, max_y_int+1, step_y):
                    ylabels.append(i)
                    
                # points    
                xpoints = []
                ypoints = []
                # first point
                f_x = min_x*0.25 - 95.0 
                xpoints.append(f_x)
                ypoints.append(0.0)
                for j in range(400):
                    if j >= min_x and j <= max_x:
                        f_x = j*0.25 - 95.0 
                        xpoints.append(f_x)
                        if type_hist == 1:
                            ypoints.append(bin_affinity[lig_i][j])
                        else:
                            ypoints.append(cluster_bin_affinity[lig_i][numcluster][j])
                # last point
                f_x = max_x*0.25 - 95.0 
                xpoints.append(f_x)
                ypoints.append(0.0)
                
                #=========================================================================================
                # canvas
                canvas = Canvas(top, width=wid, height=hei, bg='white')
                canvas.grid(row=1, column=0, sticky=S)

                # x-axis + labels
                nextspot = xmin
		for label in xlabels:
		    canvas.create_line((nextspot, ymax+5,nextspot, ymax-5),fill="black",width=2)
		    canvas.create_text(nextspot, ymax+15, text=label)
		    if len(xlabels) == 1:
			nextspot = xmax
		    else:
		        nextspot = nextspot + (xmax - xmin)/ (len(xlabels) - 1)
                canvas.create_line(xmin, ymax, xmax, ymax)
                
                # y-axis + labels
		nextspot = ymax
    		for label in ylabels:
		    canvas.create_line((xmin+5,nextspot,xmin-5,nextspot),fill="black",width=2)
		    canvas.create_text(xmin-20,nextspot,text=label)
		    if len(ylabels) == 1:
			nextspot = ymin
		    else:
		        nextspot = nextspot - (ymax - ymin)/ (len(ylabels) - 1)
                canvas.create_line(xmin, ymin, xmin, ymax)
                
                # points connected by lines
    		for ci in range(len(xpoints)):
                    label_x = xpoints[ci]
                    label_y = ypoints[ci]
                    xp = xmin + (xmax - xmin)*(label_x - min_x_int)/(max_x_int - min_x_int)
                    yp = ymax - (ymax - ymin)*(label_y - min_y_int)/(max_y_int - min_y_int)
                    if ci > 0:
                        canvas.create_line((xp_old,yp_old,xp,yp),fill="black",width=2)
                    xp_old = xp
                    yp_old = yp
                
                # Create axis label
                x_label_text = (xmin+xmax)/2
                canvas.create_text(x_label_text, ymax+30, text="Score [kcal/mol]")
                canvas.create_text(5, 10, text="Frequency", anchor=W)
                # create title
                if type_hist == 1:
                    canvas.create_text(x_label_text, 10, text="all frames", anchor=W)
                else:
                    canvas.create_text(x_label_text, 10, text=sol_name[lig_i][cluster_number[lig_i][numcluster]], anchor=W)


                b = Button(top, text="Close", padx=5, command=self.ok)        
                b.grid(row=10, column=0, sticky=S)

        def ok(self):
                self.top.destroy()


#######################################################################################
class displayResultsDialog:
        def __init__(self, parent, type):
                top = self.top = Toplevel(parent)
        
		self.ok_flag = 0

#		self.b1 = [[] for i in range(100000)]
#		self.e1 = [[] for i in range(100000)]
#                self.e2 = [[] for i in range(100000)]        
                
                Label(top, text="AutoDock results", font=("Helvetica", "16")).grid(row=0, columnspan=3, sticky=N)
                Label(top, text="  ").grid(row=1, columnspan=6, sticky=N)
                
                Label(top, text="Ligand name        ").grid(row=10, column=0, sticky=N)
                Label(top, text="Binding free energy ").grid(row=10, column=1, sticky=N)
                


                ci = 0
		for na in range(num_ligands):
                    Label(top, text=g_ligand_name[order2[ci]]).grid(row=ci+11, column=0, sticky=W)
                    if mean_affinity[order2[ci]] < 5.0:
                        if type == 1:
                            Label(top, text= "%.2f" % (mean_affinity[order2[ci]])).grid(row=ci+11, column=1, sticky=E)
                            action = lambda x=ci: self.showdetails1(x)
                            Button(top, text="Show details", padx=3, command=action).grid(row=ci+11, column=2, sticky=E)
                        if type == 2:
                            Label(top, text= "%.2f (+/- %.2f)" % (mean_affinity[order2[ci]], stdv_affinity[order2[ci]])).grid(row=ci+11, column=1, sticky=E)
                            action = lambda x=ci: self.showdetails2(x)
                            Button(top, text="Show details", padx=3, command=action).grid(row=ci+11, column=2, sticky=E)
                    else:
                        Label(top, text="No solution found").grid(row=ci+11, column=1, sticky=E)
                        

                    ci += 1
        
                
                b = Button(top, text="Close", padx=5, command=self.ok)        
                b.grid(row=num_ligands+30, column=1, sticky=S)


        def ok(self):
                self.top.destroy()

        def showdetails1(self, numlig):
                d_details = detailsDialog1(self.top, numlig)

        def showdetails2(self, numlig):
                d_details = detailsDialog2(self.top, numlig)

#######################################################################################
class saveLigandDialog(tkSimpleDialog.Dialog):
        def body(self, master):
                self.parent = master

                self.ok_flag = 0
		self.project_dir = ""
		self.base_dir = ""

		Label(master, text="Library name ").grid(row=1, column=0, sticky=W)

		curdir = os.getcwd()

		self.e2 = Entry(master, width=100)
		self.e2.insert(0, "library_name")
		self.e2.grid(row=1, column=1)

		return self.e2 # initial focus

	def validate(self):
 		try:
			return 1

		except ValueError:
			tkMessageBox.showwarning(
				"Bad input",
				"Illegal values, please try again"
			)
			return 0

	def apply(self):
		self.project_dir = ""
		self.project_dir = "%s/%s" % (autodock_dir, self.e2.get())
                self.project_dir_short = "%s" % (self.e2.get())
		self.base_dir = ""
		self.base_dir = "%s" % (autodock_dir)

		self.ok_flag = 1

	def cancel(self, event=None):
		self.parent.focus_set()
		self.destroy()



#######################################################################################
class importDialog(tkSimpleDialog.Dialog):
    def body(self, master):
        self.ok_flag = 0

        Label(master, text="The following ligand libraries are available.\nPlease choose library to import.").pack(fill='x',expand=1,padx=5,pady=5)

        self.v_visual = IntVar()
        self.cb_visual = Checkbutton(master, text="Open ligands in PyMOL", variable=self.v_visual)
        self.cb_visual.pack(fill='x',expand=1,padx=5,pady=5)
        self.v_visual.set(0)

        filenam_autodock  = autodock_dir + "/AUTODOCK_DICTONARY"

        self.rb_top = [[] for i in range(100000)]
        self.v0 = StringVar()

        # overwrite existing duplicate entry
        if os.path.isfile(filenam_autodock):
            fi = open(filenam_autodock, "r")
            ci = 0
            for i in fi:
                tmp1, tmp2 = i.split(None, 1)
                na = str(tmp2.strip())
                self.rb_top[ci] = Radiobutton(master, text=na, variable=self.v0, value=na, indicatoron="false")
                self.rb_top[ci].pack(fill='x',expand=1,padx=5,pady=5)
                ci += 1
        else:
            tkMessageBox.showwarning(
                "No libaries",
                "Cannot find any ligand libraries in folder %s" % (autodock_dir)
            )
            self.ok_flag = 0
            return 0
            


        return self.rb_top[0] # initial focus
        
    def validate(self):
        try:
            return 1

        except ValueError:
            tkMessageBox.showwarning(
                "Bad input",
                "Illegal values, please try again"
            )
            return 0

    def apply(self):
        self.ok_flag = 1
        self.library_name = str(self.v0.get())
#        print self.library_name

    def cancel(self, event=None):
        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

#######################################################################################
class selectProteinDialog(tkSimpleDialog.Dialog):
    def body(self, master):
        self.ok_flag = 0

        Label(master, text="Select protein for SLIDE").grid(row=0, columnspan=2, sticky=N)

        self.rb_top0 = [[] for i in range(10000)]
        # top level radio box to select between ligand only, ligand + protein(fixed), ligand + protein(zone), all
        self.v0 = StringVar()
        ci = 0
        for na in cmd.get_names("objects"):
            cmd.select('pro', na)
            cmd.select('pro2', 'pro and resn ALA+ARG+ASH+ASN+ASP+CYS+CYX+CY1+GLH+GLN+GLU+GLY+HIS+HIE+HID+HIP+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL+HEM+HEO+NDP+NAP')
            L = cmd.count_atoms('pro2')
            if L > 1:
                self.rb_top0[ci] = Radiobutton(master, text=na, variable=self.v0, value=na)
                self.rb_top0[ci].grid(row=ci+1, column=0, sticky=W)
                ci += 1
            cmd.delete("pro2")
            cmd.delete("pro")
            

        self.rb_top0[0].select()
#        self.v0.set(0)
        
    def validate(self):
        try:
            return 1

        except ValueError:
            tkMessageBox.showwarning(
                "Bad input",
                "Illegal values, please try again"
            )
            return 0

    def apply(self):
        self.object_name = self.v0.get()
        self.ok_flag = 1

    def cancel(self, event=None):
        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

#######################################################################################
class prepareSystem4Autodock:
    def __init__(self, top):

        self.dialog = Pmw.Dialog(top,
                                 buttons = ('OK','Cancel'),
                                 defaultbutton = 'OK',
                                 title = 'AutoDockVina Plugin',
                                 command = self.apply)
                                 
        master = self.dialog.interior()
        self.applic = master

        l0 = Tkinter.Label(master, text='Choose target protein and ligand library for AutoDock simulation')
        l0.pack(expand = 1, fill = 'both', padx = 5, pady = 5)

        group1 = Pmw.Group(master, tag_text='Project definition name')

        f1 = Tkinter.Frame(group1.interior())

        l1 = Tkinter.Label(f1, text='Base directory:')
        l1.grid(row=0, column=0, sticky=W, padx = 5, pady=1)

        curdir = os.getcwd()
        self.eproj = Tkinter.Entry(f1, width=80)
        self.eproj.insert(0, curdir)
        self.eproj.grid(row=0, column=1, sticky=W, padx = 5, pady=1)
        self.bproj = Tkinter.Button(f1, text='Browse', command=self.read_base_dir)
        self.bproj.grid(row=0, column=2, sticky=W, padx = 5, pady=1)

        l2 = Tkinter.Label(f1, text='Project subdirectory:')
        l2.grid(row=1, column=0, sticky=W, padx = 5, pady=1)
        
        self.eproj2 = Tkinter.Entry(f1, width=80)
        self.eproj2.insert(0, 'project_name')
        self.eproj2.grid(row=1, column=1, sticky=W, padx = 5, pady=1)

        f1.pack(fill = 'x', padx = 5, pady=1)

        group1.pack(fill = 'both', expand = 0, padx = 5, pady = 5)



        group2 = Pmw.Group(master, tag_text='Ligand library')

        l3 = Tkinter.Label(group2.interior(), text='Path to ligand library:')
        l3.pack(side = LEFT, fill = 'x', expand = 0, padx = 5, pady = 1)

        self.e1 = Tkinter.Entry(group2.interior(), width=80)
        self.e1.insert(0, "Full path to ligand library")
        self.e1.pack(side = LEFT, fill = 'x', expand = 0, padx = 5, pady = 1)
        self.e3 = Tkinter.Button(group2.interior(), text='Search and Import', command=self.searchLibrary)
        self.e3.pack(side = LEFT, fill = 'none', expand = 0, padx = 5, pady = 1)

        group2.pack(fill = 'both', expand = 0, padx = 5, pady = 5)



        group3 = Pmw.Group(master, tag_text='Protein selection')

        self.v0 = Tkinter.StringVar()

        f2 = Tkinter.Frame(group3.interior())
        group3a = Pmw.Group(f2,                                                                   
                            tag_pyclass = Tkinter.Radiobutton,
                            tag_text = 'Protein from file [PDB,MOL2,MOL]',
                            tag_value = 'file_single',
                            tag_variable = self.v0
                            )

        curdir = os.getcwd()
        l4 = Tkinter.Label(group3a.interior(), text='Path to protein file:')
        l4.pack(side = LEFT, fill = 'x', expand = 0, padx = 5, pady = 1)
        self.e1_pdb = Tkinter.Entry(group3a.interior(), width=80)
        self.e1_pdb.insert(0, "Full path to protein file")
        self.e1_pdb.pack(side = LEFT, fill = 'x', expand = 0, padx = 5, pady = 1)
        self.e3_pdb = Tkinter.Button(group3a.interior(), text='Search and Import', command=self.searchPDB)
        self.e3_pdb.pack(side = LEFT, fill = 'none', expand = 0, padx = 5, pady = 1)

        group3a.pack(fill = 'both', expand = 0, padx = 5, pady = 5)

        f2.pack(fill = 'both', expand = 0, padx = 5, pady = 5)


        f3 = Tkinter.Frame(group3.interior())
        group3b = Pmw.Group(f3,                                                                   
                            tag_pyclass = Tkinter.Radiobutton,
                            tag_text = 'File containing multiple protein structures [NMR-PDB]',
                            tag_value = 'file_folder',
                            tag_variable = self.v0
                            )

        l5 = Tkinter.Label(group3b.interior(), text='Path to protein file:')
        l5.pack(side = LEFT, fill = 'x', expand = 0, padx = 5, pady = 1)
        self.e1_folder = Tkinter.Entry(group3b.interior(), width=80)
        self.e1_folder.insert(0, "Full path to protein file")
        self.e1_folder.pack(side = LEFT, fill = 'x', expand = 0, padx = 5, pady = 1)
        self.e3_folder = Tkinter.Button(group3b.interior(), text='Search and Import', command=self.searchFOLDER)
        self.e3_folder.pack(side = LEFT, fill = 'none', expand = 0, padx = 5, pady = 1)

        group3b.pack(fill = 'both', expand = 0, padx = 5, pady = 5)

        f3.pack(fill = 'both', expand = 0, padx = 5, pady = 5)


        f4 = Tkinter.Frame(group3.interior())
        group3c = Pmw.Group(f4,                                                                   
                            tag_pyclass = Tkinter.Radiobutton,
                            tag_text = 'Trajectory from Amber',
                            tag_value = 'file_multiple',
                            tag_variable = self.v0
                            )

        l6 = Tkinter.Label(group3c.interior(), text='Path to topology file:  ')
        l6.grid(row=0, column=0, sticky=W, padx = 5, pady=1)
        self.e1_top = Tkinter.Entry(group3c.interior(), width=80)
        self.e1_top.insert(0, "Full path to .top file")
        self.e1_top.grid(row=0, column=1, sticky=W, padx = 5, pady=1)
        self.e3_top = Tkinter.Button(group3c.interior(), text='Search and Import', command=self.searchTOP)
        self.e3_top.grid(row=0, column=2, sticky=W, padx = 5, pady=1)
        l7 = Tkinter.Label(group3c.interior(), text='Path to trajectory file:')
        l7.grid(row=1, column=0, sticky=W, padx = 5, pady=1)
        self.e1_trj = Tkinter.Entry(group3c.interior(), width=80)
        self.e1_trj.insert(0, "Full path to .trj file")
        self.e1_trj.grid(row=1, column=1, sticky=W, padx = 5, pady=1)
        self.e3_trj = Tkinter.Button(group3c.interior(), text='Search and Import', command=self.searchTRJ)
        self.e3_trj.grid(row=1, column=2, sticky=W, padx = 5, pady=1)

        group3c.pack(fill = 'both', expand = 0, padx = 5, pady = 5)

        f4.pack(fill = 'both', expand = 0, padx = 5, pady = 5)


        f5 = Tkinter.Frame(group3.interior())
        group3d = Pmw.Group(f5,                                                                   
                            tag_pyclass = Tkinter.Radiobutton,
                            tag_text = 'Protein in current PyMol session',
                            tag_value = 'visual',
                            tag_variable = self.v0
                            )

        l7 = Tkinter.Label(group3d.interior(), text='Objects to select from:  ')
        l7.grid(row=0, column=0, sticky=W, padx = 5, pady=1)

        self.rb_sub0 = [[] for i in range(10000)]
        self.v_sub = Tkinter.StringVar()
        ci = 0
        for na in cmd.get_names("objects"):
            cmd.select('pro', na)
            cmd.select('pro2', 'pro and resn ALA+ARG+ASH+ASN+ASP+CYS+CYX+CY1+GLH+GLN+GLU+GLY+HIS+HIE+HID+HIP+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL+HEM+HEO+NDP+NAP')
            L = cmd.count_atoms('pro2')
            if L > 1:
                self.rb_sub0[ci] = Tkinter.Radiobutton(group3d.interior(), text=na, variable=self.v_sub, value=na)
                self.rb_sub0[ci].grid(row=ci, column=1, sticky=W)
                ci += 1
            cmd.delete("pro2")
            cmd.delete("pro")
                   
        if ci > 0:
            self.rb_sub0[0].select()
        else:
            l8 = Tkinter.Label(group3d.interior(), text='No objects in current PyMOL session')
            l8.grid(row=0, column=0, sticky=W, padx = 5, pady=1)


        group3d.pack(fill = 'both', expand = 0, padx = 5, pady = 5)
        f5.pack(fill = 'both', expand = 0, padx = 5, pady = 5)


        if ci > 0:
            self.v0.set('visual')
        else:
            self.v0.set('file_single')


        group3.pack(fill = 'both', expand = 0, padx = 5, pady = 5)


        group4 = Pmw.Group(master, tag_text='Protein preparation')
        self.v_hyd = Tkinter.IntVar()
        self.cb_hyd = Tkinter.Checkbutton(group4.interior(), text="Let autodock change protonation states", variable=self.v_hyd)
        self.cb_hyd.grid(row=0, column=0, sticky=W, padx = 5, pady=1)

        self.v_hyd.set(0)
        group4.pack(fill = 'both', expand = 0, padx = 5, pady = 5)

        
        self.dialog.activate(geometry = 'centerscreenalways')

    def read_base_dir(self):
        base_dir = ""
        base_dir = askdirectory(title='Select base directory', mustexist=1)
        self.eproj.delete(0, END)
        self.eproj.insert(0, base_dir)

    def searchLibrary(self):
        import_ligands(self.applic)
        self.e1.delete(0, END)
        self.e1.insert(0, g_library_name)

    def searchPDB(self):
	ftypes=(('pdb file', '*.pdb'), ('mol2 file', '*.mol2'), ('mol file', '*.mol'), ("All files", "*"))
	indir = os.getcwd()
	self.g_pdb_name = askopenfilename(initialdir=indir, filetypes=ftypes)
	if self.g_pdb_name:
            self.e1_pdb.delete(0, END)
            self.e1_pdb.insert(0, self.g_pdb_name)
            cmd.load(self.g_pdb_name)
        
    def searchFOLDER(self):
        self.folder_PDB = ""
        self.folder_PDB = askdirectory(title="Select directory", mustexist=1)
	if self.folder_PDB:
            self.e1_folder.delete(0, END)
            self.e1_folder.insert(0, self.folder_PDB)

    def searchTOP(self):
	ftypes=(('top file', '*.top'), ("All files", "*"))
	indir = os.getcwd()
	self.g_top_name = askopenfilename(initialdir=indir, filetypes=ftypes)
	if self.g_top_name:
            self.e1_top.delete(0, END)
            self.e1_top.insert(0, self.g_top_name)
#            cmd.load(self.g_top_name)
        
    def searchTRJ(self):
	ftypes=(('trj file', '*.trj'), ("All files", "*"))
	indir = os.getcwd()
	self.g_trj_name = askopenfilename(initialdir=indir, filetypes=ftypes)
	if self.g_trj_name:
            self.e1_trj.delete(0, END)
            self.e1_trj.insert(0, self.g_trj_name)

    def apply(self, result):
        if result == 'OK':
            self.project_dir = ""
            self.project_dir_short = ""
            self.project_dir = "%s/%s" % (self.eproj.get(), self.eproj2.get())
            self.project_dir_short = "%s" % (self.eproj2.get())
            self.object_name = self.v0.get()
            self.library_dir = self.e1.get()
            self.add_hydrogens = self.v_hyd.get()
            self.ok_flag = 1
            self.dialog.deactivate()
            self.dialog.withdraw()
        else:
            self.ok_flag = 0
            self.dialog.deactivate()
            self.dialog.withdraw()


#######################################################################################
class sizeSearchSpaceDialog:
    def __init__(self, top):

        self.dialog = Pmw.Dialog(top,
                                 buttons = ('OK','Cancel'),
                                 defaultbutton = 'OK',
                                 title = 'AutoDockVina Plugin',
                                 command = self.apply)
                                 
        master = self.dialog.interior()

        nb1 = Pmw.NoteBook(master)

        # search volume pages
        search_volume_page = nb1.add('Search volume')

        for na in cmd.get_names("objects"):
            if na.find("preparedProtein") >= 0:
                cmd.select('usersel', na)
            
                model = cmd.get_model('usersel')
                x, y, z = 0,0,0
                min_x, min_y, min_z = 99999,99999,99999
                max_x, max_y, max_z = -99999,-99999,-99999
                for a in model.atom:
                    x += a.coord[0]
                    y += a.coord[1]
                    z += a.coord[2]
                    if a.coord[0] < min_x:
                        min_x = a.coord[0]
                    if a.coord[0] > max_x:
                        max_x = a.coord[0]
                    if a.coord[1] < min_y:
                        min_y = a.coord[1]
                    if a.coord[1] > max_y:
                        max_y = a.coord[1]
                    if a.coord[2] < min_z:
                        min_z = a.coord[2]
                    if a.coord[2] > max_z:
                        max_z = a.coord[2]
                com_x = x/len(model.atom)
                com_y = y/len(model.atom)
                com_z = z/len(model.atom)
                min_x = math.floor(5*math.floor(min_x/5))
                min_y = math.floor(5*math.floor(min_y/5))
                min_z = math.floor(5*math.floor(min_z/5))
                max_x = math.ceil(5*math.ceil(max_x/5))
                max_y = math.ceil(5*math.ceil(max_y/5))
                max_z = math.ceil(5*math.ceil(max_z/5))
                
                cmd.delete('usersel')

        group_center_box = Pmw.Group(search_volume_page, tag_text='Center of box')
        fcb = Tkinter.Frame(group_center_box.interior())

        group_center_xyz = Pmw.Group(fcb, tag_text='Definition based on x,y,z coordinates')
        lcb1 = Tkinter.Label(group_center_xyz.interior(), text='Center (x):')
        lcb1.grid(row=0, column=0, sticky=W, padx = 5, pady=1)
        self.scale_center_x = Tkinter.Scale(group_center_xyz.interior(), orient=HORIZONTAL, length=500, resolution=0.1, from_= min_x, to= max_x, tickinterval=10, command=self.updateBox)
        self.scale_center_x.set(com_x)
        self.scale_center_x.grid(row=0, column=1, sticky=W)
        lcb2 = Tkinter.Label(group_center_xyz.interior(), text='Center (y):')
        lcb2.grid(row=1, column=0, sticky=W, padx = 5, pady=1)
        self.scale_center_y = Tkinter.Scale(group_center_xyz.interior(), orient=HORIZONTAL, length=500, resolution=0.1, from_= min_y, to= max_y, tickinterval=10, command=self.updateBox)
        self.scale_center_y.set(com_y)
        self.scale_center_y.grid(row=1, column=1, sticky=W)
        lcb3 = Tkinter.Label(group_center_xyz.interior(), text='Center (z):')
        lcb3.grid(row=2, column=0, sticky=W, padx = 5, pady=1)
        self.scale_center_z = Tkinter.Scale(group_center_xyz.interior(), orient=HORIZONTAL, length=500, resolution=0.1, from_= min_z, to= max_z, tickinterval=10, command=self.updateBox)
        self.scale_center_z.set(com_z)
        self.scale_center_z.grid(row=2, column=1, sticky=W)
        group_center_xyz.pack(fill = 'both', expand = 0, padx = 5, pady = 5)

        group_center_selection = Pmw.Group(fcb, tag_text='Definition based on center-of-mass of user selection')
        Tkinter.Label(group_center_selection.interior(), text='Selection').grid(row=0, column=0, sticky=W)
        self.esele = Tkinter.Entry(group_center_selection.interior(), width=20)
        self.esele.insert(0, 'sele')
        self.esele.grid(row=0, column=1)
        self.b1 = Tkinter.Button(group_center_selection.interior(), text='Determine new box coordinates', command=self.com_calc)
        self.b1.grid(row=0, column=2)
        group_center_selection.pack(fill = 'both', expand = 0, padx = 5, pady = 5)

        fcb.pack(fill = 'x', padx = 5, pady=1)
        group_center_box.pack(fill = 'both', expand = 0, padx = 5, pady = 5)

        group_size_box = Pmw.Group(search_volume_page, tag_text='Size of box')
        fsb = Tkinter.Frame(group_size_box.interior())

        group_size_xyz = Pmw.Group(fsb, tag_text='Definition based on x,y,z coordinates')
        lsb1 = Tkinter.Label(group_size_xyz.interior(), text='Box length (x):')
        lsb1.grid(row=0, column=0, sticky=W, padx = 5, pady=1)
        self.box_length_x = Tkinter.Scale(group_size_xyz.interior(), orient=HORIZONTAL, length=500, resolution=0.1, from_= 0, to= 100, tickinterval=10, command=self.updateBox)
        self.box_length_x.set(25)
        self.box_length_x.grid(row=0, column=1, sticky=W)
        lsb2 = Tkinter.Label(group_size_xyz.interior(), text='Box length (y):')
        lsb2.grid(row=1, column=0, sticky=W, padx = 5, pady=1)
        self.box_length_y = Tkinter.Scale(group_size_xyz.interior(), orient=HORIZONTAL, length=500, resolution=0.1, from_= 0, to= 100, tickinterval=10, command=self.updateBox)
        self.box_length_y.set(25)
        self.box_length_y.grid(row=1, column=1, sticky=W)
        lsb3 = Tkinter.Label(group_size_xyz.interior(), text='Box length (z):')
        lsb3.grid(row=2, column=0, sticky=W, padx = 5, pady=1)
        self.box_length_z = Tkinter.Scale(group_size_xyz.interior(), orient=HORIZONTAL, length=500, resolution=0.1, from_= 0, to= 100, tickinterval=10, command=self.updateBox)
        self.box_length_z.set(25)
        self.box_length_z.grid(row=2, column=1, sticky=W)
        group_size_xyz.pack(fill = 'both', expand = 0, padx = 5, pady = 5)
        
        group_size_selection = Pmw.Group(fsb, tag_text='Definition based on user selection')
        Tkinter.Label(group_size_selection.interior(), text='Selection').grid(row=0, column=0, sticky=W)
        self.esele2 = Tkinter.Entry(group_size_selection.interior(), width=20)
        self.esele2.insert(0, 'sele')
        self.esele2.grid(row=0, column=1)
        Tkinter.Label(group_size_selection.interior(), text='Radius around selection').grid(row=0, column=2, sticky=W)
        self.erad = Tkinter.Entry(group_size_selection.interior(), width=10)
        self.erad.insert(0, '5.0')
        self.erad.grid(row=0, column=3)
        bbd = Tkinter.Button(group_size_selection.interior(), text='Determine new box dimensions', command=self.size_calc)
        bbd.grid(row=0, column=4)
        group_size_selection.pack(fill = 'both', expand = 0, padx = 5, pady = 5)

        fsb.pack(fill = 'x', padx = 5, pady=1)
        group_size_box.pack(fill = 'both', expand = 0, padx = 5, pady = 5)

        # flexible residues page
        flexible_residues_page = nb1.add('Flexible residues')

        self.v0 = StringVar()
        group_no_flexible_residues = Pmw.Group(flexible_residues_page,                                                                   
                                    tag_pyclass = Tkinter.Radiobutton,
                                    tag_text = 'No flexible residues',
                                    tag_value = 'none',
                                    tag_variable = self.v0
                                    )
        group_no_flexible_residues.pack(fill = 'both', expand = 0, padx = 5, pady = 5)

        group_box_flexible_residues = Pmw.Group(flexible_residues_page,                                                                   
                                    tag_pyclass = Tkinter.Radiobutton,
                                    tag_text = 'All residues fully contained in search volume',
                                    tag_value = 'box',
                                    tag_variable = self.v0
                                    )
        group_box_flexible_residues.pack(fill = 'both', expand = 0, padx = 5, pady = 5)

        group_cap_flexible_residues = Pmw.Group(flexible_residues_page,                                                                   
                                    tag_pyclass = Tkinter.Radiobutton,
                                    tag_text = 'Residues partially contained within sphere around any atom of selection',
                                    tag_value = 'cap',
                                    tag_variable = self.v0
                                    )
        Tkinter.Label(group_cap_flexible_residues.interior(), text='Selection: ').grid(row=0, column=0, sticky=W)
        self.e1_sel = Tkinter.Entry(group_cap_flexible_residues.interior(), width=20)
        self.e1_sel.insert(0, 'sele')
        self.e1_sel.grid(row=0, column=1)
        self.v0.set('none')
        Tkinter.Label(group_cap_flexible_residues.interior(), text='Radius:').grid(row=0, column=2, sticky=N)
        self.e1_rad = Tkinter.Entry(group_cap_flexible_residues.interior(), width=10)
        self.e1_rad.insert(0, '5.0')
        self.e1_rad.grid(row=0, column=3)
        group_cap_flexible_residues.pack(fill = 'both', expand = 0, padx = 5, pady = 5)

        e3_pdb = Tkinter.Button(flexible_residues_page, text='Define and visualize flexible residues', command=self.searchFlexRes)
        e3_pdb.pack(fill = 'x', expand = 1, padx = 5, pady = 5)


        # output setting page
        output_settings_page = nb1.add('Output settings')

        if g_frames == 'visual' or g_frames == 'file_single':
            Tkinter.Label(output_settings_page, text='max. number of binding modes displayed').grid(row=0, column=0, sticky=W)
            self.enbm = Tkinter.Entry(output_settings_page, width=10)
            self.enbm.insert(0, '10')
            self.enbm.grid(row=0, column=1, sticky=W)

            Tkinter.Label(output_settings_page, text='max. energy difference between best binding mode\n and worst one displayed (kcal/mol)').grid(row=1, column=0, sticky=W)
            self.eenergy = Tkinter.Entry(output_settings_page, width=10)
            self.eenergy.insert(0, '5.0')
            self.eenergy.grid(row=1, column=1, sticky=W)
        else:
            Tkinter.Label(output_settings_page, text='For Autodock simulations\non multiple snapshots (= ensemble docking)\nonly one solution per snapshots will be generated').grid(row=0, column=0, sticky=W)
            Tkinter.Label(output_settings_page, text='').grid(row=0, column=0, sticky=W+E+S+N)
        

        nb1.pack(fill='both',expand=1,padx=5,pady=5)
        nb1.setnaturalsize()

        self.dialog.activate(geometry = 'centerscreenalways')

    def updateScale(self):
        self.scale_center_x.set(self.entry_scale_x.get())
    
    def updateBox(self, tmpf):
        box = makeBox(self.scale_center_x.get(), self.scale_center_y.get(), self.scale_center_z.get(), self.box_length_x.get(), self.box_length_y.get(), self.box_length_z.get())
        drawBox(box)

    def com_calc(self):
        cmd.select('tmp_sel', self.esele.get())
        
        model = cmd.get_model('tmp_sel')
        x, y, z = 0,0,0
        for a in model.atom:
            x += a.coord[0]
            y += a.coord[1]
            z += a.coord[2]
        com_x = x/len(model.atom)
        com_y = y/len(model.atom)
        com_z = z/len(model.atom)
        cmd.delete('tmp_sel')
        self.scale_center_x.set(com_x)
        self.scale_center_y.set(com_y)
        self.scale_center_z.set(com_z)

    def size_calc(self):
        cmd.select('tmp_sel', self.esele2.get())
        
        model = cmd.get_model('tmp_sel')
        min_x, min_y, min_z = 99999,99999,99999
        max_x, max_y, max_z = -99999,-99999,-99999
        for a in model.atom:
            if a.coord[0] < min_x:
                min_x = a.coord[0]
            if a.coord[0] > max_x:
                max_x = a.coord[0]
            if a.coord[1] < min_y:
                min_y = a.coord[1]
            if a.coord[1] > max_y:
                max_y = a.coord[1]
            if a.coord[2] < min_z:
                min_z = a.coord[2]
            if a.coord[2] > max_z:
                max_z = a.coord[2]
        siz_x = math.ceil(max_x - min_x + 2.0*float(self.erad.get()))
        siz_y = math.ceil(max_y - min_y + 2.0*float(self.erad.get()))
        siz_z = math.ceil(max_z - min_z + 2.0*float(self.erad.get()))
        cmd.delete('tmp_sel')
        self.box_length_x.set(siz_x)
        self.box_length_y.set(siz_y)
        self.box_length_z.set(siz_z)

    def searchFlexRes(self):
        if self.v0.get() == "cap":
            cmd.select('tmp_sel', self.esele.get())
            radius = float(self.e1_rad.get())
            for na in cmd.get_names("objects"):
                if na.find("preparedProtein") >= 0:
                    cmd.select('tmp_pro', na)
            cmd.select('flexible_residues', 'tmp_pro and byres (tmp_sel around %f)' % radius)
            cmd.delete('tmp_sel')
            cmd.delete('tmp_pro')
        elif self.v0.get() == "box":
            min_box_x = self.scale_center_x.get() - self.box_length_x.get()/2.0
            min_box_y = self.scale_center_y.get() - self.box_length_y.get()/2.0
            min_box_z = self.scale_center_z.get() - self.box_length_z.get()/2.0
            max_box_x = self.scale_center_x.get() + self.box_length_x.get()/2.0
            max_box_y = self.scale_center_y.get() + self.box_length_y.get()/2.0
            max_box_z = self.scale_center_z.get() + self.box_length_z.get()/2.0
#            print min_box_x, min_box_y, min_box_z, max_box_x, max_box_y, max_box_z
            for na in cmd.get_names("objects"):
                if na.find("preparedProtein") >= 0:
                    cmd.select('tmp_sel', na)
            model = cmd.get_model('tmp_sel')
            oldresi = -9999
            flag_in = 0
            ci = 0
            sel_string = ""
            for a in model.atom:
                if a.resi != oldresi:
                    if flag_in == 1:
                        if ci == 0:
                            sel_string += oldresi
                        else:
                            sel_string += "+" + oldresi
                        ci += 1
                    oldresi = a.resi
                    flag_in = 1
                if a.coord[0] < min_box_x or a.coord[1] < min_box_y or a.coord[2] < min_box_z or a.coord[0] > max_box_x or a.coord[1] > max_box_y or a.coord[2] > max_box_z:
                    flag_in = 0 
            if flag_in == 1:
                if ci == 0:
                    sel_string += oldresi
                else:
                    sel_string += "+" + oldresi
                ci += 1
            cmd.select('flexible_residues', 'tmp_sel and resi %s' % sel_string)
            cmd.delete('tmp_sel')
        else:
            pass
#            cmd.select('flexible_residues', '')
                    
    def validate(self):
        try:
            return 1

        except ValueError:
            tkMessageBox.showwarning(
                "Bad input",
                "Illegal values, please try again"
            )
            return 0

    def apply(self, result):
        if result == 'OK':
            self.cb_x = self.scale_center_x.get()
            self.cb_y = self.scale_center_y.get()
            self.cb_z = self.scale_center_z.get()
            self.bl_x = self.box_length_x.get()
            self.bl_y = self.box_length_y.get()
            self.bl_z = self.box_length_z.get()
            self.flex_residues = ""
            if g_frames == "visual" or g_frames == "file_single":
                self.num_binding_modes = int(self.enbm.get())
                self.energy_binding_modes = float(self.eenergy.get())
            else:
                self.num_binding_modes = 1
                self.energy_binding_modes = 1.0
            if self.v0.get() != "none":
                str_flex = ""
                cmd.select('tmp_sel', 'flexible_residues and name CA')
                model = cmd.get_model('tmp_sel')
                ci = 0
                for a in model.atom:
#                    print a.name, a.resn, a.resi
                    if a.resn != "TYR":
                        if ci == 0:
                            str_flex += "%s%s" % (a.resn, a.resi)
                        else:
                            str_flex += "_%s%s" % (a.resn, a.resi)
                        ci += 1
                cmd.delete('tmp_sel')
                self.flex_residues = str_flex
                print self.flex_residues 

            self.ok_flag = 1
            self.dialog.deactivate()
            self.dialog.withdraw()
        else:
            self.ok_flag = 0
            self.dialog.deactivate()
            self.dialog.withdraw()

#######################################################################################
class affinityDialog(tkSimpleDialog.Dialog):
	def body(self, master):
		self.var_check = []

		self.ok_flag = 0

		self.e1 = [[] for i in range(10000)]
		self.e2 = [[] for i in range(10000)]
		self.e3 = [[] for i in range(10000)]
		Label(master, text="Ligand name ").grid(row=0, column=0, sticky=W)
		Label(master, text="Affinity [nanomolar] ").grid(row=0, column=1, sticky=W)

                cm = 0
                cc = 0
                
                for cj in range(num_ligands):
                        
                        if cm >= 25:
                            cc += 1
                            cm = 0
                        Label(master, text="Ligand name ").grid(row=0, column=3*cc+0, sticky=W)
                        Label(master, text="Affinity [nanomolar] ").grid(row=0, column=3*cc+1, sticky=W)

                        self.e1[cj] = Entry(master, width=12)
			self.e1[cj].insert(0, g_ligand_name[cj])
			self.e1[cj].grid(row=cm+1, column=3*cc+0)

			self.e2[cj] = Entry(master, width=12)
			self.e2[cj].insert(0, g_dg[cj])
			self.e2[cj].grid(row=cm+1, column=3*cc+1)

                        cm += 1


		return self.e1[0] # initial focus

	def validate(self):
 		try:
		#			self.resname= self.e1.get()
			return 1

		except ValueError:
			tkMessageBox.showwarning(
				"Bad input",
				"Illegal values, please try again"
			)
			return 0

	def apply(self):
                for cj in range(num_ligands):
			g_ligand_name.insert(cj, self.e1[cj].get())
			g_dg.insert(cj, self.e2[cj].get())

		self.ok_flag = 1

	def cancel(self, event=None):
		# put focus back to the parent window
		self.parent.focus_set()
		self.destroy()

#######################################################################################
class numPosesDialog(tkSimpleDialog.Dialog):
	def body(self, master):
		self.var_check = []

		self.ok_flag = 0
                
                Label(master, text="Maximum number of top-ranked poses used for alignment").grid(row=0, column=0, sticky=W)
                self.e1 = Entry(master, width=30)
                self.e1.insert(0, res_num_poses)
                self.e1.grid(row=0, column=1)

                Label(master, text="Name of output folder").grid(row=1, column=0, sticky=W)
                self.e2 = Entry(master, width=30)
                self.e2.insert(0, "autodock_results")
                self.e2.grid(row=1, column=1)

		return self.e1 # initial focus

	def validate(self):
 		try:
		#			self.resname= self.e1.get()
			return 1

		except ValueError:
			tkMessageBox.showwarning(
				"Bad input",
				"Illegal values, please try again"
			)
			return 0

	def apply(self):
                self.numPoses4Symposar = self.e1.get()
                self.results_folder = self.e2.get()

		self.ok_flag = 1

	def cancel(self, event=None):
		# put focus back to the parent window
		self.parent.focus_set()
		self.destroy()

#######################################################################################
class testtrainDialog(tkSimpleDialog.Dialog):
	def body(self, master):
		self.ok_flag = 0

		Label(master, text="Select ligands (click on ligands in main display) and press one of the following buttons ").grid(row=0, column=0, sticky=W)

		self.e1 = Button(master, text='Select as training set', command=self.select_train)
		self.e1.grid(row=1, column=0)

		self.e2 = Button(master, text='Select as test set', command=self.select_test)
		self.e2.grid(row=2, column=0)

		return self.e1 # initial focus

	def validate(self):
 		try:
#			self.resname= self.e1.get()
			return 1

		except ValueError:
			tkMessageBox.showwarning(
				"Bad input",
				"Illegal values, please try again"
			)
			return 0

	def apply(self):
		self.ok_flag = 1

	def cancel(self, event=None):
		# put focus back to the parent window
		self.parent.focus_set()
		self.destroy()

	def select_train(self):
		stored.resnam = []
		cmd.select('sele2', 'sele and index 1')
		cmd.iterate ('sele2', "stored.resnam.append(resn)")
		ci = 0
		for na in cmd.get_names("objects"):
                    cmd.select('pro', na)
                    cmd.select('pro2', 'pro and resn ALA+ARG+ASH+ASN+ASP+CYS+CYX+CY1+GLH+GLN+GLU+GLY+HIS+HIE+HID+HIP+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL+HEM+HEO+NDP+NAP')
                    L = cmd.count_atoms('pro2')
                    if L == 0:
			g_training[ci] = 0
			for sa in stored.resnam:
				if na == sa:
					g_training[ci] = 1
#			cmd.select('pro', na)
                    ci = ci + 1
                    cmd.delete("pro2")
                    cmd.delete("pro")
		cmd.delete('sele2')

# color training green, test red
		ci = 0
		for na in cmd.get_names("objects"):
                    cmd.select('pro', na)
                    cmd.select('pro2', 'pro and resn ALA+ARG+ASH+ASN+ASP+CYS+CYX+CY1+GLH+GLN+GLU+GLY+HIS+HIE+HID+HIP+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL+HEM+HEO+NDP+NAP')
                    L = cmd.count_atoms('pro2')
                    if L == 0:
			if g_training[ci] == 0:
				util.cbam(na)
			else:
				util.cbag(na)
                    ci = ci + 1
                    cmd.delete("pro2")
                    cmd.delete("pro")

	def select_test(self):
		stored.resnam = []
		cmd.select('sele2', 'sele and index 1')
		cmd.iterate ('sele2', "stored.resnam.append(resn)")
		ci = 0
		for na in cmd.get_names("objects"):
                    cmd.select('pro', na)
                    cmd.select('pro2', 'pro and resn ALA+ARG+ASH+ASN+ASP+CYS+CYX+CY1+GLH+GLN+GLU+GLY+HIS+HIE+HID+HIP+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL+HEM+HEO+NDP+NAP')
                    L = cmd.count_atoms('pro2')
                    if L == 0:
			g_training[ci] = 1
			for sa in stored.resnam:
				if na == sa:
					g_training[ci] = 0
#			cmd.select('pro', na)
                    ci = ci + 1
                    cmd.delete("pro2")
                    cmd.delete("pro")
		cmd.delete('sele2')

# color training green, test red
		ci = 0
		for na in cmd.get_names("objects"):
                    cmd.select('pro', na)
                    cmd.select('pro2', 'pro and resn ALA+ARG+ASH+ASN+ASP+CYS+CYX+CY1+GLH+GLN+GLU+GLY+HIS+HIE+HID+HIP+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL+HEM+HEO+NDP+NAP')
                    L = cmd.count_atoms('pro2')
                    if L == 0:
			if g_training[ci] == 0:
				util.cbam(na)
			else:
				util.cbag(na)
                    ci = ci + 1
                    cmd.delete("pro2")
                    cmd.delete("pro")


#######################################################################################
def write_mol2(sel1, filename):
# find name and object number of selected molecule
	cmd.select('sele_t', sel1)
	cmd.select('sele2', 'sele_t and index 1')
	stored.resnam = []
	cmd.iterate ('sele2', "stored.resnam.append(resn)")
        ci = 0
        for na in cmd.get_names("objects"):
		for sa in stored.resnam:
			if na == sa:
                            name = ci
                ci += 1
        cmd.delete("sele_t")
        cmd.delete("sele2")

# save pdb and mol file
	cmd.save("%s.pdb" % (filename), sel1, 0, "pdb")
	cmd.save("%s.mol" % (filename), sel1, 0, "mol")
# convert file using babel
        os.system("babel -imol %s.mol -omol2 %s.mol2\n" % (filename, filename))
# change atom names in mol2 file to those in pdb file
	f_pdb = open("%s.pdb" % (filename), "r")
	atnam = []
	resnam = []
	i = 0
	for line in f_pdb:
		if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
                    # Br atoms are written one column to the left by pymol
                    if line[12:13] == "B":
                        atnam.append(line[12:15])
                    elif line[12:13] == "C":
                        atnam.append(line[12:15])
                    else:
			atnam.append(line[13:16])
                    resnam.append(line[17:21])
                    i = i + 1
	f_pdb.close()
        

	f_old = open("%s.mol2" % (filename), "r")
	f_new = open("new.mol2", "w")

        
	flag = 0
	for line in f_old:
		if line[0:13] == '@<TRIPOS>BOND':
			flag = 0
		if flag == 0:
			f_new.write(line)
		else:
			new_line = []
			new_line.append(line[0:8])
			new_line.append(str(atnam.pop(0)))
			new_line.append(line[11:58])
			new_line.append(str(resnam.pop(0)))
			new_line.append(line[62:])
			f_new.write(''.join(new_line))
		if line[0:13] == '@<TRIPOS>ATOM':
			flag = 1

	f_old.close()
	f_new.close()

        os.remove("%s.mol2" % (filename))
        os.rename("new.mol2", "%s.mol2" % (filename))
        os.remove("%s.mol" % (filename))
        os.remove("%s.pdb" % (filename))

#######################################################################################
def export_ligands(app):
	curdir = os.getcwd()

        # store ligands in g_ligand_protein
        ci = 0
	for na in cmd.get_names("objects"):
            cmd.select('pro', na)
            cmd.select('pro2', 'pro and resn ALA+ARG+ASH+ASN+ASP+CYS+CYX+CY1+GLH+GLN+GLU+GLY+HIS+HIE+HID+HIP+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL+HEM+HEO+NDP+NAP')
            L = cmd.count_atoms('pro2')
            if L == 0:
                g_ligand_protein[ci] = 1
            cmd.delete("pro2")
            cmd.delete("pro")
            ci += 1
        

        # check for residue names == amino acids        
        cj = 0
        for na in cmd.get_names("objects"):
            cmd.select('pro', na)
            cmd.select('pro2', 'pro and resn ALA+ARG+ASH+ASN+ASP+CYS+CYX+CY1+GLH+GLN+GLU+GLY+HIS+HIE+HID+HIP+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL+HEM+HEO+NDP+NAP')
            L = cmd.count_atoms('pro2')
            if L > 0:
                cj += 1
            cmd.delete("pro2")
            cmd.delete("pro")
          
        # found ligands with amino acid names
        if cj > 0:
            aa_sel2 = ligand_select_dialog(app.root)
            if aa_sel2.ok_flag == 0:
		return

            ci = 0
            cj = 0
            for na in cmd.get_names("objects"):
                cmd.select('pro', na)
                cmd.select('pro2', 'pro and resn ALA+ARG+ASH+ASN+ASP+CYS+CYX+CY1+GLH+GLN+GLU+GLY+HIS+HIE+HID+HIP+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL+HEM+HEO+NDP+NAP')
                L = cmd.count_atoms('pro2')
                if L > 0:
                    if aa_sel2.v_aa[cj].get() == 1:
                        g_ligand_protein[ci] = 1
                    cj += 1
                ci += 1
        # did only find ligands with name different from amino acid
        else:
            cj = 0
            for na in cmd.get_names("objects"):
                g_ligand_protein[cj] = 1
                cj += 1
            


	z = namesDialog(app.root)

	if z.ok_flag == 0:
		return
        
        # assign names to residues
	ci = 0
        c_lig = 0
	for na in cmd.get_names("objects"):
            if g_ligand_protein[ci] == 1:
		aaa = str(g_resname[ci])
		cmd.alter(na, "resn=\'%s\'" % (aaa))
		cmd.set_name(na, aaa)
                c_lig += 1
            ci += 1

        
        # check if all residues have name with max. length of 3
#        ci = 0
#	for na in cmd.get_names("objects"):
#            if g_ligand_protein[ci] == 1:
#                if len(na) > 3:
#                    tkMessageBox.showwarning(
#				"Bad object name",
#				"The residue/object name of all compounds should be 3 characters long.\nUse 'Change ligand properties' to modify name."
#                    )
#                    return 0
#            ci += 1


        # library name and directory
        np = saveLigandDialog(app.root)
	if np.ok_flag == 0:
		return

	proj_dir = np.project_dir
        proj_dir_short = np.project_dir_short
        base_dir = np.base_dir
        if not os.path.isdir(base_dir):
            os.mkdir(base_dir)

	if os.path.isdir(proj_dir):
            if tkMessageBox.askokcancel("Library exists","Library with name %s already exists\nDo you want to overwrite library?" % (proj_dir)) == 0:
                return 0
            
            for af in os.listdir(proj_dir):
                possible_dir = proj_dir + "/%s" % (af)
                if os.path.isdir(possible_dir):
                    for bf in os.listdir(possible_dir):
                        os.remove("%s/%s" % (possible_dir, bf))
                    os.rmdir(possible_dir)
                else:
                    os.remove("%s/%s" % (proj_dir, af))
            os.rmdir(proj_dir)
        os.mkdir(proj_dir)

        # generate or append library to dictonary
        if not os.path.isdir(autodock_dir):
            os.mkdir(autodock_dir)
        os.chdir(autodock_dir)
        filenam_autodock  = autodock_dir + "/AUTODOCK_DICTONARY"
        if not os.path.isfile(filenam_autodock):
            fi = open(filenam_autodock, "w")
            fi.close()
            

        # overwrite existing duplicate entry
        fout = open("tmp.tmp", "w")
        if os.path.isfile(filenam_autodock):
            fi = open(filenam_autodock, "r")
            for i in fi:
                tmp1, tmp3 = i.split(None, 1)
                tmp2 = tmp3.strip()
                tmp4 = proj_dir_short.strip()
                if not (tmp2.find(tmp4) == 0 and tmp4.find(tmp2) == 0):
                    fout.write(i)
            fi.close()
        fout.write("%d   %s\n" % (c_lig, proj_dir_short.strip()))

                                
        fout.close()
        if os.path.isfile(filenam_autodock):
            os.remove(filenam_autodock)
        os.rename("tmp.tmp", filenam_autodock)
        

        # write files to dictionary
	os.chdir(proj_dir)
        os.mkdir("%s/PDBs" % proj_dir)
#	print proj_dir
#	print os.getcwd()
        

        # generate omega2 output        
        fcom = open("RunCommands", "w")

        os.system("cp %s/prepare_ligand4.py ." % (autodock_exe_dir))
        ci = 0
	for na in cmd.get_names("objects"):
            if g_ligand_protein[ci] == 1:
                cmd.select('pro', na)
                os.chdir("%s/PDBs" % proj_dir)
                cmd.save("%s.pdb" % (na), 'pro', 0, "pdb")
                os.chdir(proj_dir)
                write_mol2('pro', "%s" % (na))
#                fcom.write("pythonsh ./prepare_ligand4.py -l %s.mol2 -A 'hydrogens'\n" % (na))
                fcom.write("pythonsh ./prepare_ligand4.py -l %s.mol2\n" % (na))
                fcom.write("if [ ! -f %s.pdbqt ]; then\n" % na)
#                fcom.write("   ./prepare_ligand4.py -l %s.mol2 -A 'hydrogens'\n" % (na))
                fcom.write("   ./prepare_ligand4.py -l %s.mol2\n" % (na))
                fcom.write("fi\n")
                fcom.write("rm %s.mol2\n" % (na))
            ci += 1
        
        fcom.close()
        # translate mol into mol2 file
	os.system("mkdir tmp_autodock")
	os.system("cp * ./tmp_autodock/")
	os.system("cd ./tmp_autodock; dos2unix RunCommands; chmod a+x RunCommands; ./RunCommands; cd ..")
	os.system("cp ./tmp_autodock/*.pdbqt .")
	os.system("rm -rf tmp_autodock")
#        return 0

        os.remove("RunCommands")
        os.remove("prepare_ligand4.py")

        ci = 0
	for na in cmd.get_names("objects"):
            if g_ligand_protein[ci] == 1:
                os.remove("%s.mol2" % (na))
            ci += 1
            
        print "New library has been prepared."

	os.chdir(curdir)

#######################################################################################
def import_ligands(app):
        global g_library_name
        global g_library_name_short
        
	curdir = os.getcwd()
        ni = importDialog(app)
        if ni.ok_flag == 0:
            return 0
            
        g_library_name = autodock_dir + "/" + ni.library_name
        g_library_name_short = ni.library_name
#        print g_library_name
        
        if ni.v_visual.get() == 1:
#            ci = 0
#            for na in cmd.get_names("objects"):
#                ci += 1
                
#            if ci > 0:
#                if tkMessageBox.askokcancel("Read project", "Reading project will delete all exisiting objects in current session.") == 0:
#                    return 0
#            # remove existing objects
#            for na in cmd.get_names("objects"):
#                cmd.remove(na)
#                cmd.delete(na)

            for i in os.listdir(g_library_name):
                j = g_library_name + "/" + i
                cmd.load(j)

#######################################################################################
def prepare_protein(app):
        global g_protein_name_slide
        global object_name
        global flag_metabol
        global metabol_xyz
        global g_frames
        global num_frames
        
        flag_metabol = 0
        metabol_xyz = [0.0 for i in range(3)]
        
        curdir = os.getcwd()
        
        ci = 0
        for na in cmd.get_names("objects"):
            cmd.select('pro', na)
            cmd.select('pro2', 'pro and resn ALA+ARG+ASH+ASN+ASP+CYS+CYX+CY1+GLH+GLN+GLU+GLY+HIS+HIE+HID+HIP+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL+HEM+HEO+NDP+NAP')
            L = cmd.count_atoms('pro2')
            if L > 1:
                object_name = na
                ci += 1
            cmd.delete("pro2")
            cmd.delete("pro")
        
        nprep = prepareSystem4Autodock(app.root)
        if nprep.ok_flag == 0:
            return 0

	md_outfile = nprep.project_dir + "_client"
        lib_infile = md_outfile + "/in"
        lib_outfile = md_outfile + "/out"
        mult_frames_file = md_outfile + "/frames"
	md_save_project = 1

        # create project directory and delete old files (if existing)
        for root, dirs, files in os.walk(md_outfile, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        if os.path.isdir(md_outfile) == 0:
            os.mkdir(md_outfile)

	curdir = os.getcwd()

        # change into project directory
        os.chdir(md_outfile)

        os.system("cp %s/prepare_receptor4.py %s/prepare_flexreceptor4.py ." % (autodock_exe_dir, autodock_exe_dir))
        os.system("cp -r %s ." % (nprep.library_dir))
        os.system("mv %s in" % g_library_name_short)

        # change into frames directory
        os.mkdir(mult_frames_file)
        os.chdir(mult_frames_file)
        
        g_frames = nprep.object_name

        if nprep.object_name == "visual":
            object_name = nprep.v_sub.get()
            for na in cmd.get_names("objects"):
                if na == object_name:
                    # remove hydrogens from protein
                    cmd.select('pro', na)
                    cmd.save("protein4autodock.pdb", 'pro', 0, "pdb")
                    modifyPDB2standard("protein4autodock.pdb")
                    if nprep.add_hydrogens == 1:
                        os.system("pythonsh ../prepare_receptor4.py -r protein4autodock.pdb -o preparedProtein.pdbqt -A 'checkhydrogens' -U nphs\n")
                    else:
                        os.system("pythonsh ../prepare_receptor4.py -r protein4autodock.pdb -o preparedProtein.pdbqt -A 'None' -U nphs\n")
                    pdbqt2pdb("preparedProtein.pdbqt", "preparedProtein.pdb")
                    cmd.delete("pro")
        elif nprep.object_name == "file_single":
            os.system("cp %s ./protein4autodock.pdb" % nprep.g_pdb_name)
            modifyPDB2standard("protein4autodock.pdb")
            if nprep.add_hydrogens == 1:
                os.system("pythonsh ../prepare_receptor4.py -r protein4autodock.pdb -o preparedProtein.pdbqt -A 'checkhydrogens' -U nphs\n")
            else:
                os.system("pythonsh ../prepare_receptor4.py -r protein4autodock.pdb -o preparedProtein.pdbqt -A 'None' -U nphs\n")
            pdbqt2pdb("preparedProtein.pdbqt", "preparedProtein.pdb")
        elif nprep.object_name == "file_folder":
            # copy protein files
            j = os.listdir(nprep.folder_PDB)
            j.sort()
            for na in j:
                if na.find(".pdb") >= 0 or na.find(".PDB") >= 0:
                    # generate pdbqt file for first frame
                    os.system("cp %s/%s %s/" % (nprep.folder_PDB, na, mult_frames_file))
            # take first file
            j = os.listdir(mult_frames_file)
            j.sort()
            for na in j:
                if na.find(".pdb") >= 0 or na.find(".PDB") >= 0:
                    # generate pdbqt file for first frame
                    if nprep.add_hydrogens == 1:
                        os.system("cd %s; pythonsh ../prepare_receptor4.py -r %s -o preparedProtein.pdbqt -A 'checkhydrogens' -U nphs" % (mult_frames_file, na))
                    else:
                        os.system("cd %s; pythonsh ../prepare_receptor4.py -r %s -o preparedProtein.pdbqt -A 'None' -U nphs" % (mult_frames_file, na))
                    break
            # 
            pdbqt2pdb("preparedProtein.pdbqt", "preparedProtein.pdb")
            cmd.load("preparedProtein.pdb")
            # define search volume
            nsize = sizeSearchSpaceDialog(app.root)
            if nsize.ok_flag == 0:
                return 0
            os.system("rm preparedProtein.pdbqt")
            os.system("rm preparedProtein.pdb")
            
            fcom = open("RunCommands", "w")
            flog = open("Map_pdb_names.txt", "w")
            j = os.listdir(nprep.folder_PDB)
            j.sort()
            ci = 0
            for na in j:
                if na.find(".pdb") >= 0 or na.find(".PDB") >= 0 or na.find(".mol2") >= 0 or na.find(".MOL2") >= 0 or na.find(".mol") >= 0 or na.find(".MOL") >= 0:
                    ktmpa = (ci+1) % 100000
                    k100000 = ((ci+1) - ktmpa)/100000
                    ktmpb = ktmpa % 10000
                    k10000 = (ktmpa - ktmpb)/10000
                    ktmpa = ktmpb % 1000
                    k1000 = (ktmpb - ktmpa)/1000
                    ktmpb = ktmpa % 100
                    k100 = (ktmpa - ktmpb)/100
                    ktmpa = ktmpb % 10
                    k10 = (ktmpb - ktmpa)/10
                    ktmpb = ktmpa % 1
                    k1 = (ktmpa - ktmpb)
                    if nprep.add_hydrogens == 1:
                        fcom.write("pythonsh ../prepare_receptor4.py -r %s -o preparedProtein_%d%d%d%d%d%d.pdbqt -A 'checkhydrogens' -U nphs\n" % (na, k100000, k10000, k1000, k100, k10, k1))
                    else:
                        fcom.write("pythonsh ../prepare_receptor4.py -r %s -o preparedProtein_%d%d%d%d%d%d.pdbqt -A 'None' -U nphs\n" % (na, k100000, k10000, k1000, k100, k10, k1))                        
                    flog.write("%s preparedProtein_%d%d%d%d%d%d.pdbqt\n" % (na, k100000, k10000, k1000, k100, k10, k1))
                    ci += 1
            fcom.close()
            flog.close()
        elif nprep.object_name == "file_multiple":
            os.system("cp %s %s/protein4autodock.top\n" % (nprep.g_top_name, mult_frames_file))
            os.system("cp %s %s/protein4autodock.trj\n" % (nprep.g_trj_name, mult_frames_file))
            # generate first frame for defining search volume
            fcom = open("RunCommands", "w")
            fcom.write("%s/exe/ptraj protein4autodock.top << EOF\n" % (amber_dir))
            fcom.write("trajin protein4autodock.trj 1 1 1\n")
            fcom.write("strip :WAT\n")
            fcom.write("strip \"@Cl-\"\n")
            fcom.write("strip \"@Na+\"\n")
            fcom.write("strip !:ALA,ARG,ASH,ASN,ASP,CYS,CYX,CY1,GLH,GLN,GLU,GLY,HIS,HIE,HID,HIP,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL,HEM,HEO,FAD,NDP,NAP\n")
            fcom.write("rms first @CA,C,N\n")
            fcom.write("trajout protein4autodock.pdb pdb\n")
            fcom.write("EOF\n")
            fcom.close()
            os.system("cd %s; chmod a+x RunCommands; ./RunCommands\n" % mult_frames_file)
            os.system("rm RunCommands\n")
            # generate pdbqt file for first frame
            os.system("mv protein4autodock.pdb.1 tmp.pdb")
            if nprep.add_hydrogens == 1:
                os.system("pythonsh ../prepare_receptor4.py -r tmp.pdb -o preparedProtein.pdbqt -A 'checkhydrogens' -U nphs")
            else:
                os.system("pythonsh ../prepare_receptor4.py -r tmp.pdb -o preparedProtein.pdbqt -A 'None' -U nphs")
            os.system("rm tmp.pdb")
            # 
            pdbqt2pdb("preparedProtein.pdbqt", "preparedProtein.pdb")
            cmd.load("preparedProtein.pdb")
            # define search volume
            nsize = sizeSearchSpaceDialog(app.root)
            if nsize.ok_flag == 0:
                return 0
            os.system("rm preparedProtein.pdbqt")
            os.system("rm preparedProtein.pdb")

            # generate all pdbqt files
            fcom = open("RunCommands", "w")
            fcom.write("%s/exe/ptraj protein4autodock.top << EOF\n" % (amber_dir))
            fcom.write("trajin protein4autodock.trj\n")
            fcom.write("strip :WAT\n")
            fcom.write("strip \"@Cl-\"\n")
            fcom.write("strip \"@Na+\"\n")
            fcom.write("strip !:ALA,ARG,ASH,ASN,ASP,CYS,CYX,CY1,GLH,GLN,GLU,GLY,HIS,HIE,HID,HIP,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL,HEM,HEO,FAD,NDP,NAP\n")
            fcom.write("rms first @CA,C,N\n")
#            fcom.write("cluster out testcluster representative pdb average pdb averagelinkage epsilon 1.0 clusters 100 rms @*\n")
            fcom.write("trajout protein4autodock.pdb pdb\n")
            fcom.write("EOF\n")
            fcom.close()
            os.system("cd %s; chmod a+x RunCommands; ./RunCommands\n" % mult_frames_file)
            os.system("rm protein4autodock.top protein4autodock.trj RunCommands\n")
            
            # cluster trajectory
            num_frames = len(os.listdir(mult_frames_file))
            d_traj_clust = clusterTrajDialog(app.root)
            qt_traj_clustering(mult_frames_file, d_traj_clust.rmsd_max, d_traj_clust.min_in_cluster, d_traj_clust.aim_cluster_num, nsize.cb_x - 0.5*nsize.bl_x, nsize.cb_x + 0.5*nsize.bl_x, nsize.cb_y - 0.5*nsize.bl_y, nsize.cb_y + 0.5*nsize.bl_y, nsize.cb_z - 0.5*nsize.bl_z, nsize.cb_z + 0.5*nsize.bl_z)
            os.system("rm protein4autodock.pdb.*\n")
            
            fcom = open("RunCommands", "w")
#            fcom.write("cd %s\n" % mult_frames_file)
            j = os.listdir(mult_frames_file)
            j.sort()
            ci = 0
            for i in j:
                try:
                    tmpc, tmpc, a = i.split(".", 2)
                    ai = int(a)
                    ktmpa = ai % 100000
                    k100000 = (ai - ktmpa)/100000
                    ktmpb = ktmpa % 10000
                    k10000 = (ktmpa - ktmpb)/10000
                    ktmpa = ktmpb % 1000
                    k1000 = (ktmpb - ktmpa)/1000
                    ktmpb = ktmpa % 100
                    k100 = (ktmpa - ktmpb)/100
                    ktmpa = ktmpb % 10
                    k10 = (ktmpb - ktmpa)/10
                    ktmpb = ktmpa % 1
                    k1 = (ktmpa - ktmpb)
                    fcom.write("mv %s tmp.pdb\n" % i)
    #                os.system("mv %s tmp.pdb" % i)
                    if nprep.add_hydrogens == 1:
                        fcom.write("pythonsh ../prepare_receptor4.py -r tmp.pdb -o preparedProtein_%d%d%d%d%d%d.pdbqt -A 'checkhydrogens' -U nphs\n" % (k100000, k10000, k1000, k100, k10, k1))
                    else:
                        fcom.write("pythonsh ../prepare_receptor4.py -r tmp.pdb -o preparedProtein_%d%d%d%d%d%d.pdbqt -A 'None' -U nphs\n" % (k100000, k10000, k1000, k100, k10, k1))
                    fcom.write("rm tmp.pdb\n")
    #                os.system("rm tmp.pdb")
                    ci += 1
                except:
                    pass
            fcom.close()
            
#        os.system("rm ../prepare_receptor4.py")
        if nprep.object_name == "visual" or nprep.object_name == "file_single":
            cmd.load("preparedProtein.pdb")
            # define search volume
            nsize = sizeSearchSpaceDialog(app.root)
            if nsize.ok_flag == 0:
                return 0
                
            
        # change into project directory
        os.chdir(md_outfile)

        fconf = open("vina.config", "w")
        fconf.write("\n")
        fconf.write("center_x =  %f\n" % (nsize.cb_x))
        fconf.write("center_y =  %f\n" % (nsize.cb_y))
        fconf.write("center_z =  %f\n" % (nsize.cb_z))
        fconf.write("\n")
        fconf.write("size_x = %f\n" % (nsize.bl_x))
        fconf.write("size_y = %f\n" % (nsize.bl_y))
        fconf.write("size_z = %f\n" % (nsize.bl_z))
        fconf.write("\n")
        fconf.write("energy_range = %d\n" % nsize.energy_binding_modes)
        fconf.write("num_modes = %d\n" % nsize.num_binding_modes)
        fconf.close()
        
        # write monitor file
        fmon = open("%s/Monitor.aut" % (md_outfile), "w")
        fmon.write("%s\n" % (nprep.project_dir_short))
        if nprep.object_name == "visual" or nprep.object_name == "file_single":
            fmon.write("single 1\n")
            fmon.write("%d\n" % nsize.num_binding_modes)
        else:
            fmon.write("multiple %d\n" % (ci))
            fmon.write("%d\n" % nsize.num_binding_modes)

        
        fcom = open("RunScript", "w")
        fcom.write("mkdir out\n")
        if nprep.object_name == "visual" or nprep.object_name == "file_single":
            if len(nsize.flex_residues) > 1:
                fcom.write("pythonsh ./prepare_flexreceptor4.py -r frames/preparedProtein.pdbqt -s %s -g frames/rigidProtein.pdbqt -x frames/flexProtein.pdbqt\n" % nsize.flex_residues)
            j = os.listdir(lib_infile)
            j.sort()
            for i in j:
                if i.find(".pdbqt") >= 0:
                    a, b = i.split(".", 1)
                    fmon.write("%s\n" % a)
                    fcom.write("cd out; mkdir %s; cd ..\n" % a)
                    if len(nsize.flex_residues) > 1:
                        fcom.write("%s/vina --ligand in/%s --receptor frames/rigidProtein.pdbqt --flex frames/flexProtein.pdbqt --config vina.config --out out/%s/%s --log out/%s/%s.log\n" % (vina_dir, i, a, i, a, a))
                    else:
                        fcom.write("%s/vina --ligand in/%s --receptor frames/preparedProtein.pdbqt --config vina.config --out out/%s/%s --log out/%s/%s.log\n" % (vina_dir, i, a, i, a, a))
        else:
            j = os.listdir(lib_infile)
            j.sort()
            for i in j:
               if i.find(".pdbqt") >= 0:
                    a, b = i.split(".", 1)
                    fmon.write("%s\n" % a)
                    fcom.write("cd out; mkdir %s; cd ..\n" % a)
                    for ai in range(1, ci+1):
                        ktmpa = ai % 100000
                        k100000 = (ai - ktmpa)/100000
                        ktmpb = ktmpa % 10000
                        k10000 = (ktmpa - ktmpb)/10000
                        ktmpa = ktmpb % 1000
                        k1000 = (ktmpb - ktmpa)/1000
                        ktmpb = ktmpa % 100
                        k100 = (ktmpa - ktmpb)/100
                        ktmpa = ktmpb % 10
                        k10 = (ktmpb - ktmpa)/10
                        ktmpb = ktmpa % 1
                        k1 = (ktmpa - ktmpb)
                        fcom.write("%s/vina --ligand in/%s --receptor frames/preparedProtein_%d%d%d%d%d%d.pdbqt --config vina.config --out out/%s/%s_%d%d%d%d%d%d.pdbqt --log out/%s/%s_%d%d%d%d%d%d.log\n" % (vina_dir, i, k100000, k10000, k1000, k100, k10, k1, a, a, k100000, k10000, k1000, k100, k10, k1, a, a, k100000, k10000, k1000, k100, k10, k1))

        fcom.write("cd %s; touch autodock.done\n" % (md_outfile))
        fcom.close()
            
        fmon.close()

        
        if nprep.object_name == "visual" or nprep.object_name == "file_single":
            childp = subprocess.Popen("cd %s; chmod a+x RunScript; nohup ./RunScript" % (md_outfile), shell=True)
        else:
            childp = subprocess.Popen("cd %s; chmod a+x RunCommands; nohup ./RunCommands; cd %s; chmod a+x RunScript; nohup ./RunScript" % (mult_frames_file, md_outfile), shell=True)
	tkMessageBox.showinfo('Vina calculation has been started',"Data can be read in by using 'Vina, Import Results'.")

        os.chdir(curdir)
        

#######################################################################################
def pdbqt2pdb(infile, outfile):
    fi = open(infile, "r")
    fo = open(outfile, "w")
    for i in fi:
        if i.find("ATOM") >= 0 or i.find("HETAT<") >= 0:
            fo.write("%s\n" % i[:67])
    fi.close()
    fo.close()
    

#######################################################################################
def rmsd_calc(j, k, numatoms):
    rmsd = 0.0;
    for ci in range(numatoms):
        dx = coord[j][ci][0] - coord[k][ci][0]
        dy = coord[j][ci][1] - coord[k][ci][1]
        dz = coord[j][ci][2] - coord[k][ci][2]
        rmsd += dx*dx + dy*dy + dz*dz
    rmsd /= numatoms
    rmsd = math.sqrt(rmsd)
    
    return rmsd
        
              
#######################################################################################
def qt_traj_clustering(indir, rmsd_max, min_in_cluster, aim_cluster_num, min_x, max_x, min_y, max_y, min_z, max_z):
    # allocate memory
    global coord 
    global num_clusters
    global in_cluster
    global num_in_cluster
    global cluster_number
    global cluster_mean_affinity
    global cluster_min_affinity
    global cluster_max_affinity
    global cluster_bin_affinity
    global cluster_stdv_affinity
    
    coord = [[] for k in range(num_frames)]
    flag_in_cluster = [0 for k in range(num_frames)]
    in_cluster = [[] for k in range(num_frames)]
    num_in_cluster = [0 for k in range(num_frames)]
    cluster_number = [0 for k in range(num_frames)]
    cluster_file = [[] for k in range(num_frames)]
    tmp_cluster = [0 for k in range(num_frames)]
    tmp_cluster2 = [0 for k in range(num_frames)]
    rmsd_value = [[] for k in range(num_frames)]
    for j in range(num_frames):
        rmsd_value[j] = [0.0 for k in range(num_frames)]
    # read pdbqt files and store coordinates
    lf = os.listdir(indir)
    lf.sort()
    j = 0
    for nn in lf:
        cluster_file[j] = nn
        in_cluster[j] = [0 for k in range(num_frames)]
        # find number of atoms
        if j == 0:
            resid_list = []
            resnam_list = []
            old_resid = -99999
            old_resnam = ""
            flag = 0
        
            try:
                fi = open(nn, "r")
            except:
                print "Cannot open file %s/\n" % nn
                return 0
            ci = 0
            lines = fi.readlines()
            # residue list
            for line in lines:
                if line.find("ATOM") >= 0:
                    if line[13:14] != "H" or (line[17:20] == "SER" and line[13:15] == "HG"):
                        resid = int(line[21:26].strip())
                        resnam = line[17:20]
                        if resid != old_resid or resnam != old_resnam:
                            if flag == 1:
                                resid_list.append(resid)
                                resnam_list.append(resnam)
                            flag = 0
                            old_resid = resid
                            old_resnam = resnam
                        c_x = float(line[30:38])  
                        c_y = float(line[38:46])
                        c_z = float(line[46:54])
                        if c_x < max_x and c_x > min_x and c_y < max_y and c_y > min_y and c_z < max_z and c_z > min_z:
                            flag = 1
            # last residue
            if flag == 1:
                resid_list.append(resid)
                resnam_list.append(resnam)
            
            # number of atoms in residue list            
            for line in lines:
                if line.find("ATOM") >= 0:
                    if line[13:14] != "H" or (line[17:20] == "SER" and line[13:15] == "HG"):
                        resid = int(line[21:26].strip())
                        resnam = line[17:20]
                        flag = 0
                        ri = 0
                        for rr in resid_list:
                            rn = resnam_list[ri]
                            ri += 1
                            if resid == rr and resnam ==  rn:
                                ci += 1
            numatoms = ci
            fi.close()
        # allocate memory
        coord[j] = [[] for k in range(numatoms)]
        # read coordinates
        try:
            fi = open(nn, "r")
        except:
            print "Cannot open file %s/\n" % nn
            return 0
        ci = 0
        lines = fi.readlines()
        for line in lines:
            if line.find("ATOM") >= 0:
                if line[13:14] != "H" or (line[17:20] == "SER" and line[13:15] == "HG"):
                    resid = int(line[21:26].strip())
                    resnam = line[17:20]
                    ri = 0
                    for rr in resid_list:
                        rn = resnam_list[ri]
                        ri += 1
                        if resid == rr and resnam ==  rn:
                            coord[j][ci] = [0.0 for k in range(3)]
                            c_x = float(line[30:38])  
                            c_y = float(line[38:46])
                            c_z = float(line[46:54])
                            coord[j][ci][0] = c_x
                            coord[j][ci][1] = c_y
                            coord[j][ci][2] = c_z
                            ci += 1
        fi.close()
        j += 1
        
    # pre-calculate rmsd values between frames
    for j in range(num_frames):
        rmsd_value[j][j] = 0.0
        for k in range(j+1, num_frames):
            rmsd_value[j][k] = rmsd_calc(j, k, numatoms)
            rmsd_value[k][j] = rmsd_value[j][k]
            


    # find largest clusters step-by-step
    num_clusters = 0
    while num_clusters > 1.5*aim_cluster_num or num_clusters < 0.5*aim_cluster_num:
        for k in range(num_frames):
            flag_in_cluster[k] = 0
            tmp_cluster[k] = 0
            tmp_cluster2[k] = 0
            for m in range(num_frames):
                in_cluster[k][m] = 0
        cluster_found = 1
        cc = 0
        while cluster_found == 1:
            # find next largest cluster
            num_tmp_cluster2 = 0
            for j in range(num_frames):
                # if not already assigned to any cluster
                if flag_in_cluster[j] == 0:
                    # determine all conformations with RMSD smaller than given value
                    num_tmp_cluster = 0
                    for k in range(num_frames):
                        # if not already assigned to any cluster
                        if flag_in_cluster[k] == 0:
                            if rmsd_value[j][k] < rmsd_max:
                                tmp_cluster[num_tmp_cluster] = k
                                num_tmp_cluster += 1
                    # save tmp_cluster if currently largest
                    if num_tmp_cluster > num_tmp_cluster2:
                        for k in range(num_tmp_cluster):
                            tmp_cluster2[k] = tmp_cluster[k]
                        num_tmp_cluster2 = num_tmp_cluster
                        tmp_cluster2_number = j
            # if actually largest cluster is larger than threshold than store information and go to next cluster
            if num_tmp_cluster2 >= min_in_cluster:
                for k in range(num_tmp_cluster2):
                    in_cluster[cc][k] = tmp_cluster2[k]
                    flag_in_cluster[tmp_cluster2[k]] = 1
                num_in_cluster[cc] = num_tmp_cluster2
                cluster_number[cc] = tmp_cluster2_number            # center of cluster
#                os.system("cp %s/%s %s/cluster.pdb.%d" % (indir, cluster_file[cluster_number[cc]], indir, cc+1))
#                nn = "cluster.pdb.%d" % (cc+1)
#                cmd.load(nn)
                cc += 1                                             # current cluster number
            else:
                cluster_found = 0
        num_clusters = cc
        print "num_clusters: %d (%.2f)" % (num_clusters, rmsd_max)
        if num_clusters > 1.5*aim_cluster_num:
            rmsd_max += 0.01
        elif num_clusters < 0.5*aim_cluster_num:
            rmsd_max -= 0.01

    for k in range(num_clusters):
        os.system("cp %s/%s %s/cluster.pdb.%d" % (indir, cluster_file[cluster_number[k]], indir, k+1))
        nn = "cluster.pdb.%d" % (k+1)
        cmd.load(nn)
#        print "cluster found: %d : %d %d" % (k, cluster_number[k], num_in_cluster[k])
#        print in_cluster[k]
        
            
                

                    
                
        
    
#######################################################################################
def qt_clustering(indir, rmsd_max, min_in_cluster):
    # allocate memory
    global coord 
    global num_clusters
    global in_cluster
    global num_in_cluster
    global cluster_number
    global cluster_mean_affinity
    global cluster_min_affinity
    global cluster_max_affinity
    global cluster_bin_affinity
    global cluster_stdv_affinity
    
    num_clusters = [0 for k in range(num_ligands)]
    in_cluster = [[] for k in range(num_ligands)]
    num_in_cluster = [[] for k in range(num_ligands)]
    cluster_number = [[] for k in range(num_ligands)]
    cluster_mean_affinity = [[] for k in range(num_ligands)]
    cluster_min_affinity = [[] for k in range(num_ligands)]
    cluster_max_affinity = [[] for k in range(num_ligands)]
    cluster_bin_affinity = [[] for k in range(num_ligands)]
    cluster_stdv_affinity = [[] for k in range(num_ligands)]
    
    for lig_i in range(num_ligands):
        if num_solutions[lig_i] > 0:
            coord = [[] for k in range(num_solutions[lig_i])]
            flag_in_cluster = [0 for k in range(num_solutions[lig_i])]
            in_cluster[lig_i] = [[] for k in range(num_solutions[lig_i])]
            num_in_cluster[lig_i] = [0 for k in range(num_solutions[lig_i])]
            cluster_number[lig_i] = [0 for k in range(num_solutions[lig_i])]
            tmp_cluster = [0 for k in range(num_solutions[lig_i])]
            tmp_cluster2 = [0 for k in range(num_solutions[lig_i])]
            # read pdbqt files and store coordinates
            for j in range(num_solutions[lig_i]):
                in_cluster[lig_i][j] = [0 for k in range(num_solutions[lig_i])]
                lig_name = sol_name[lig_i][j]
                os.chdir("%s/out/%s" % (indir, g_ligand_name[lig_i]))
                nn = "%s/out/%s/%s.pdbqt" % (indir, g_ligand_name[lig_i], lig_name)
                # find number of atoms
                if j == 0:
                    fi = open(nn, "r")
                    ci = 0
                    while fi:
                        line = fi.readline()
                        if line.find("ENDMDL") >= 0:
                            break
                        if line.find("ATOM") >= 0:
                            ci += 1
                    numatoms = ci
                    fi.close()
                # allocate memory
                coord[j] = [[] for k in range(numatoms)]
                # read coordinates
                fi = open(nn, "r")
                ci = 0
                while fi:
                    line = fi.readline()
                    if line.find("ENDMDL") >= 0:
                        break
                    if line.find("ATOM") >= 0:
                        coord[j][ci] = [0.0 for k in range(3)]
                        c_x = float(line[30:38])  
                        c_y = float(line[38:46])
                        c_z = float(line[46:54])
                        coord[j][ci][0] = c_x
                        coord[j][ci][1] = c_y
                        coord[j][ci][2] = c_z
                        ci += 1
                fi.close()


            # find largest clusters step-by-step
            cluster_found = 1
            cc = 0
            while cluster_found == 1:
                # find next largest cluster
                num_tmp_cluster2 = 0
                for j in range(num_solutions[lig_i]):
                    # if not already assigned to any cluster
                    if flag_in_cluster[j] == 0:
                        # determine all conformations with RMSD smaller than given value
                        num_tmp_cluster = 0
                        for k in range(num_solutions[lig_i]):
                            # if not already assigned to any cluster
                            if flag_in_cluster[j] == 0:
                                if rmsd_calc(j, k, numatoms) < rmsd_max:
                                    tmp_cluster[num_tmp_cluster] = k
                                    num_tmp_cluster += 1
                        # save tmp_cluster if currently largest
                        if num_tmp_cluster > num_tmp_cluster2:
                            for k in range(num_tmp_cluster):
                                tmp_cluster2[k] = tmp_cluster[k]
                            num_tmp_cluster2 = num_tmp_cluster
                            tmp_cluster2_number = j
                # if actually largest cluster is larger than threshold than store information and go to next cluster
                if num_tmp_cluster2 >= min_in_cluster:
                    for k in range(num_tmp_cluster2):
                        in_cluster[lig_i][cc][k] = tmp_cluster2[k]
                        flag_in_cluster[tmp_cluster2[k]] = 1
                    num_in_cluster[lig_i][cc] = num_tmp_cluster2
                    cluster_number[lig_i][cc] = tmp_cluster2_number            # center of cluster
                    nn = "%s/out/%s/%s.pdbqt" % (indir, g_ligand_name[lig_i], sol_name[lig_i][cluster_number[lig_i][cc]])
                    cmd.load(nn)
                    cc += 1                                             # current cluster number
                    print "cluster found: %d : %d %d" % (cc, cluster_number[lig_i][cc-1], num_in_cluster[lig_i][cc-1])
                else:
                    cluster_found = 0
            num_clusters[lig_i] = cc
                    
            # determine mean value and standard deviation for each cluster
            cluster_mean_affinity[lig_i] = [0.0 for k in range(num_clusters[lig_i])]
            cluster_min_affinity[lig_i] = [0.0 for k in range(num_clusters[lig_i])]
            cluster_max_affinity[lig_i] = [0.0 for k in range(num_clusters[lig_i])]
            cluster_bin_affinity[lig_i] = [[] for k in range(num_clusters[lig_i])]
            cluster_stdv_affinity[lig_i] = [0.0 for k in range(num_clusters[lig_i])]
            for j in range(num_clusters[lig_i]):
                cluster_bin_affinity[lig_i][j] = [0 for k in range(400)]
                for k in range(num_in_cluster[lig_i][j]):
                    tmp3 = sol_affinity[lig_i][in_cluster[lig_i][j][k]] 
                    cluster_mean_affinity[lig_i][j] += tmp3  
                    if tmp3 < cluster_min_affinity[lig_i][j]:
                        cluster_min_affinity[lig_i][j] = tmp3
                    if tmp3 > cluster_max_affinity[lig_i][j]:
                        cluster_max_affinity[lig_i][j] = tmp3                                
                    val1 = int((tmp3 + 95.0)/0.25)
                    cluster_bin_affinity[lig_i][j][val1] += 1
                cluster_mean_affinity[lig_i][j] /= num_in_cluster[lig_i][j]
                # determine standard deviation
                stdv = 0.0
                for k in range(num_in_cluster[lig_i][j]):
                    df = cluster_mean_affinity[lig_i][j] - sol_affinity[lig_i][in_cluster[lig_i][j][k]]
                    stdv += df*df
                stdv /= num_in_cluster[lig_i][j]
                cluster_stdv_affinity[lig_i][j] = math.sqrt(stdv)

            
                

                    
                
        
    
#######################################################################################
def import_results(app):
    global res_num_poses
    global num_ligands
    global mean_affinity
    global stdv_affinity
    global min_affinity
    global max_affinity
    global bin_affinity
    global g_ligand_name
    global order2
    global sol_name
    global sol_affinity
    global sol_vdw
    global sol_coul
    global sol_rf
    global sol_cav
    global sol_const
    global order
    global num_solutions
    global num_frames
    global import_dir
    
    ftypes=(('aut file', '*.aut'), ("All files", "*"))
    curdir = os.getcwd()
    openfile = askopenfilename(initialdir=curdir,
                               filetypes=ftypes)
    if openfile:
        import_dir = os.path.dirname(openfile)
        os.chdir(import_dir)
        
        fin = open(openfile, "r")
        tmp = fin.readline()
        proj_dir = tmp.strip()
        tmp = fin.readline()
        tmp2 = tmp.strip()
        multi_or_single, tmp3 = tmp2.split(None, 1)
        num_frames = int(tmp3)
        tmp = fin.readline()
        res_num_poses= int(tmp)
        
        
        flag = 0
        ci = 0
        while flag == 0:
            tmp = fin.readline()
            tmp2 = tmp.strip()
            if len(tmp2) == 0:
                flag = 1
            else:
                g_ligand_name[ci] = tmp2
                ci += 1
            
        num_ligands = ci
        
        flag = 0
        for ent in os.listdir(import_dir):
            if ent.find('autodock.done') >= 0:
                flag = 1
        if flag == 0:
            tkMessageBox.showwarning(
                "Could not find autodock.done",
                "AutoDock docking might be still in progress or didn't finish successfully"
            )
            return 0
                
        co = 0
        for na in cmd.get_names("objects"):
            co += 1
            
        if co > 0:
            if tkMessageBox.askokcancel("Read project", "Reading results will delete all exisiting objects in current session.") == 0:
                return 0
        # remove existing objects
        for na in cmd.get_names("objects"):
            try:
                cmd.remove(na)
            except:
                pass
            cmd.delete(na)

    
    if multi_or_single == "single":
        sol_name = [[] for i in range(num_ligands)]
        sol_affinity = [[] for i in range(num_ligands)]
        mean_affinity = [0.0 for i in range(num_ligands)]
        num_solutions = [0 for i in range(num_ligands)]
        order2 = [0 for i in range(ci)]
        for j in range(num_ligands):
            sol_name[j] = [[] for k in range(res_num_poses)]
            sol_affinity[j] = [0.0 for k in range(res_num_poses)]
            lig_name = g_ligand_name[j]
            fn = import_dir + "/out/" + lig_name + "/" + lig_name + ".log"
            fi = open(fn, "r")
            fi.seek(0)
            flag = 0
            ci = 0
            while fi:
                line = fi.readline()
                if line.find("Writing output ... done") >= 0:
                    flag = 0
                    break
                if flag == 1:
                    tmp, tmp2, tmp = line.split(None, 2)
                    sol_name[j][ci] = lig_name
                    sol_affinity[j][ci] = float(tmp2)
                    
                    os.chdir("%s/out/%s" % (import_dir, lig_name))
                    if ci == 0:
                        nn = import_dir + "/out/" + lig_name + "/" + lig_name + ".pdbqt"
                        cmd.load(nn)
                        mean_affinity[j] = sol_affinity[j][0]                
                    ci += 1
                if line.find("-----+") >= 0:
                    flag = 1
            num_solutions[j] = ci
            if ci == 0:
                mean_affinity[j] = 999.9
            
        # sort ligand
        for k in range(num_ligands):
            order2[k] = k
        for k in range(num_ligands-1):
            for m in range(k+1, num_ligands):
                if mean_affinity[order2[m]] < mean_affinity[order2[k]]:
                    temp = order2[k]
                    order2[k] = order2[m]
                    order2[m] = temp
    
        d_res = displayResultsDialog(app.root, 1)

        fi.close()
        
    else:
        d_clust = clusterDialog(app.root)
        
        sol_name = [[] for i in range(num_ligands)]
        sol_affinity = [[] for i in range(num_ligands)]
        mean_affinity = [0.0 for i in range(num_ligands)]
        num_solutions = [0 for i in range(num_ligands)]
        bin_affinity = [[] for i in range(num_ligands)]
        min_affinity = [9999.9 for i in range(num_ligands)]
        max_affinity = [-9999.9 for i in range(num_ligands)]
        stdv_affinity = [0.0 for i in range(num_ligands)]
        order2 = [0 for i in range(ci)]
        for j in range(num_ligands):
            sol_name[j] = [[] for k in range(num_frames)]
            sol_affinity[j] = [0.0 for k in range(num_frames)]
            bin_affinity[j] = [0 for k in range(400)] # from -95 to +5 in 0.25 kcal/mol bins
            lig_name = g_ligand_name[j]
            ai = 1
            for k in range(num_frames):
                aj = k+1
                ktmpa = aj % 100000
                k100000 = (aj - ktmpa)/100000
                ktmpb = ktmpa % 10000
                k10000 = (ktmpa - ktmpb)/10000
                ktmpa = ktmpb % 1000
                k1000 = (ktmpb - ktmpa)/1000
                ktmpb = ktmpa % 100
                k100 = (ktmpa - ktmpb)/100
                ktmpa = ktmpb % 10
                k10 = (ktmpb - ktmpa)/10
                ktmpb = ktmpa % 1
                k1 = (ktmpa - ktmpb)
#                fn = import_dir + "/out/" + lig_name + "/" + lig_name + "_" + ".log"
                fn = "%s/out/%s/%s_%d%d%d%d%d%d.log" % (import_dir, lig_name, lig_name, k100000, k10000, k1000, k100, k10, k1)
                fi = open(fn, "r")
                fi.seek(0)
                flag = 0
                while fi:
                    line = fi.readline()
                    if line.find("Writing output ... done") >= 0:
                        flag = 0
                        break
                    if flag == 1:
                        tmp, tmp2, tmp = line.split(None, 2)
                        tmp3 = float(tmp2)
                        if tmp3 < 5.0:
                            sol_affinity[j][ai-1] = tmp3
                            sol_name[j][ai-1] = "%s_%d%d%d%d%d%d" % (lig_name, k100000, k10000, k1000, k100, k10, k1)
                            
#                            os.chdir("%s/out/%s" % (import_dir, lig_name))
#                            nn = "%s/out/%s/%s_%d%d%d%d%d%d.pdbqt" % (import_dir, lig_name, lig_name, k100000, k10000, k1000, k100, k10, k1)
#                            cmd.load(nn)
                            # calculate mean energy
                            mean_affinity[j] += tmp3   
                            if tmp3 < min_affinity[j]:
                                min_affinity[j] = tmp3
                            if tmp3 > max_affinity[j]:
                                max_affinity[j] = tmp3                                
                            val1 = int((tmp3 + 95.0)/0.25)
                            bin_affinity[j][val1] += 1
                            ai += 1
                        break
                    if line.find("-----+") >= 0:
                        flag = 1
            num_solutions[j] = ai-1
            # determine mean value
            if ai > 1:
                mean_affinity[j] /= (ai-1)
                # determine standard deviation
                stdv = 0.0
                for k in range(num_solutions[j]):
                    df = mean_affinity[j] - sol_affinity[j][k]
                    stdv += df*df
                stdv /= num_solutions[j]
                stdv_affinity[j] = math.sqrt(stdv)
            else:
                mean_affinity[j] = 999.9
                stdv_affinity[j] = 0.0
        print "qt_clustering: start"        
        qt_clustering(import_dir, d_clust.rmsd_max, d_clust.min_in_cluster)
        print "qt_clustering: end"        

            
        # sort ligand
        for k in range(num_ligands):
            order2[k] = k
        for k in range(num_ligands-1):
            for m in range(k+1, num_ligands):
                if mean_affinity[order2[m]] < mean_affinity[order2[k]]:
                    temp = order2[k]
                    order2[k] = order2[m]
                    order2[m] = temp
#        for k in range(num_ligands):
#            print "%d %s %f" % (k, sol_name[order2[k]][0], mean_affinity[order2[k]])
    
        d_res = displayResultsDialog(app.root, 2)

        fi.close()
                    

#######################################################################################
def pdbqt_to_pdb(MAX_COUNTER, results_folder):

    TARGET_PATH = import_dir
    SYMPOSAR_PATH = import_dir + "/SYMPOSAR"
    SYMPOSAR_RESULTS_PATH = SYMPOSAR_PATH + "/SYMPOSAR_RESULTS"
    PDB_DIR = import_dir + "/in/PDBs"
    QT_DIR = import_dir + "/out"
#    MAX_COUNTER = 10

    if os.path.isdir(PDB_DIR) == 0 or os.path.isdir(QT_DIR) == 0:
        tkMessageBox.showwarning(
                "Could not find directories",
                "Could not find directories %s and %s." % (PDB_DIR, QT_DIR)
        )
        return 0

    pdb_list=[]
    pdb_atom=list()
    pdb_connect=list()
  
 
    PDB_PATH=os.path.join(TARGET_PATH, "Solutions_PDB")
    MOL2_PATH=os.path.join(SYMPOSAR_RESULTS_PATH, results_folder)

    if os.path.isdir(PDB_PATH):
        commandstring = "rm -r "+PDB_PATH
        os.system(commandstring)
    if os.path.isdir(MOL2_PATH):
        commandstring = "rm -r "+MOL2_PATH
        os.system(commandstring) 
    
    if os.path.isdir(SYMPOSAR_PATH) == 0:
        os.mkdir(SYMPOSAR_PATH)
        os.mkdir(SYMPOSAR_RESULTS_PATH)
    os.mkdir(PDB_PATH)
    os.mkdir(MOL2_PATH)

    for root,dirs,files in os.walk(PDB_DIR):
        name_count = 0
            
        for name in files:
        
            if name[-4:] == ".pdb":

                pdb_list.append(name)
                pdb_file = open(os.path.join(root,name), 'r')
                pdb_file2 = open(os.path.join(root,name), 'r')
                filelines = pdb_file.readlines()
                filelines2 = pdb_file2.readlines()
                pdb_atom.append(filelines)
                pdb_connect.append(filelines2)
                for i in range(len(pdb_atom[name_count])-1,-1,-1):
                    string = pdb_atom[name_count][i]
                    pdb_atom[name_count].pop(i)
#                    if cmp(string[0:6],"HETATM") != 0:
#                        pdb_atom[name_count].pop(i)
                  
                for i in range(len(pdb_connect[name_count])-1,-1,-1):
                    string = pdb_connect[name_count][i]
                    if cmp(string[0:6],"CONECT") != 0:
                        pdb_connect[name_count].pop(i)
                         
                name_count += 1
         
                pdb_file.close()
                pdb_file2.close()



   
    for root,dirs,files in os.walk(QT_DIR):              
        for name in files:
            name_incrementer = 97
            connect_conversion_table=list()
            counter = 0
            converted_list = list()
            if cmp(name[-6:],".pdbqt") == 0:
                for i in range(len(pdb_list)):
                    str1 = pdb_list[i]
                    if cmp(str1[:-4],name[:-6])==0:
                        qt_file = open(os.path.join(root,name), 'r')
                        for line in qt_file.readlines():
                            if (cmp(line[0:4],"ROOT") == 0 and int(counter)<int(MAX_COUNTER)):
                                file_string ="/"+str1[:-4] + chr(name_incrementer)
                                out_file = open(PDB_PATH+file_string+".pdb", 'w')
                            if (cmp(line[0:4],"ATOM") == 0 and int(counter)<int(MAX_COUNTER)):
                                out_string = "HETATM  "+line[8:21]+"    1    "+line[30:68]+"                   \n"
                                out_string = out_string[0:77]+line[13]+out_string[77:]
                                out_file.write(out_string)
                                connect_conversion_table.append(line[12:16])
   
                            if (cmp(line[0:6],"ENDMDL") == 0 and int(counter)<int(MAX_COUNTER)):
                                for j in range(len(pdb_connect[i])):
                                    connect_check = 1
                                    connection_str = pdb_connect[i][j].split()
                                    connection=list()
                                    converted_list = list()
                                    for k in range(len(connection_str)):
                                        if k > 0:
                                                connection.append(connection_str[k])
                                    for l in range(len(connection)):
                                        string = pdb_atom[i][int(connection[l])-1]
                                        if string[12:16] in connect_conversion_table:
                                            converted_list.append(connect_conversion_table.index(string[12:16]))                  
                                        else:
                                            connect_check = 0
                                 
                                    if connect_check == 1:
                                        outline = "CONECT"
                                        for n in range(len(connection)):
                                            string = str(converted_list[n]+1)
                                            if converted_list[n]+1 < 10:
                                                outline += "    "+string[-2:]
                                            elif converted_list[n]+1 < 100:
                                                outline += "   "+string[-2:]
                                            else:
                                                outline += "  "+string[-2:]
                                        outline += "\n"
                                        out_file.write(outline)
                                counter += 1
                                out_file.write("END")
                                out_file.close()
                                name_incrementer += 1
                                #print counter
                                commandstring = "babel -h -i pdb " +PDB_PATH+file_string+".pdb"+" -o mol2 "+MOL2_PATH+file_string+".mol2"
                                os.system(commandstring)
                        qt_file.close()

                                

#######################################################################################
def modifyPDB2standard(pdbfile):
    fi = open(pdbfile, "r")
    fo = open("tmp.tmp", "w")
    for i in fi:
        if i.find("ATOM") >= 0 or i.find("HETATM") >= 0:
            j = i[17:20]
            if j.find("HIE") >= 0:
                fo.write("%sHIS%s" % (i[:17], i[20:]))
            elif j.find("HID") >= 0:
                fo.write("%sHIS%s" % (i[:17], i[20:]))
            elif j.find("HIP") >= 0:
                fo.write("%sHIS%s" % (i[:17], i[20:]))
#            elif j.find("CYX") >= 0:
#                fo.write("%sCYS%s" % (i[:17], i[20:]))
            else:
                fo.write(i)
        else:
            fo.write(i)
    fi.close()
    fo.close()
    os.rename("tmp.tmp", pdbfile)

#######################################################################################
def results_to_raptor(app):
    global g_dg
    
    dlgnum = numPosesDialog(app.root)
    # results have to be read in first
    pdbqt_to_pdb(dlgnum.numPoses4Symposar, dlgnum.results_folder)
    SYMPOSAR_PATH = import_dir + "/SYMPOSAR"
    SYMPOSAR_RESULTS_PATH = SYMPOSAR_PATH + "/SYMPOSAR_RESULTS"
    MOL2_PATH = os.path.join(SYMPOSAR_RESULTS_PATH, dlgnum.results_folder)
    os.chdir(MOL2_PATH)
    # add affinity data
    g_dg = [0.0 for i in range(num_ligands)]
    g_solv = [0.0 for i in range(num_ligands)]
    g_tds = [0.0 for i in range(num_ligands)]
    g_int = [0.0 for i in range(num_ligands)]
    dlgaff = affinityDialog(app.root)
    # rename hydrogen atoms
    for af in os.listdir(MOL2_PATH):
        for i in range(num_ligands):
            if af.find(g_ligand_name[i]) >= 0:
                name = i
                break
        f_new = open("%s/tmp.mol2" % MOL2_PATH, "w")
        f_new.write("@<TRIPOS>RAPTOR_DATA\n")
        tm = float(g_dg[name])
        if tm > 1e-15:
            tmp = 0.58216*math.log(tm) - 12.064256
        else:
            tmp = 0.0
        f_new.write("G    %f\n" % (tmp))
        
        tm = float(g_solv[name])
        if tm > 1e-15:
            tmp = 0.58216*math.log(tm) - 12.064256
        else:
            tmp = 0.0
        f_new.write("Solv %f\n" % (tmp))
        
        tm = float(g_tds[name])
        if tm > 1e-15:
            tmp = 0.58216*math.log(tm) - 12.064256
        else:
            tmp = 0.0
        f_new.write("TdS  %f\n" % (tmp))
        
        tm = float(g_int[name])
        if tm > 1e-15:
            tmp = 0.58216*math.log(tm) - 12.064256
        else:
            tmp = 0.0
        f_new.write("Int  %f\n" % (tmp))
        
        flag = 0
        countH = 0
        fname = "%s/%s" % (MOL2_PATH, af)
        fi = open(fname, "r")
        for dd in fi:
            if dd.find("@<TRIPOS>BOND") >= 0:
                flag = 0
            if flag == 1:
                a, b, c = dd.split(None, 2)
                b2 = b.strip()
                if b2[0] == "H":
                    countH += 1
                    f_new.write("%8s H%-2d%s" % (dd[:8], countH, dd[12:]))
                else:
                    f_new.write(dd)
            else:
                f_new.write(dd)
            if dd.find("@<TRIPOS>ATOM") >= 0:
                flag = 1
        f_new.close()
        os.system("mv %s/tmp.mol2 %s/%s" % (MOL2_PATH, MOL2_PATH, af))
        
    # define training and test set
#    dlgtrte = testtrainDialog(app.root)  

    # write additional Symposar files
    os.chdir(SYMPOSAR_PATH)
    os.system("touch Symposar.done")
    
    fout = open("ligand_data.cmp", "w")
    ci = 0
    j = os.listdir(MOL2_PATH)
    j.sort()
    for af in j:
        if af.find(".mol2") >= 0:
            cf, tmp = af.split(".")
            for i in range(num_ligands):
                if cf.find(g_ligand_name[i]) >= 0:
                    name = i
                    break
            v_train = 'n'
            fout.write('%5d %4s %14.7f %14.7f %14.7f %14.7f %3s\n' % (ci, cf, float (g_dg[name]), float (g_solv[name]), float (g_tds[name]), float (g_int[name]), v_train))
            ci += 1
    fout.close()
                    
    fout = open("Monitor.sym", "w")
    fout.write("%s\n" % import_dir)
    fout.write("SYMPOSAR_RESULTS\n")
    fout.close()

    os.chdir(SYMPOSAR_RESULTS_PATH)
    os.system("touch Symposar.done")

#######################################################################################
def makeBox(cm_x, cm_y, cm_z, dx, dy, dz):
        v = [[]]*8
 
        #make a cuboid with the lower corner on the origin
        v[0] = [cm_x - dx/2, cm_y - dy/2, cm_z - dz/2]                # [0, 0, 0]
        v[1] = [cm_x + dx/2, cm_y - dy/2, cm_z - dz/2]                # [1, 0, 0]
        v[2] = [cm_x + dx/2, cm_y + dy/2, cm_z - dz/2]                # [1, 1, 0]
        v[3] = [cm_x + dx/2, cm_y - dy/2, cm_z + dz/2]                # [1, 0, 1]
        v[4] = [cm_x - dx/2, cm_y + dy/2, cm_z - dz/2]                # [0, 1, 0]
        v[5] = [cm_x - dx/2, cm_y - dy/2, cm_z + dz/2]                # [0, 0, 1]
        v[6] = [cm_x - dx/2, cm_y + dy/2, cm_z + dz/2]                # [0, 1, 1]
        v[7] = [cm_x + dx/2, cm_y + dy/2, cm_z + dz/2]                # [1, 1, 1]
 
        bot =  [v[0], v[1], v[2], v[4]] # O, x, xy, y
        top =  [v[7], v[3], v[5], v[6]] # xyz, xz, z, yz
        minL = [v[0], v[4], v[6], v[5]] # O, y, yz, z
        minR = [v[0], v[5], v[3], v[1]] # O, z, xz, x
        maxL = [v[4], v[2], v[7], v[6]] # y, xy, xyz, yz
        maxR = [v[3], v[1], v[2], v[7]] # xz, x, xy, xyz
        box = [bot, minR, minL, maxR, maxL, top]
 
        return box
 
#######################################################################################
def drawBox(box):
        from pymol.cgo import *
        
        for na in cmd.get_names("objects"):
            if na.find("Box") >= 0:
                cmd.delete(na)
        colourRGB=[1.0,0.85,0.0]
        boxObj = []
        for side in box:
                boxObj.append(BEGIN)
                boxObj.append(LINE_STRIP)
                boxObj.append(COLOR)
                boxObj.extend(colourRGB)
                for point in side:
                        boxObj.append(VERTEX)
                        boxObj.extend(point)
                boxObj.append(END)
 
        cmd.set('auto_zoom', 0)
        cmd.load_cgo(boxObj, "Box")
#        print boxObj
        cmd.set('auto_zoom', 1)

#######################################################################################
class modifySettingsFile:
	def __init__(self, top):
		self.dialog = Pmw.Dialog(top,
								buttons = ('Exit','Save to file','Reset to defaults'),
								defaultbutton = 'Exit',
								title = 'Generate/Modify Settings_Linux.txt file',
								command = self.apply)

		self.curdir = os.getcwd()
 	                                
		master = self.dialog.interior()

		Tkinter.Label(master, text='Generate new and modify existing Settings_Linux.txt file').pack(expand = 1, fill = 'both', padx = 5, pady = 5)

		nb1 = Pmw.NoteBook(master)

		# username page
		user_page = nb1.add('Username')

		Tkinter.Label(user_page, text='username:').grid(row=0, column=0, sticky=W, padx = 5, pady=1)
		self.e_username = Tkinter.Entry(user_page, width=70)
		self.e_username.insert(0, username)
		self.e_username.grid(row=0, column=1, sticky=W, padx = 5, pady=1)


		# MM page
		MM_page = nb1.add('Amber')

		Tkinter.Label(MM_page, text='library_dir (library containing cofactors for Amber):').grid(row=0, column=0, sticky=W, padx = 5, pady=1)
		self.e_library_dir = Tkinter.Entry(MM_page, width=70)
		self.e_library_dir.insert(0, library_dir)
		self.e_library_dir.grid(row=0, column=1, sticky=W, padx = 5, pady=1)
		b_library_dir = Tkinter.Button(MM_page, text='Browse', command = lambda: self.searchFOLDER(library_dir,self.e_library_dir)).grid(row=0, column=2, sticky=W, padx = 5, pady=1)

		Tkinter.Label(MM_page, text='amber_dir (main amber directory = $AMBERHOME):').grid(row=1, column=0, sticky=W, padx = 5, pady=1)
		self.e_amber_dir = Tkinter.Entry(MM_page, width=70)
		self.e_amber_dir.insert(0, amber_dir)
		self.e_amber_dir.grid(row=1, column=1, sticky=W, padx = 5, pady=1)
		b_amber_dir = Tkinter.Button(MM_page, text='Browse', command = lambda: self.searchFOLDER(amber_dir,self.e_amber_dir)).grid(row=1, column=2, sticky=W, padx = 5, pady=1)

		Tkinter.Label(MM_page, text='sie_dir (directory containing SIE executable):').grid(row=2, column=0, sticky=W, padx = 5, pady=1)
		self.e_sie_dir = Tkinter.Entry(MM_page, width=70)
		self.e_sie_dir.insert(0, sie_dir)
		self.e_sie_dir.grid(row=2, column=1, sticky=W, padx = 5, pady=1)
		b_sie_dir = Tkinter.Button(MM_page, text='Browse', command = lambda: self.searchFOLDER(sie_dir,self.e_sie_dir)).grid(row=2, column=2, sticky=W, padx = 5, pady=1)
		
		Tkinter.Label(MM_page, text='changepdb_dir (directory containing changepdb executable):').grid(row=3, column=0, sticky=W, padx = 5, pady=1)
		self.e_changepdb_dir = Tkinter.Entry(MM_page, width=70)
		self.e_changepdb_dir.insert(0, changepdb_dir)
		self.e_changepdb_dir.grid(row=3, column=1, sticky=W, padx = 5, pady=1)
		b_changepdb_dir = Tkinter.Button(MM_page, text='Browse', command = lambda: self.searchFOLDER(changepdb_dir,self.e_changepdb_dir)).grid(row=3, column=2, sticky=W, padx = 5, pady=1)

		Tkinter.Label(MM_page, text='changecrd_dir (directory containing changecrd executable):').grid(row=4, column=0, sticky=W, padx = 5, pady=1)
		self.e_changecrd_dir = Tkinter.Entry(MM_page, width=70)
		self.e_changecrd_dir.insert(0, changecrd_dir)
		self.e_changecrd_dir.grid(row=4, column=1, sticky=W, padx = 5, pady=1)
		b_changecrd_dir = Tkinter.Button(MM_page, text='Browse', command = lambda: self.searchFOLDER(changecrd_dir,self.e_changecrd_dir)).grid(row=4, column=2, sticky=W, padx = 5, pady=1)
		
		Tkinter.Label(MM_page, text='reduce_exe_dir (directory containing reduce executable):').grid(row=5, column=0, sticky=W, padx = 5, pady=1)
		self.e_reduce_exe_dir = Tkinter.Entry(MM_page, width=70)
		self.e_reduce_exe_dir.insert(0, reduce_exe_dir)
		self.e_reduce_exe_dir.grid(row=5, column=1, sticky=W, padx = 5, pady=1)
		b_reduce_exe_dir = Tkinter.Button(MM_page, text='Browse', command = lambda: self.searchFOLDER(reduce_exe_dir,self.e_reduce_exe_dir)).grid(row=5, column=2, sticky=W, padx = 5, pady=1)


		# Vina page
		Vina_page = nb1.add('AutoDock Vina')

		Tkinter.Label(Vina_page, text='autodock_dir (library for ligand datasets):').grid(row=0, column=0, sticky=W, padx = 5, pady=1)
		self.e_autodock_dir = Tkinter.Entry(Vina_page, width=70)
		self.e_autodock_dir.insert(0, autodock_dir)
		self.e_autodock_dir.grid(row=0, column=1, sticky=W, padx = 5, pady=1)
		b_autodock_dir = Tkinter.Button(Vina_page, text='Browse', command = lambda: self.searchFOLDER(autodock_dir,self.e_autodock_dir)).grid(row=0, column=2, sticky=W, padx = 5, pady=1)

		Tkinter.Label(Vina_page, text='autodock_exe_dir (directory containing AutoDock Tools):').grid(row=1, column=0, sticky=W, padx = 5, pady=1)
		self.e_autodock_exe_dir = Tkinter.Entry(Vina_page, width=70)
		self.e_autodock_exe_dir.insert(0, autodock_exe_dir)
		self.e_autodock_exe_dir.grid(row=1, column=1, sticky=W, padx = 5, pady=1)
		b_autodock_exe_dir = Tkinter.Button(Vina_page, text='Browse', command = lambda: self.searchFOLDER(autodock_exe_dir,self.e_autodock_exe_dir)).grid(row=1, column=2, sticky=W, padx = 5, pady=1)

		Tkinter.Label(Vina_page, text='vina_dir (directory containing AutoDock Vina executable):').grid(row=2, column=0, sticky=W, padx = 5, pady=1)
		self.e_vina_dir = Tkinter.Entry(Vina_page, width=70)
		self.e_vina_dir.insert(0, vina_dir)
		self.e_vina_dir.grid(row=2, column=1, sticky=W, padx = 5, pady=1)
		b_vina_dir = Tkinter.Button(Vina_page, text='Browse', command = lambda: self.searchFOLDER(vina_dir,self.e_vina_dir)).grid(row=2, column=2, sticky=W, padx = 5, pady=1)


		# Slide page
		Slide_page = nb1.add('Slide docking')

		Tkinter.Label(Slide_page, text='slide_dir (library for ligand datasets):').grid(row=0, column=0, sticky=W, padx = 5, pady=1)
		self.e_slide_dir = Tkinter.Entry(Slide_page, width=70)
		self.e_slide_dir.insert(0, slide_dir)
		self.e_slide_dir.grid(row=0, column=1, sticky=W, padx = 5, pady=1)
		b_slide_dir = Tkinter.Button(Slide_page, text='Browse', command = lambda: self.searchFOLDER(slide_dir,self.e_slide_dir)).grid(row=0, column=2, sticky=W, padx = 5, pady=1)

		Tkinter.Label(Slide_page, text='slide_main_dir (base directory for Slide):').grid(row=1, column=0, sticky=W, padx = 5, pady=1)
		self.e_slide_main_dir = Tkinter.Entry(Slide_page, width=70)
		self.e_slide_main_dir.insert(0, slide_main_dir)
		self.e_slide_main_dir.grid(row=1, column=1, sticky=W, padx = 5, pady=1)
		b_slide_main_dir = Tkinter.Button(Slide_page, text='Browse', command = lambda: self.searchFOLDER(slide_main_dir,self.e_slide_main_dir)).grid(row=1, column=2, sticky=W, padx = 5, pady=1)

		Tkinter.Label(Slide_page, text='slide_exe_dir (directory containing Slide executable):').grid(row=2, column=0, sticky=W, padx = 5, pady=1)
		self.e_slide_exe_dir = Tkinter.Entry(Slide_page, width=70)
		self.e_slide_exe_dir.insert(0, slide_exe_dir)
		self.e_slide_exe_dir.grid(row=2, column=1, sticky=W, padx = 5, pady=1)
		b_slide_exe_dir = Tkinter.Button(Slide_page, text='Browse', command = lambda: self.searchFOLDER(slide_exe_dir,self.e_slide_exe_dir)).grid(row=2, column=2, sticky=W, padx = 5, pady=1)

		Tkinter.Label(Slide_page, text='slide_main_dir_multiple (base directory for Slide (multiple poses)):').grid(row=3, column=0, sticky=W, padx = 5, pady=1)
		self.e_slide_main_dir_multiple = Tkinter.Entry(Slide_page, width=70)
		self.e_slide_main_dir_multiple.insert(0, slide_main_dir_multiple)
		self.e_slide_main_dir_multiple.grid(row=3, column=1, sticky=W, padx = 5, pady=1)
		b_slide_main_dir_multiple = Tkinter.Button(Slide_page, text='Browse', command = lambda: self.searchFOLDER(slide_main_dir_multiple,self.e_slide_main_dir_multiple)).grid(row=3, column=2, sticky=W, padx = 5, pady=1)

		Tkinter.Label(Slide_page, text='slide_exe_dir_multiple (directory containing Slide executable (multiple poses)):').grid(row=4, column=0, sticky=W, padx = 5, pady=1)
		self.e_slide_exe_dir_multiple = Tkinter.Entry(Slide_page, width=70)
		self.e_slide_exe_dir_multiple.insert(0, slide_exe_dir_multiple)
		self.e_slide_exe_dir_multiple.grid(row=4, column=1, sticky=W, padx = 5, pady=1)
		b_slide_exe_dir_multiple = Tkinter.Button(Slide_page, text='Browse', command = lambda: self.searchFOLDER(slide_exe_dir_multiple,self.e_slide_exe_dir_multiple)).grid(row=4, column=2, sticky=W, padx = 5, pady=1)

		Tkinter.Label(Slide_page, text='slide_data_dir (directory for storing data):').grid(row=5, column=0, sticky=W, padx = 5, pady=1)
		self.e_slide_data_dir = Tkinter.Entry(Slide_page, width=70)
		self.e_slide_data_dir.insert(0, slide_data_dir)
		self.e_slide_data_dir.grid(row=5, column=1, sticky=W, padx = 5, pady=1)
		b_slide_data_dir = Tkinter.Button(Slide_page, text='Browse', command = lambda: self.searchFOLDER(slide_data_dir,self.e_slide_data_dir)).grid(row=5, column=2, sticky=W, padx = 5, pady=1)

		Tkinter.Label(Slide_page, text='slide_sort (directory containing sort program for Slide output):').grid(row=6, column=0, sticky=W, padx = 5, pady=1)
		self.e_slide_sort = Tkinter.Entry(Slide_page, width=70)
		self.e_slide_sort.insert(0, slide_sort)
		self.e_slide_sort.grid(row=6, column=1, sticky=W, padx = 5, pady=1)
		b_slide_sort = Tkinter.Button(Slide_page, text='Browse', command = lambda: self.searchFOLDER(slide_sort,self.e_slide_sort)).grid(row=6, column=2, sticky=W, padx = 5, pady=1)


		# QSAR page
		QSAR_page = nb1.add('Symposar/Raptor')

		Tkinter.Label(QSAR_page, text='symposar_exe_dir (library for ligand datasets):').grid(row=0, column=0, sticky=W, padx = 5, pady=1)
		self.e_symposar_exe_dir = Tkinter.Entry(QSAR_page, width=70)
		self.e_symposar_exe_dir.insert(0, symposar_exe_dir)
		self.e_symposar_exe_dir.grid(row=0, column=1, sticky=W, padx = 5, pady=1)
		b_symposar_exe_dir = Tkinter.Button(QSAR_page, text='Browse', command = lambda: self.searchFOLDER(symposar_exe_dir,self.e_symposar_exe_dir)).grid(row=0, column=2, sticky=W, padx = 5, pady=1)

		Tkinter.Label(QSAR_page, text='raptor_exe_dir (library for ligand datasets):').grid(row=1, column=0, sticky=W, padx = 5, pady=1)
		self.e_raptor_exe_dir = Tkinter.Entry(QSAR_page, width=70)
		self.e_raptor_exe_dir.insert(0, raptor_exe_dir)
		self.e_raptor_exe_dir.grid(row=1, column=1, sticky=W, padx = 5, pady=1)
		b_raptor_exe_dir = Tkinter.Button(QSAR_page, text='Browse', command = lambda: self.searchFOLDER(raptor_exe_dir,self.e_raptor_exe_dir)).grid(row=1, column=2, sticky=W, padx = 5, pady=1)


		# Infrastructure page
		Cluster_page = nb1.add('Cluster infrastructure')
	
		Tkinter.Label(Cluster_page, text='Cluster #:').grid(row=0, column=0, sticky=W, padx = 5, pady=1)
		Tkinter.Label(Cluster_page, text='Cluster name/address:').grid(row=0, column=1, sticky=W, padx = 5, pady=1)
		Tkinter.Label(Cluster_page, text='max. # \nof \nprocessors:').grid(row=0, column=2, sticky=W, padx = 5, pady=1)
		Tkinter.Label(Cluster_page, text='SSH \nport:').grid(row=0, column=3, sticky=W, padx = 5, pady=1)
		Tkinter.Label(Cluster_page, text='Queue:').grid(row=0, column=4, sticky=W, padx = 5, pady=1)
		Tkinter.Label(Cluster_page, text='Home directory on cluster:').grid(row=0, column=5, sticky=W, padx = 5, pady=1)
		Tkinter.Label(Cluster_page, text='$AMBERHOME on cluster:').grid(row=0, column=6, sticky=W, padx = 5, pady=1)


		self.cb_aa = [[] for i in range(10)]
		self.e_SERVER_nodename = ["" for i in range(10)]
		self.e_SERVER_max_proc = ["0" for i in range(10)]
		self.e_SERVER_ssh_port = ["" for i in range(10)]
		self.e_SERVER_home_dir = ["" for i in range(10)]
		self.e_SERVER_amber_dir = ["" for i in range(10)]
		self.e_SERVER_queue_name = ["" for i in range(10)]
		self.v_aa = []
		cj = 0
		for i in range(10):
			if cj < SERVER_number_of:
				var = IntVar()
				self.cb_aa[cj] = Checkbutton(Cluster_page, text=str(i+1), variable=var)
				self.cb_aa[cj].grid(row=cj+1, column=0, sticky=W, padx=5)
				self.v_aa.append(var)
				self.cb_aa[cj].select()

				self.e_SERVER_nodename[cj] = Tkinter.Entry(Cluster_page, width=25)
				self.e_SERVER_nodename[cj].insert(0, SERVER_nodename[cj])
				self.e_SERVER_nodename[cj].grid(row=cj+1, column=1, sticky=W, padx = 5, pady=1)

				self.e_SERVER_max_proc[cj] = Tkinter.Entry(Cluster_page, width=5)
				self.e_SERVER_max_proc[cj].insert(0, SERVER_max_proc[cj])
				self.e_SERVER_max_proc[cj].grid(row=cj+1, column=2, sticky=W, padx = 5, pady=1)

				self.e_SERVER_ssh_port[cj] = Tkinter.Entry(Cluster_page, width=5)
				self.e_SERVER_ssh_port[cj].insert(0, SERVER_ssh_port[cj])
				self.e_SERVER_ssh_port[cj].grid(row=cj+1, column=3, sticky=W, padx = 5, pady=1)

				self.e_SERVER_queue_name[cj] = Tkinter.Entry(Cluster_page, width=15)
				self.e_SERVER_queue_name[cj].insert(0, SERVER_queue_name[cj])
				self.e_SERVER_queue_name[cj].grid(row=cj+1, column=4, sticky=W, padx = 5, pady=1)

				self.e_SERVER_home_dir[cj] = Tkinter.Entry(Cluster_page, width=25)
				self.e_SERVER_home_dir[cj].insert(0, SERVER_home_dir[cj])
				self.e_SERVER_home_dir[cj].grid(row=cj+1, column=5, sticky=W, padx = 5, pady=1)

				self.e_SERVER_amber_dir[cj] = Tkinter.Entry(Cluster_page, width=30)
				self.e_SERVER_amber_dir[cj].insert(0, SERVER_amber_dir[cj])
				self.e_SERVER_amber_dir[cj].grid(row=cj+1, column=6, sticky=W, padx = 5, pady=1)
			else:
				var = IntVar()
				self.cb_aa[cj] = Checkbutton(Cluster_page, text=str(i+1), variable=var)
				self.cb_aa[cj].grid(row=cj+1, column=0, sticky=W, padx=5)
				self.v_aa.append(var)
				self.cb_aa[cj].deselect()

				self.e_SERVER_nodename[cj] = Tkinter.Entry(Cluster_page, width=25)
				self.e_SERVER_nodename[cj].insert(0, "")
				self.e_SERVER_nodename[cj].grid(row=cj+1, column=1, sticky=W, padx = 5, pady=1)

				self.e_SERVER_max_proc[cj] = Tkinter.Entry(Cluster_page, width=5)
				self.e_SERVER_max_proc[cj].insert(0, "")
				self.e_SERVER_max_proc[cj].grid(row=cj+1, column=2, sticky=W, padx = 5, pady=1)

				self.e_SERVER_ssh_port[cj] = Tkinter.Entry(Cluster_page, width=5)
				self.e_SERVER_ssh_port[cj].insert(0, "")
				self.e_SERVER_ssh_port[cj].grid(row=cj+1, column=3, sticky=W, padx = 5, pady=1)

				self.e_SERVER_queue_name[cj] = Tkinter.Entry(Cluster_page, width=15)
				self.e_SERVER_queue_name[cj].insert(0, "")
				self.e_SERVER_queue_name[cj].grid(row=cj+1, column=4, sticky=W, padx = 5, pady=1)

				self.e_SERVER_home_dir[cj] = Tkinter.Entry(Cluster_page, width=25)
				self.e_SERVER_home_dir[cj].insert(0, "")
				self.e_SERVER_home_dir[cj].grid(row=cj+1, column=5, sticky=W, padx = 5, pady=1)

				self.e_SERVER_amber_dir[cj] = Tkinter.Entry(Cluster_page, width=30)
				self.e_SERVER_amber_dir[cj].insert(0, "")
				self.e_SERVER_amber_dir[cj].grid(row=cj+1, column=6, sticky=W, padx = 5, pady=1)
			cj += 1
			



		nb1.pack(fill='both',expand=1,padx=5,pady=5)
		nb1.setnaturalsize()

#		self.dialog.active()
		self.dialog.activate(geometry = 'centerscreenalways')



	def reset(self, textcur, filecur):
		textcur.delete(0, END)
		textcur.insert(0, filecur)

	def writeNewFile(self):
		path_home = os.environ.get('HOME')
		filename = "%s/Settings_Linux.txt" % path_home
            
		# write Settings_Linux.txt file
		fi = open(filename, "w")
		fi.write("USER:\n")
		fi.write("username                           %s\n" % username)
		fi.write("\nCLIENT:\n")
		fi.write("library_dir                        %s\n" % library_dir)
		fi.write("amber_dir                          %s\n" % amber_dir)
		fi.write("slide_main_dir                     %s\n" % slide_main_dir)
		fi.write("slide_exe_dir                      %s\n" % slide_exe_dir)
		fi.write("slide_main_dir_multiple            %s\n" % slide_main_dir_multiple)
		fi.write("slide_exe_dir_multiple             %s\n" % slide_exe_dir_multiple)
		fi.write("slide_data_dir                     %s\n" % slide_data_dir)
		fi.write("slide_dir                          %s\n" % slide_dir)
		fi.write("slide_sort                         %s\n" % slide_sort)
		fi.write("sie_dir                            %s\n" % sie_dir)
		fi.write("autodock_dir                       %s\n" % autodock_dir)
		fi.write("autodock_exe_dir                   %s\n" % autodock_exe_dir)
		fi.write("vina_dir                           %s\n" % vina_dir)
		fi.write("changepdb_dir                      %s\n" % changepdb_dir)
		fi.write("changecrd_dir                      %s\n" % changecrd_dir)
		fi.write("symposar_exe_dir                   %s\n" % symposar_exe_dir)
		fi.write("raptor_exe_dir                     %s\n" % raptor_exe_dir)
		fi.write("reduce_exe_dir                     %s\n" % reduce_exe_dir)
		fi.write("\n")
		fi.write("STANDARD SERVER SETTINGS:\n")
		fi.write("SERVER_number_of                   %d\n" % SERVER_number_of)
		fi.write("SERVER_nodename")
		for cj in range(SERVER_number_of):
			fi.write("     %s" % SERVER_nodename[cj])
		fi.write("\n")
		fi.write("SERVER_max_proc")
		for cj in range(SERVER_number_of):
			fi.write("     %s" % SERVER_max_proc[cj])
		fi.write("\n")
		fi.write("SERVER_ssh_port")
		for cj in range(SERVER_number_of):
			fi.write("     %s" % SERVER_ssh_port[cj])
		fi.write("\n")
		fi.write("SERVER_queue_name")
		for cj in range(SERVER_number_of):
			fi.write("     %s" % SERVER_queue_name[cj])
		fi.write("\n")
		fi.write("SERVER_home_dir")
		for cj in range(SERVER_number_of):
			fi.write("     %s" % SERVER_home_dir[cj])
		fi.write("\n")
		fi.write("SERVER_amber_dir")
		for cj in range(SERVER_number_of):
			fi.write("     %s" % SERVER_amber_dir[cj])
		fi.write("\n")
		fi.close()
        
	def searchFOLDER(self, filecur, textcur):
		self.folder_select = ""
		self.folder_select = askdirectory(title="Select directory", initialdir=filecur, mustexist=1)
		if self.folder_select:
			textcur.delete(0, END)
			textcur.insert(0, self.folder_select)
			filecur = self.folder_select

	def apply(self, result):
		global username
    
		global library_dir
		global amber_dir
		global gro_amb_dir
		global slide_main_dir
		global slide_exe_dir
		global slide_main_dir_multiple
		global slide_exe_dir_multiple
		global slide_data_dir
		global slide_dir
		global slide_sort
		global autodock_dir
		global autodock_exe_dir
		global vina_dir
		global sie_dir
		global changepdb_dir
		global changecrd_dir
		global symposar_exe_dir
		global raptor_exe_dir
		global reduce_exe_dir
 	   
		global SERVER_number_of
		global SERVER_nodename
		global SERVER_max_proc
		global SERVER_ssh_port
		global SERVER_queue_name
		global SERVER_home_dir
		global SERVER_amber_dir

		if result == 'Exit':
			self.dialog.deactivate()
			self.dialog.withdraw()
		elif result == 'Save to file':
			username                 = self.e_username.get()
			library_dir              = self.e_library_dir.get()
			amber_dir                = self.e_amber_dir.get()
#			gro_amb_dir              = self.e_gro_amb_dir.get()
			gro_amb_dir              = ""
			slide_main_dir           = self.e_slide_main_dir.get()
			slide_exe_dir            = self.e_slide_exe_dir.get()
			slide_main_dir_multiple  = self.e_slide_main_dir_multiple.get()
			slide_exe_dir_multiple   = self.e_slide_exe_dir_multiple.get()
			slide_data_dir           = self.e_slide_data_dir.get()
			slide_dir                = self.e_slide_dir.get()
			slide_sort               = self.e_slide_sort.get()
			sie_dir                  = self.e_sie_dir.get()
			autodock_dir             = self.e_autodock_dir.get()
			autodock_exe_dir         = self.e_autodock_exe_dir.get()
			vina_dir                 = self.e_vina_dir.get()
			changepdb_dir            = self.e_changepdb_dir.get()
			changecrd_dir            = self.e_changecrd_dir.get()
			symposar_exe_dir         = self.e_symposar_exe_dir.get()
			raptor_exe_dir           = self.e_raptor_exe_dir.get()
			reduce_exe_dir           = self.e_reduce_exe_dir.get()
			# number of clusters
			cj = 0
			for i in range(10):
				if self.v_aa[i].get() == 1:
					cj+=1
			SERVER_number_of = cj
			for cj in range(10):
				SERVER_nodename[cj] = self.e_SERVER_nodename[cj].get()
				SERVER_max_proc[cj] = self.e_SERVER_max_proc[cj].get()
				SERVER_ssh_port[cj] = self.e_SERVER_ssh_port[cj].get()
				SERVER_queue_name[cj] = self.e_SERVER_queue_name[cj].get()
				SERVER_home_dir[cj] = self.e_SERVER_home_dir[cj].get()
				SERVER_amber_dir[cj] = self.e_SERVER_amber_dir[cj].get()

			self.writeNewFile()			
		else:
			path_home                = os.environ.get('HOME')
			username                 = os.getlogin()
			library_dir              = "/usr/local/AMBER_library"
			amber_dir                = "/usr/local/amber10"
			gro_amb_dir              = ""
			slide_main_dir           = "/usr/local/slide"
			slide_exe_dir            = "/usr/local/slide/bin"
			slide_main_dir_multiple  = "/usr/local/slide_all"
			slide_exe_dir_multiple   = "/usr/local/slide_all/bin"
			slide_data_dir           = "%s/slide" % path_home
			slide_dir                = "%s/SLIDE_library" % path_home
			slide_sort               = "/usr/local/sortSlideResults"
			sie_dir                  = "/usr/local/Brimm"
			autodock_dir             = "%s/AUTODOCK_library" % path_home
			autodock_exe_dir         = "/usr/local/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24"
			vina_dir                 = "/usr/local/autodock_vina_1_1_2_linux_x86/bin"
			changepdb_dir            = "/usr/local/Ryde/changepdb"
			changecrd_dir            = "/usr/local/Ryde/changecrd"
			symposar_exe_dir         = "/usr/local/symposar"
			raptor_exe_dir           = "/usr/local/raptor"
			reduce_exe_dir           = "/usr/local/reduce/reduce.3.14.080821.src/reduce_src"
			SERVER_number_of         = 2    
			SERVER_nodename          = ["mycluster1.myuniversity.edu", "mycluster2.myuniversity.edu", "", "", "", "", "", "", "", ""]
			SERVER_max_proc          = ["8", "32", "", "", "", "", "", "", "", ""]
			SERVER_ssh_port          = ["22", "22", "", "", "", "", "", "", "", ""]
			SERVER_queue_name          = ["queue1", "queue2", "", "", "", "", "", "", "", ""]
			SERVER_home_dir          = ["/home_on_cluster1/username", "/home_on_cluster2/username", "", "", "", "", "", "", "", ""]
			SERVER_amber_dir         = ["/application_dir_on_cluster1/amber10", "/application_dir_on_cluster2/amber10", "", "", "", "", "", "", "", ""]
			self.reset(self.e_username, username)
			self.reset(self.e_library_dir, library_dir)
			self.reset(self.e_amber_dir, amber_dir)
			self.reset(self.e_slide_main_dir, slide_main_dir)
			self.reset(self.e_slide_exe_dir, slide_exe_dir)
			self.reset(self.e_slide_main_dir_multiple, slide_main_dir_multiple)
			self.reset(self.e_slide_exe_dir_multiple, slide_exe_dir_multiple)
			self.reset(self.e_slide_data_dir, slide_data_dir)
			self.reset(self.e_slide_dir, slide_dir)
			self.reset(self.e_slide_sort, slide_sort)
			self.reset(self.e_sie_dir, sie_dir)
			self.reset(self.e_autodock_dir, autodock_dir)
			self.reset(self.e_autodock_exe_dir, autodock_exe_dir)
			self.reset(self.e_vina_dir, vina_dir)
			self.reset(self.e_changepdb_dir, changepdb_dir)
			self.reset(self.e_changecrd_dir, changecrd_dir)
			self.reset(self.e_symposar_exe_dir, symposar_exe_dir)
			self.reset(self.e_raptor_exe_dir, raptor_exe_dir)
			self.reset(self.e_reduce_exe_dir, reduce_exe_dir)
			for cj in range(SERVER_number_of):
				self.reset(self.e_SERVER_nodename[cj], SERVER_nodename[cj])
				self.reset(self.e_SERVER_max_proc[cj], SERVER_max_proc[cj])
				self.reset(self.e_SERVER_ssh_port[cj], SERVER_ssh_port[cj])
				self.reset(self.e_SERVER_queue_name[cj], SERVER_queue_name[cj])
				self.reset(self.e_SERVER_home_dir[cj], SERVER_home_dir[cj])
				self.reset(self.e_SERVER_amber_dir[cj], SERVER_amber_dir[cj])

  
#######################################################################################
def modify_settings(app):
	mfs = modifySettingsFile(app.root)


#######################################################################################
def read_settings(app):
	global username
    
	global library_dir
	global amber_dir
	global gro_amb_dir
	global slide_main_dir
	global slide_exe_dir
	global slide_main_dir_multiple
	global slide_exe_dir_multiple
	global slide_data_dir
	global slide_dir
	global slide_sort
	global autodock_dir
	global autodock_exe_dir
	global vina_dir
	global sie_dir
	global changepdb_dir
	global changecrd_dir
	global symposar_exe_dir
	global raptor_exe_dir
	global reduce_exe_dir
    
	global SERVER_number_of

	path1 = os.environ.get('HOME')
	filename = "%s/Settings_Linux.txt" % (path1)
	try:
		fi = open(filename, "r")
	except:
		tkMessageBox.showwarning("Settings file not found", "Settings_Linux.txt not found. Please, copy download default file and copy to $HOME directory.")
		sys.exit()

#    print "Read settings for SLIDE"
	for i in fi:
		if i.find("username") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			username = dat2.rstrip()

		if i.find("library_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			library_dir = dat2.rstrip()
		if i.find("amber_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			amber_dir = dat2.rstrip()
		if i.find("gro_amb_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			gro_amb_dir = dat2.rstrip()
		if i.find("slide_main_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			slide_main_dir = dat2.rstrip()
		if i.find("slide_exe_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			slide_exe_dir = dat2.rstrip()
		if i.find("slide_main_dir_multiple") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			slide_main_dir_multiple = dat2.rstrip()
		if i.find("slide_exe_dir_multiple") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			slide_exe_dir_multiple = dat2.rstrip()
		if i.find("slide_data_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			slide_data_dir = dat2.rstrip()
		if i.find("slide_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			slide_dir = dat2.rstrip()
		if i.find("slide_sort") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			slide_sort = dat2.rstrip()
		if i.find("autodock_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			autodock_dir = dat2.rstrip()
		if i.find("autodock_exe_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			autodock_exe_dir = dat2.rstrip()
		if i.find("vina_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			vina_dir = dat2.rstrip()
		if i.find("sie_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			sie_dir = dat2.rstrip()
		if i.find("changepdb_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			changepdb_dir = dat2.rstrip()
		if i.find("changecrd_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			changecrd_dir = dat2.rstrip()
		if i.find("symposar_exe_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			symposar_exe_dir = dat2.rstrip()
		if i.find("raptor_exe_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			raptor_exe_dir = dat2.rstrip()
		if i.find("reduce_exe_dir") == 0:
			tmpc, dat = i.split(None, 1)
			dat2 = str(dat)
			reduce_exe_dir = dat2.rstrip()

		if i.find("SERVER_number_of") == 0:
			tmpc, dat = i.split(None, 1)
			SERVER_number_of = int(dat)

	global SERVER_nodename
	global SERVER_max_proc
	global SERVER_ssh_port
	global SERVER_queue_name
	global SERVER_home_dir
	global SERVER_amber_dir
	SERVER_nodename = ["" for i in range(10)]
	SERVER_max_proc = ["0" for i in range(10)]
	SERVER_ssh_port = ["" for i in range(10)]
	SERVER_queue_name = ["" for i in range(10)]
	SERVER_home_dir = ["" for i in range(10)]
	SERVER_amber_dir = ["" for i in range(10)]
 
	fi.seek(0)
	for i in fi:
		if i.find("SERVER_nodename") == 0:
			i2 = i.strip()
			datarray = i2.split(None)
			datarray.pop(0)
			cj = 0
			for x in datarray:
				SERVER_nodename[cj] = str(x)
				cj += 1
		if i.find("SERVER_max_proc") == 0:
			i2 = i.strip()
			datarray = i2.split(None)
			datarray.pop(0)
			cj = 0
			for x in datarray:
				SERVER_max_proc[cj] = str(x)
				cj += 1
		if i.find("SERVER_ssh_port") == 0:
			i2 = i.strip()
			datarray = i2.split(None)
			datarray.pop(0)
			cj = 0
			for x in datarray:
				SERVER_ssh_port[cj] = str(x)
				cj += 1
		if i.find("SERVER_queue_name") == 0:
			i2 = i.strip()
			datarray = i2.split(None)
			datarray.pop(0)
			cj = 0
			for x in datarray:
				SERVER_queue_name[cj] = str(x)
				cj += 1
		if i.find("SERVER_home_dir") == 0:
			i2 = i.strip()
			datarray = i2.split(None)
			datarray.pop(0)
			cj = 0
			for x in datarray:
				SERVER_home_dir[cj] = str(x)
				cj += 1
		if i.find("SERVER_amber_dir") == 0:
			i2 = i.strip()
			datarray = i2.split(None)
			datarray.pop(0)
			cj = 0
			for x in datarray:
				SERVER_amber_dir[cj] = str(x)
				cj += 1
    
	fi.close()



#read_settings(self)
