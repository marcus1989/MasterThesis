#!/usr/bin/python
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from itertools import product, combinations
import os
import re
import math
from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain
import box
import glob
from os.path import basename
#import invokepores

print "		step 3.1 : Calculating the center of Mass of the pore and creating the BOX"

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
rate_path = "%s/RESULT/MODULE2" %curdir_up

for curdir, subdirs, files in os.walk(rate_path):
	
	for file in files[:]:
		
		if '.pdb' in file:
			
			new = file.split(".")[0]

			pdb = open(os.path.join(curdir , file),'r').readlines()

			charged_res = []
			x =[]
			z =[]
			y =[]
			com =[]
			z_cutoff = []
			
			for item in pdb:
				
				if item.startswith('ATOM'):

					xcoor = item[31:38]
					ycoor = item[39:46]
					zcoor = item[47:54]
					x.append(float(xcoor))
					y.append(float(ycoor))
					z.append(float(zcoor))
					
			for item in z[:]:
				
				if -15.0<= item and item <=15.0:
					z_cutoff.append(item)
					
					
			x_cord = sum(x)/len(x)
			
			y_cord = sum(y)/len(y)
			
			try:
				
				z_cord = sum(z_cutoff)/len(z_cutoff)
				
			except ZeroDivisionError:
				
				pass
			
			com.append([x_cord,y_cord,z_cord])
			out1 = open(os.path.join(curdir,new + 'out' + '_com ' + '.txt'),'w')
			out1.write(str(com))
			out1.close()
			
			one_comp = box.box((x_cord,y_cord,z_cord))
			out2 = open(os.path.join(curdir,new + 'out' + '_box ' + '.txt'),'w')
			out2.write(str(one_comp))
			out2.close()
	
		
			###########################3d-plot###################################
			
			centre = ax.scatter([x_cord],[y_cord],[z_cord],color="g",s=100)
			points = box.box((x_cord,y_cord,z_cord))
			scatterfile = ax.scatter(points[0], points[1], points[2],color = "r")
			#plt.show()
			
#invokepores.invoke()



