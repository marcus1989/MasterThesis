import csv
import os
from scipy.stats.stats import pearsonr 
from scipy.stats import hypergeom
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pylab import *
from matplotlib.font_manager import FontProperties
import math
import numpy 
import seaborn as sns
import pandas as pd

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect("equal")

folder_path = '/home/jacob/Desktop/RESULT/'

groups = []

for root, subfolders, files in os.walk(folder_path):
	
	for folders in subfolders:
		groups.append(folders)


group1_larger = folder_path + groups[2]
group2_smaller = folder_path + groups[5]

grp_smaller = []

for root, subfolders, files in os.walk(group2_smaller):
	
	for fils in files:
		
		if 'plr_pos' in fils and fils != 'plr_pos.png':
			
			map_files_smaller = root + '/' + fils
			grp_smaller.append(map_files_smaller)
			


grp_larger = []
for root, subfolders , files in os.walk(group1_larger):
	
	for fils in files:
		
		if 'plr_pos' in fils and fils != 'plr_pos.png' :
			
			map_files_larger = root + '/' + fils
			grp_larger.append(map_files_larger)

whole_set = grp_smaller +grp_larger
print len(grp_smaller),len(grp_larger),len(whole_set)


pvalue = []
hypergom = []
all_names = []

def comp_hyper(pop, plrs, cons, plr_con):
	result = 0
	for i in range(plr_con,pop+1):
		rv = hypergeom(pop, plrs, cons)
		result += rv.pmf(i)
	return result 

for plr in whole_set:
	
	plrfiles = plr
	plr_pos = open(plrfiles).readlines()

	xy = plr_pos[-3:-1]
	pdb  =plr_pos[-1][:-1]

	row_c =[]
	row_p =[]
	len_r = []
	no_c = []
	no_p = []
	no_cp = []
	
	names = plrfiles.split('/')[8]
	all_names.append(names)
	for i,item in enumerate(pdb):
		
		len_r.append(i)
		c,p = xy[0][i], xy[1][i]
		#print c,p,item,plrfiles
		if c == 'c':
			row_c.append(int(1))
			no_c.append(int(1))
	
		elif c == '*':
			row_c.append(int(0))
									
		if p == 'p':
			row_p.append(int(1))
			no_p.append(int(1))
				
		elif p == '*':
			row_p.append(int(0))
		
		if c == 'c' and p == 'p':
			no_cp.append(int(1))
		
		else:
			
			"nothing"
	
	#print len(len_r), len(no_p), len(no_c),len(no_cp),plrfiles
	hypergom.append(comp_hyper(len(len_r), len(no_p), len(no_c), len(no_cp)))
	##print row_c,row_p
	#pearson_corr = pearsonr(row_c,row_p)
	#pearson.append(pearson_corr[0])
	#pvalue.append(pearson_corr[1])
	#print pearson_corr, plrfiles

print(hypergom)
#print(all_names)
#plt.title('Hypergeometric test')
#plt.xlabel('CORRELATION VALUES', fontsize=10)
#plt.ylabel('P-VALUES', fontsize=10)
#ax.axhline(0, color='red', lw=0.5)
sns.set_style("whitegrid")

data=pd.DataFrame({'Protein':all_names, 'P-value':hypergom})
g=sns.FacetGrid(data,aspect=2,size=6)
g.map(sns.barplot, 'Protein','P-value')

plt.show()
	
		
