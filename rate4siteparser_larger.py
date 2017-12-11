import os
import re
import warnings
import scipy.stats as stats
import matplotlib.pyplot as plt
import pylab as plt2
from pylab import *
import numpy as np
from decimal import *
from scipy.stats import norm
from sklearn import preprocessing
from math import floor, ceil 


print "		step 2.4 : Rate4Site Parser and Plotting the Bins"
   
def rate4site():
	
	curdir = os.getcwd()
	curdir_up = '/'.join(curdir.split('/')[:-1])
	rate_path = "%s/RESULT/MODULE2-LARGER" %curdir_up
	
	
	for curdir, subdirs, files in os.walk(rate_path): 
		
		for file in files[:]:
			
			
			if '_rate' in file: 
					
				
				query = file.split('_')[0]
				outpath = os.path.join(curdir, file)
				rate_file = open(outpath)
				confidence = []
				data_reqs = []
				fil_data = []
				scores = [] 
				residues = []
				range1 = []
				range2=[]
				y_count = []
				
				for i,interval in enumerate(rate_file):
					
					if i == 0:
						
						True
						
					else:
						
						if not interval.startswith('#') and interval != '\n':
							
							interval = interval.split(',')
							parts1 = interval[0].split()
							seq = parts1[0]
							res = parts1[1]
							residues.append(res)
							y_count.append(seq)
							score = parts1[2]				
							scores.append(float(score))					
							range1.append(parts1[3][1:])
							range1 = [item for item in range1 if item.strip()]					
							parts2 = interval[-1].split()
							range2.append(parts2[0][:-1])
							range2 = [item for item in range2 if item.strip()]
							range_diff = range2 + range1
							confidence.append(interval)
							data_reqs.append([score,range_diff])
					
					
				fig, ax = plt.subplots(1, 1)
				mean, var, skew, kurt = norm.stats(moments='mvsk')
				rang = [min(scores), max(scores)]
				Long = len(scores)
				Maxim = max(scores) #MaxValue
				Minim = min(scores) #MinValue
				av = np.mean(scores) #Average
				StDev = np.std(scores) #Standard Dev.
				
				
				x = np.linspace(Minim, Maxim, Long)
				ax.plot(x, norm.pdf(x, av, StDev),'r-', lw=3, alpha=0.9, label='RATE SCORES')
				
				weights = np.ones_like(scores)/len(scores)
				normalized = [(s-min(scores))/(max(scores)-min(scores)) for s in scores]
				newpath = os.path.join(curdir)


				ax.hist(normalized, weights = weights, normed=True, histtype='stepfilled', alpha=0.2,label='NORMALIZED RATE SCORES')
				plt.title('%s' %query +'_Normalized_Rate4Site_Scores')
				plt.xlabel('Rate4Site_Scores', fontsize=14)
				plt.ylabel('Sequence_Count', fontsize=14)
				plt.legend(loc='upper right')
				fig.savefig(newpath +'_%s_normalized.png'%query)
				plt.close("all")

				y_count = map(int, y_count)
				
				color_path = os.path.join(curdir)
				
				bins = np.arange(floor(min(normalized[:])), ceil(max(normalized[:])), 0.10)
				
				colors = ('blue', 'red', 'green','cyan','purple','pink','violet','lime','aqua')
				
				# get the max count for a particular bin for all classes combined
				max_bin = max(np.histogram(normalized[:], bins=bins)[0])
				plt.figure()
				n, bins, patches = plt.hist(normalized[:], bins , alpha=0.3)
				
				for c, p in zip(colors, patches):
					
					plt.setp(p, 'facecolor', c)  
					
				plt.ylim([0, max_bin*1.3])
				plt.title('%s'%query + '_Normalized_Scores_In_9_Bins')
				plt.xlabel('Color_Bins', fontsize=14)
				plt.ylabel('Sequence_Count', fontsize=14)
				plt.legend(loc='upper right')
				plt.savefig(newpath + '_%s_ColorBins.png' %query)
				plt.close("all")
			
			
					
				pairs = [(x,y,z) for x,y,z in zip(normalized, residues, y_count)]
				group1 = []
				group2 = []
				group3 = []
				
				for item in pairs[:]:
					
					if item[0]<= 0.3:
						
						group1.append(item)
						
					elif 0.4 <= item[0] <=0.7:
						
						group2.append(item)
						
					else:
						
						group3.append(item)	
				
					
				conserved_residues = open(newpath + '/%s_conserved_residues.txt' %query, 'w')

				respo = []
				
				for items in group1[:]:
					
					res_po = tuple((items[2],items[1]))
					respo.append(res_po)
					
				strs =" ".join(str(x) for x in respo)
				conserved_residues.write(strs+"\n")
				conserved_residues.close()
				#return respo	
				
rate4site()



