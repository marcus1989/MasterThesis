import os
import glob
import re
import scipy.stats as stats
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
from decimal import *
from scipy.stats import norm
from sklearn import preprocessing
from math import floor, ceil 

folder_path = '/home/jacob/Desktop/RESULT/COMPARE_LARGER'


groups = []
for root, subfolders, files in os.walk(folder_path):
	
	for folders in subfolders:
		groups.append(folders)


group1 = folder_path+ '/'+ groups[0]
group2 = folder_path + '/' + groups[1]
group3 = folder_path + '/' + groups[2]

points_path = folder_path + '/Points'
points = []
for root, subfolders, files in os.walk(points_path):
	
	for folders in subfolders:
		points.append(folders)
		
points1 = points_path + '/' + points[0] 
points2 = points_path + '/' + points[1]
points3 = points_path + '/' + points[2]
scatter = [points1,points2,points3]


conserved_kyte_scores = []
conserved_rate_scores = []
for items in scatter[:]:
	
	for root, folders, fil in os.walk(items):
		
		for files in fil:
			
			#print files
			counter_1 = 0
			if 'PLR' in files:
				
				plrfile = open(items + '/' + files).readlines()
				con_kyte = []
				con_rate = []
				for con in plrfile:
					
					consk = con.split()[2]
					consk = consk.strip(',')
					conskyte = float(consk)
					con_kyte.append(conskyte)
					consrate = con.split()[3]
					consr = consrate.strip(',')
					consrate = float(consr)
					con_rate.append(consrate)
			
				conserved_kyte_scores.append(con_kyte)	
				conserved_rate_scores.append(con_rate)	
				
#print conserved_kyte_scores[0]

box_items = [group1, group2, group3]
#print box_items
kyte_scaled_scores = []
rate_scaled_scores=[]
kyte_mean_list = []
kyte_sd_list = []
rate_mean_list = []
rate_sd_list = []
for items in box_items[:]:
	#print items
	for root, folders, fil in os.walk(items):
		
		for files in fil:
			
			#print files
			counter_1 = 0
			if 'kyte' in files:
				
				#print files
				
				kytefile = open(items + '/' + files).readlines()
				#kytefile = open(files).readlines()
				#print items + '/' + files
				#print kytefile
				#break
				description = []
				for newitems in kytefile:
					#print newitems
					if newitems.startswith("Position"):
						
						descriptions = newitems.strip(" ").split()
						description.append(descriptions)
				#print description
				##break		
				
				kytescores = []
				for kscores in description[:]:
					
					kdscores = float(kscores[3])
					kytescores.append(kdscores)
				
				#print min(kytescores), max(kytescores)
				##break
				
				newscores = [(x - min(kytescores))/(max(kytescores) - min(kytescores)) for x in kytescores]	
				kyte_scaled_scores.append(kytescores)
				#print kyte_scaled_scores
				##break				
				av_kyte = np.mean(kytescores)
				kyte_mean_list.append(av_kyte)
				StDev_kyte = np.std(kytescores)
				kyte_sd_list.append(StDev_kyte)
				
				#print av_kyte,StDev_kyte
				##break
			else:
				pass

#print box_items
for items in box_items[:]:
	#print items
	for root, folders, fil in os.walk(items):
		
		for files in fil:
			
			#print files
			counter_1 = 0
			if 'RATE' in files:	
							
				rate_file = open(items+ '/' + files).readlines()
				##print rate_file
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
							##print res
							y_count.append(seq)
							score = parts1[2]		
							scores.append(float(score))					
				#print scores
				
				#print min(scores),max(scores)
				normalized = [(s-min(scores))/(max(scores)-min(scores)) for s in scores]
				rate_scaled_scores.append(scores)



				av_rate = np.mean(scores) 
				rate_mean_list.append(av_rate)
				StDev_rate = np.std(scores) #Standard Dev.
				rate_sd_list.append(StDev_rate)
				
				#print av_rate,StDev_rate

def setBoxColors(bp):
	
    setp(bp['boxes'][0], color='blue')
    setp(bp['caps'][0], color='blue')
    setp(bp['caps'][1], color='blue')
    setp(bp['whiskers'][0], color='blue')
    setp(bp['whiskers'][1], color='blue')
    setp(bp['fliers'][0], color='blue')
    setp(bp['fliers'][1], color='blue')
    setp(bp['medians'][0], color='blue')

    setp(bp['boxes'][1], color='red')
    setp(bp['caps'][2], color='red')
    setp(bp['caps'][3], color='red')
    setp(bp['whiskers'][2], color='red')
    setp(bp['whiskers'][3], color='red')
    setp(bp['fliers'][2], color='red')
    setp(bp['fliers'][3], color='red')
    setp(bp['medians'][1], color='red')

ABC = [kyte_scaled_scores[0],rate_scaled_scores[0]]
FATP = [kyte_scaled_scores[1],rate_scaled_scores[1]]
PATP = [kyte_scaled_scores[2],rate_scaled_scores[2]] 

#print PATP
data_to_plot = [ABC,FATP,PATP]

fig = plt.figure()
ax = axes()
hold = True

bp = boxplot(ABC,positions = [1.5, 2.5])
plt.scatter([1.5], [kyte_mean_list[0]],color = 'blue')
plt.scatter([1.5],[kyte_mean_list[0]-kyte_sd_list[0]],s=80,color = 'blue', marker=(5, 2))
plt.scatter([1.5],[kyte_mean_list[0]+kyte_sd_list[0]],s=80,color = 'blue', marker=(5, 2))



plt.scatter([2.5], [rate_mean_list[0]],color = 'red')
plt.scatter([2.5],[rate_mean_list[0]-rate_sd_list[0]],s=80,color = 'red', marker=(5, 2))
plt.scatter([2.5],[rate_mean_list[0]+rate_sd_list[0]],s=80,color = 'red', marker=(5, 2))

setBoxColors(bp)

bp = boxplot(FATP, positions = [4.5, 5.5])
plt.scatter([4.5], [kyte_mean_list[0]],color = 'blue')
plt.scatter([4.5],[kyte_mean_list[0]-kyte_sd_list[0]],s=80,color = 'blue', marker=(5, 2))
plt.scatter([4.5],[kyte_mean_list[0]+kyte_sd_list[0]],s=80,color = 'blue', marker=(5, 2))

plt.scatter([5.5], [rate_mean_list[0]],color = 'red')
plt.scatter([5.5],[rate_mean_list[0]-rate_sd_list[0]],s=80,color = 'red', marker=(5, 2))
plt.scatter([5.5],[rate_mean_list[0]+rate_sd_list[0]],s=80,color = 'red', marker=(5, 2))

setBoxColors(bp)

bp = boxplot(PATP, positions = [7.5, 8.5])
plt.scatter([7.5], [kyte_mean_list[0]],color = 'blue')
plt.scatter([7.5],[kyte_mean_list[0]-kyte_sd_list[0]],s=80,color = 'blue', marker=(5, 2))
plt.scatter([7.5],[kyte_mean_list[0]+kyte_sd_list[0]],s=80,color = 'blue', marker=(5, 2))

plt.scatter([8.5], [rate_mean_list[0]],color = 'red')
plt.scatter([8.5],[rate_mean_list[0]-rate_sd_list[0]],s=80,color = 'red', marker=(5, 2))
plt.scatter([8.5],[rate_mean_list[0]+rate_sd_list[0]],s=80,color = 'red', marker=(5, 2))

setBoxColors(bp)



for i in [1.5,2.5,4.5,5.5,7.5,8.5]:
	
	if i == 1.5:
		
		y = conserved_kyte_scores[:][0]
		x = np.random.normal(i, 0.04, size=len(y))
		colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(conserved_kyte_scores[:][0])))
		cs = [colors[i//len(x)] for i in range(len(conserved_kyte_scores[:][0])*len(x))]
		matplotlib.pyplot.scatter(x,y,color=cs)
		
	elif i == 2.5:
		
		y = conserved_rate_scores[:][0]
		x = np.random.normal(i, 0.04, size=len(y))
		colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(conserved_rate_scores[:][0])))
		cs = [colors[i//len(x)] for i in range(len(conserved_rate_scores[:][0])*len(x))]
		matplotlib.pyplot.scatter(x,y,color=cs)
		
	elif i == 4.5:
		
		y = conserved_kyte_scores[:][1]
		x = np.random.normal(i, 0.04, size=len(y))
		colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(conserved_kyte_scores[:][1])))
		cs = [colors[i//len(x)] for i in range(len(conserved_kyte_scores[:][1])*len(x))]
		matplotlib.pyplot.scatter(x,y,color=cs)
		
	elif i == 5.5:
		
		y = conserved_rate_scores[:][1]
		x = np.random.normal(i, 0.04, size=len(y))
		colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(conserved_rate_scores[:][1])))
		cs = [colors[i//len(x)] for i in range(len(conserved_rate_scores[:][1])*len(x))]
		matplotlib.pyplot.scatter(x,y,color=cs)
		
	elif i == 7.5:
		
		y = conserved_kyte_scores[:][2]
		x = np.random.normal(i, 0.04, size=len(y))
		colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(conserved_kyte_scores[:][2])))
		cs = [colors[i//len(x)] for i in range(len(conserved_kyte_scores[:][2])*len(x))]
		matplotlib.pyplot.scatter(x,y,color=cs)
		
	elif i ==8.5:
		
		y = conserved_rate_scores[:][2]
		x = np.random.normal(i, 0.04, size=len(y))
		colors = matplotlib.cm.rainbow(np.linspace(0, 1, len(conserved_rate_scores[:][2])))
		cs = [colors[i//len(x)] for i in range(len(conserved_rate_scores[:][2])*len(x))]
		matplotlib.pyplot.scatter(x,y,color=cs)
						
	else:
		
		pass

ax.set_ylim(-5.0,5.0)
ax.set_xlim(0,9)

plt.ylabel('KYTE/RATE SCORES', fontsize=14)
ax.set_xticklabels(['ABC', 'FATP', 'PATP'])
ax.set_xticks([2.0, 5.0, 8.0])

# draw temporary red and blue lines and use them to create a legend
hB, = plot([1, 1],'b-')
hR, = plot([1, 1],'r-')

fontP = FontProperties()
fontP.set_size('large')
legend((hB, hR),('KYTE-DOOLITTLE SCORES', 'RATE SCORES'), prop = fontP)
hB.set_visible(False)
hR.set_visible(False)

plt.show()
