import csv
import os
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pylab import *
from matplotlib.font_manager import FontProperties
import math

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect("equal")

os.chdir('/home/jacob/Desktop/RESULT/')

with open('DOC sheet.csv','rb') as csvfile:
	reader = csv.reader(csvfile)
	#print reader
	family = []
	conserved_plr = []
	conserved_nonplr = []
	for row in reader:
		if 'TRANSPORT FAMILY' not in row:
			new = row[2].strip()
			family.append(new)
			family = [x for x in family if x!= '']
			new1= row[6].strip()
			conserved_plr.append(new1)
			conserved_plr = [x for x in conserved_plr if x!= '']
			new2= row[10].strip()
			conserved_nonplr.append(new2)
			conserved_nonplr = [x for x in conserved_nonplr if x!= '']


transport_ids = len(family)
index = np.arange(transport_ids)

width = 0.35

con_plr = sorted([float(x) for x in  conserved_plr])
con_nonplr = sorted([float(x) for x in conserved_nonplr])
print len(con_plr), len(con_nonplr)
mean_con_plr = sum(con_plr)/float(len(con_plr))
mean_con_nonplr = sum(con_nonplr)/float(len(con_nonplr))
std_con_plr =  np.std(con_plr)
std_con_nonplr = np.std(con_nonplr)	
#distance_between_means = 0.5/np.std(con_nonplr)

rects1 = ax.bar(index, con_plr, width, color = 'black' , yerr = std_con_plr, error_kw= dict(elinewidth = 2, ecolor = 'red'))
rects2 = ax.bar(index+width, con_nonplr, width, color= 'red', yerr = std_con_nonplr, error_kw = dict(elinewidth = 2, ecolor = 'black'))
		

ax.set_xlim(-width, len(index)+width)
ax.set_ylim(0,5)

plt.title('HYPOTHESIS TESTING ')
plt.ylabel('MEAN DOC for CONSERVED PLRs/nonPLRs', fontsize=14)

xTickmarks = [x for x in  family]
ax.set_xticks(index+width)
xtickNames = ax.set_xticklabels(xTickmarks)
plt.setp(xtickNames, rotation = 90, fontsize = 10)

ax.legend((rects1[0],rects2[0]), ('CONSERVED PLRs','CONSERVED NON-PLRs'))

#plt.show()

ax.set_ylim(0,50)
plt.title('HYPOTHESIS TESTING HISTOGRAM')
plt.ylabel('MEAN DOC', fontsize=14)
plt.xlabel('Bins', fontsize = 14)
hB, = plot([1, 1],'g-')
hR, = plot([1, 1],'b-')
fontP = FontProperties()
fontP.set_size('small')
legend((hB, hR),('NONPLRs', 'PLRs'), prop = fontP)
hB.set_visible(False)
hR.set_visible(False)
con_fit = stats.norm.pdf(con_plr, np.mean(con_plr), np.std(con_plr))
noncon_fit = stats.norm.pdf(con_nonplr, np.mean(con_nonplr), np.std(con_nonplr))

con_fit = stats.wilcoxon(con_plr,con_nonplr,zero_method='wilcox', correction=False)[1]
noncon_fit = stats.wilcoxon(con_nonplr, con_plr,zero_method='wilcox', correction=False)[1]
print con_fit, noncon_fit

plt.plot(con_fit,'-o')
plt.plot(noncon_fit,'-o')
plt.hist(con_plr, bins=20,normed = True)
plt.hist(con_nonplr, bins= 20, normed= True)
#plt.scipy.stats.wilcoxon(con_plr, bins=20,normed = True)
#plt.scipy.stats.wilcoxon(con_nonplr, bins= 20, normed= True)

#plt.show()

#print mean_con_plr,mean_con_nonplr
#print std_con_plr, std_con_nonplr


#con_plr_overlap = 0.5 *(1 + math.erf((min(con_plr) - mean_con_nonplr)/(std_con_nonplr * (2**0.5))))
#con_nonplr_overlap = 1- (0.5 *(1 + math.erf((max(con_nonplr) - mean_con_plr)/(std_con_plr * (2**0.5)))))
 
#print con_plr_overlap, con_nonplr_overlap

