import os
import re
from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt 
from pylab import *
from weblogolib import *
import pylab as plt
from matplotlib.legend_handler import HandlerLine2D

fig = plt.figure()
ax = axes()
hold = True

with open('/home/jacob/Desktop/RESULT/MODULE2-LARGER/P-ATPase/SACOL2572/4bbj/plr_pos') as f:
	
	data = f.read()		
	
data = data.split()[3]
labels = [x for x in data]

kytefile = open('/home/jacob/Desktop/RESULT/MODULE2-LARGER/P-ATPase/SACOL2572/4bbj/kyte').readlines()
kyte_scores = {}
k_scores = []

for kscore in kytefile:
	
	if kscore.startswith("Position"):
		
		kytescore = float(kscore.strip(" ").split()[3])
		kyteindex = int(kscore.strip(" ").split()[1])
		k_scores.append(kytescore)
		
		if kyteindex not in kyte_scores:
			
			kyte_scores[kyteindex] = kytescore
			
		elif kyteindex in kyte_scores:
			
			kyte_scores[kyteindex].append(kytescore)
			

	
rate_file = open('/home/jacob/Desktop/RESULT/MODULE2-LARGER/P-ATPase/SACOL2572/SACOL2572_rate','r')
rate_dict = {}
rate_scores_dict ={} 
r_scores = [] 

query_len = 0
for i,interval in enumerate(rate_file):
	
	if i == 0:
		
		True
		
	else:
		
		if not interval.startswith('#') and interval != '\n':
			
			interval = interval.split(',')
			parts1 = interval[0].split()
			rate_key = 'SACOL2572'
			rate_pos = int(parts1[0])
			rate_res = parts1[1]
			score = float(parts1[2])				
			r_scores.append(float(score))	
			rate_tuple = tuple((rate_pos,rate_res))
			
			if rate_pos not in rate_dict:
				
				rate_dict[rate_pos] = rate_tuple
				
			elif rate_pos in rate_dict:
				
				rate_dict[rate_pos].append(rate_tuple)
			
			if rate_pos not in rate_scores_dict:
				
				rate_scores_dict[rate_pos] = score
				
			elif rate_pos in rate_scores_dict:
				
				rate_scores_dict[rate_pos].append(score) 	
				
#print k_scores
#print r_scores

clustal_files = open('/home/jacob/Desktop/RESULT/MODULE2-LARGER/P-ATPase/SACOL2572/4bbj/4BBJ:A|PDBID|CHAIN|SEQUENCEclustal_file').readlines()
pdb_tuples = []
for datas in clustal_files:
	
	data1 = datas.split(',')
	data2 = int(data1[0]),data1[1].strip('\n')
	pdb_tuples.append(data2)

pdb_tuple  = sorted(pdb_tuples)

rate_tuples  = sorted(rate_dict.values())

#print len(pdb_tuple),len(rate_tuples)

mapped_rate = {}
pdb_po_start = 0
rate_po_start = 0
rate_po_index = 0
stop = False

for tup in pdb_tuple:
	
	if stop == False:
		
		pdbpos,pdbres = tup
		for i,ratetup in enumerate(rate_tuples):
			ratepos, rateres = ratetup
			
			if pdbres == rateres:

				pdb_po_start = pdbpos
				rate_po_start = ratepos
				rate_po_index = i
				stop = True
				break		
				
#print pdb_tuple	
#print rate_tuples	

for tup1 in pdb_tuple:

	try: 
		rate_tuple = rate_tuples[rate_po_index]
		key = tup1[0]
		value = (tup1[0],tup1[1],rate_tuple[0],rate_tuple[1])
		mapped_rate[key] = value
		rate_po_index += 1
		
	except IndexError:
		
		break
				
#print kyte_scores
#print rate_scores_dict

graph_dict = {}

for maps in mapped_rate:
	
	vals = mapped_rate[maps]
	rate_key = vals[2]
	pdb_pos = vals[0]
	pdb_res = vals[1]
	graph_keys = pdb_pos
	
	if pdb_pos in kyte_scores:
		
		r,k= rate_scores_dict[rate_key], kyte_scores[pdb_pos]
		graph_vals = (pdb_pos,pdb_res,r,k)
		graph_dict[graph_keys] = graph_vals

#print graph_dict

fasta_filess = open('/home/jacob/Desktop/RESULT/MODULE2-LARGER/P-ATPase/SACOL2572/4bbj/4BBJ:A|PDBID|CHAIN|SEQUENCEfastafile3').readlines()
fasta_tuples = []

for fastadatas in fasta_filess:
	
	fasta_data1 = fastadatas.split(',')
	fasta_data2 = int(fasta_data1[0]),fasta_data1[1].strip('\n')
	fasta_tuples.append(fasta_data2)

#print fasta_tuples
conserved_files = open('/home/jacob/Desktop/RESULT/MODULE2-LARGER/P-ATPase/SACOL2572/4bbj/4BBJ:A|PDBID|CHAIN|SEQUENCEconserved_set_file').readlines()
conserved_tuples = []
#print conserved_files

for conserveddatas in conserved_files:
	
	conserveddata1 = conserveddatas.split(',')
	conserveddata2 = int(conserveddata1[0]),conserveddata1[1].strip('\n')
	conserved_tuples.append(conserveddata2)

plr_paired = []
not_plr_paired = []

for tup  in fasta_tuples:

	p,r = tup
	
	if '*' not in r:
		
		plr_pair = (p,r)
		plr_paired.append(plr_pair)
		
	else:
		
		not_plr_pair = (p,r)
		not_plr_paired.append(not_plr_pair)
				
#print len(plr_paired),len(not_plr_paired)

conserved = []

for tup1 in conserved_tuples:
	
	p2,r2 = tup1

	if '*' not in r2:
		
		cons_res = (p2,r2)
		conserved.append(cons_res)
		
#print conserved

conserved_plrset = []
not_conserved_plrset = []

for pls in plr_paired[:]:
	
	p_1,r_1 = pls
	
	for conss in conserved[:]:
		
		cp_1,cr_1 = conss
		
		if cp_1 == p_1:
			
			new_con_set = (cp_1,cr_1)
			conserved_plrset.append(new_con_set)
			
		else:
			
			pass	

conserved_nonplr = list(set(conserved)-set(conserved_plrset))

print  len(conserved), len(conserved_plrset), len(conserved_nonplr)

#print conserved_plrset
	
conserved_plr = open("/home/jacob/Desktop/RESULT/MODULE2-LARGER/P-ATPase/SACOL2572/conserved_plrset_4bbj.txt","w")

for point in graph_dict.keys():
	
	p1,r1,r_s,k_s = graph_dict[point]
	ps_list = []
	
	for pa in conserved_plrset:
		
		ps,rs = pa
		ps_list = []
		
		if ps == p1:
			
			ps_list.append((p1,r1,r_s,k_s))
			plt.scatter(p1,r_s,color='blue')
			plt.scatter(p1,k_s,color='green')
			for i, txt in enumerate(ps_list):
				txt2 = (list(txt[:2]) + [txt[3]])
		
				ax.annotate(txt[1],xy =(txt[0],txt[2]),xytext=(0.0,0.1),textcoords='offset points', ha='center', va='bottom',
			bbox=dict(boxstyle='round,pad=0.1', fc='yellow', alpha=0.1),size = 10)
							
				ax.annotate(txt[1],xy = (txt[0],txt[3]),xytext=(0.0,0.1),textcoords='offset points', ha='center', va='bottom',
			bbox=dict(boxstyle='round,pad=0.1', fc='yellow', alpha=0.1),size = 10)

			for dat in ps_list[:]:
				
				dat1 =', '.join(map(str, dat[:]))
				#print dat1
				conserved_plr.write("%s\n" %dat1)
			
conserved_plr.close()
residue_pos = range(1,len(data)+1)
yrange = ylim(-4.5,4.5)

x = [residue_pos]
y = [yrange]		

plt.title('ABC_4BBJ CONSERVED PLRs')
plt.xlabel('SEQUENCE POSITION', fontsize=14)
plt.ylabel('KYTE-RATE SCORES', fontsize=14)

hB, = plot([1],marker='o', label = 'KYTE-SCORE', color = 'green')
hR, = plot([1],marker='o', label = 'RATE-SCORE', color = 'blue')

plt.legend(handler_map={hB: HandlerLine2D(numpoints=2)})
plt.legend(handler_map={hR: HandlerLine2D(numpoints=2)})

hB.set_visible(False)
hR.set_visible(False)

ax.axhline(0, color='red', lw=0.5)

plt.show()
