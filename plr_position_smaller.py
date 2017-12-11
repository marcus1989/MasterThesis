import os
import re
from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt 
from pylab import *
from weblogolib import *
import pylab as plt
from matplotlib.legend_handler import HandlerLine2D

print "		step 4 : Mapping the PLRs and Conserved residues on the primary sequence of protein"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
rate_path = "%s/RESULT/MODULE2-SMALLER" %curdir_up


plr_destinations = []
clustal_files = []
pdb_files = []
fasta_files = []
conserved = []

for curdir, subdirs, files in os.walk(rate_path): 
	
		for file in files[:]:	
			
			if 'plrs.txt' in  file:
				
				plr_destination = os.path.join(curdir,file)
				plr_destinations.append(plr_destination)

			if '.aln' in file:

				clustal_file = os.path.join(curdir,file)
				clustal_files.append(clustal_file)
					
			if 'ATOM' in file:
				
				pdb_file = os.path.join(curdir,file)
				pdb_files.append(pdb_file)
				
			if '.fasta' in file:
				
				fasta_file = os.path.join(curdir,file)
				fasta_files.append(fasta_file)
			
			if '_conserved_residues.txt' in file:
			
				conserved_file = os.path.join(curdir, file)
				conserved.append(conserved_file)
		
for ipdb,pdb_file in enumerate(pdb_files):
	
	curdir2 = '/'.join(pdb_file.split('/')[:-1])
	
	fasta_file = ''
	clustal_file = ''
	conserved_file = ''
	
	plr_destination = plr_destinations[ipdb]
	protein_name = str.upper(pdb_file.split('/')[-1].strip('ATOM'))
	protein_name_lower = pdb_file.split('/')[-1].strip('ATOM')

	strain_name = pdb_file.split('/')[-3]
	
	for fastafile in fasta_files:
		
		if strain_name in fastafile:
			
			fasta_file = fastafile
			
	
	for clustalfile in clustal_files:
		
		if strain_name in clustalfile:
			
			clustal_file = clustalfile
				
	for cons_set in conserved:
		
		if strain_name in cons_set:
			
			conserved_file = cons_set
	
							
	clustals = open(clustal_file).read().split('>')
	pdbs = open(pdb_file).readlines()
	fastas = open(fasta_file).read().split('>')	
	data_conserved_file = open(conserved_file).read()
	
	clustal_dict1 = {}
	
	for qids in clustals[:]:

		if protein_name in qids :
			
			strs = qids.split()
			str_ids = strs[0]
			sequences = ''.join(strs[1:]).split()

			if str_ids not in clustal_dict1:
				
				clustal_dict1[str_ids] = sequences
				
			elif str_ids in clustal_dict1:
				
				clustal_dict1[str_ids].append(sequences)
		else:
			pass
		
	##print clustal_dict1 
	##break
	
	clustal_dict2 = {}
	
	for seqs in clustal_dict1:
		
		sequences = clustal_dict1[seqs]
		for residues in sequences:
			
			ind_resis = list(residues)
			
			pairs = []
			for iinres, ind_resi in enumerate(ind_resis):
				pairs.append([iinres+1, ind_resi])
			clustal_dict2[seqs] = pairs
	
	#print clustal_dict2, clustal_dict1
	#break
	
	pdb_dict1 = {}
	
	for maps in pdbs:
		
		residue_name = maps[17:20]
		residue_sequence_number = int(maps[23:26])
	
		if residue_sequence_number not in pdb_dict1:
			
			 pdb_dict1[residue_sequence_number] = [residue_name]
			 
		elif residue_sequence_number in pdb_dict1:
			
			pdb_dict1[residue_sequence_number].append(residue_name)
	
	#print pdb_dict1,ipdb
	#break
		
	amino_name = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]
	code = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"] 
	
	amino_dict = {}
	
	for i in range(len(amino_name)):
		
		amino_dict[amino_name[i]]=code[i]
	pdb_dict2 = {}
	
	for pos in pdb_dict1:
		
		res = pdb_dict1[pos]
		
		for resi in res:
			
			single = amino_dict[resi]
			
			if pos not in pdb_dict2:
				
				pdb_dict2[pos] = [single]
				
			elif pos in pdb_dict2:
				
				pdb_dict2[pos].append(single)	
	
	#print len(pdb_dict2),len(pdb_dict1)
	#break
	
	pdb_dict3 = {}	
	
	for single_items in pdb_dict2:
		
		pdb_dict3[single_items] = pdb_dict2[single_items][0]
	
	#print len(pdb_dict3)
	#break
	
	start_pos = min(pdb_dict3.keys())
	
	start_res = pdb_dict3[start_pos]
	
	#print start_pos, start_res
	#break
	
	new_clustal_dict1 ={}	
			
	for keys in clustal_dict2:
		
		values = clustal_dict2[keys]
		new_clustal_dict1[keys] = []
		start = False
		start_pos2 = start_pos
		
		for items in values:
			
			if items[1] == '-':
				
				pass
				
			else:	
				
				res =  items[1]
		
				if res == start_res and not start:
					
					new_pair  = start_pos2, res
					new_clustal_dict1[keys].append(new_pair)	
					start = True
					start_pos2 += 1
					
				elif start:
					
					new_pair  = start_pos2, res
					new_clustal_dict1[keys].append(new_pair)
					start_pos2 += 1
	
			
	for clus in new_clustal_dict1.keys()[:]:
		
		line = clus+ '\n'
		hits = new_clustal_dict1[clus]
		clustal_dict_file = open(curdir2 + '/%sclustal_file' % clus, 'w')
		content = ''
		
		for hit in hits:
			
			newcontents = str(hit[0]) + ',' + str(hit[1]) + '\n'
			content += newcontents
		clustal_dict_file.write(content)
		
		clustal_dict_file.close()					
			
	#print new_clustal_dict1,clustal_dict2
	#break
	
	new_clustal_dict2={}	
	plr_file = open(plr_destination).readlines()
	plrs_tuple = []
	
	for i,plr in enumerate(plr_file[:]):
		
		residue =  plr.split()[0]
		plr_s =  re.findall(r"'\s*([^']*?)\s*'",plr[1:-1])
		newplrs = [int(x) for x in plr_s]
		
		for i2,items in enumerate(newplrs[:]):
			
			plr = tuple((newplrs[i2],residue))
			plrs_tuple.append(plr)	
			
	vs = new_clustal_dict1.values()
	ks = new_clustal_dict1.keys()
	
	for j,k in enumerate(ks):
		
		new_clustal_dict2[k] = []
		for n in vs[j]:
			
			if n in plrs_tuple:
				
				n2 = (n[0],'p')
				new_clustal_dict2[k].append(n2)

			else:
				
				n = (n[0],'*')
				new_clustal_dict2[k].append(n)
			
	#print new_clustal_dict2, new_clustal_dict1
	#break
	
	fastafile1 = {}
	
	for s,seq in enumerate(fastas[:]):
		
		if protein_name in seq:
			
			fasta_strs = seq.split()
			fasta_ids = fasta_strs[0]
			fasta_seqs = ''.join(fasta_strs[1:]).split()
	
			if fasta_ids not in fastafile1:
				
				fastafile1[fasta_ids] = fasta_seqs
				
			elif fasta_ids in fastafile1:
				
				fastafile1[fasta_ids].append(fasta_seqs)
			
	#print fastafile1
	#break
	
	fastafile2 = {}
	
	for seqs in fastafile1:
		
		sequences = fastafile1[seqs]
		
		for residues in sequences:
			
			ind_resis2 = list(residues)
			pairs = []
	
			for ires, ind_residues in enumerate(ind_resis2):
				
				pairs.append([ires+1, ind_residues])
			fastafile2[seqs] = pairs
			
	#print fastafile1, fastafile2
	#break
	
	fastafile3 = {}
	pre = ''
	post = ''
	for items in new_clustal_dict2.keys():
		
		line = items + '\n'
	
		min_pdb_len = min(new_clustal_dict2[items])
		pos1,res1 = tuple(min_pdb_len)
		max_pdb_len = max(new_clustal_dict2[items])
		pos2,res2 = tuple(max_pdb_len)
		
		#print min_pdb_len,max_pdb_len	
	
		for item in fastafile2.keys():
			
			newline = item + '\n'
			min_newseqdat_len = min(fastafile2[item])
			pos3,res3 = tuple(min_newseqdat_len)
			max_newseqdat_len = max(fastafile2[item])
			pos4,res4 = tuple(max_newseqdat_len)
			
			#print min_newseqdat_len, max_newseqdat_len	
			
			min_len = pos1 - pos3
			max_len = pos4 - pos2
			
			#print min_len, max_len
			
			for plr in new_clustal_dict2.keys():
					
				pre = range(1, min_len + 1)
				post = range(1, max_len +1)
			
		for num in pre:
			
			positions = min_pdb_len[0] - num
			new_pos = (positions, '*')
			new_clustal_dict2[items].insert(0,new_pos)
		
		for num1 in post:
			
			positions1 = max_pdb_len[0] + num1
			new_pos1 = (positions1, '*')
			new_clustal_dict2[items].append(new_pos1)  
		fastafile3[items] = new_clustal_dict2[items]	
	
	#print fastafile3, fastafile2
	#print new_clustal_dict2

	data1 = data_conserved_file.split(')')
	data2 = [ item.replace('(','').strip().split(',') for item in data1[:-1]]
	data3 = [(int(item[0]), item[1].strip().strip("'")) for item in data2]
	
	res = new_clustal_dict1.values()
	cons = new_clustal_dict1.keys()
	conserved_set = {}

	for c, con in enumerate(cons):
		
		conserved_set[con] = []
		
		for s in new_clustal_dict1[con]:

			if s in data3:
				
				s2 = (s[0],'c')
				conserved_set[con].append(s2)
				
			else:
			
				s = (s[0],'*')
				conserved_set[con].append(s)
	
	#print conserved_set
	#print new_clustal_dict1
	#print fastafile3

	new_conserved_set= {}
	pre1 = ''
	post1 = ''
	
	for items in conserved_set.keys():
		
		line = items + '\n'
		min_pdb_len = min(conserved_set[items])
		pos1,res1 = tuple(min_pdb_len)
		max_pdb_len = max(conserved_set[items])
		pos2,res2 = tuple(max_pdb_len)
		
		#print min_pdb_len,max_pdb_len	
	
		for item in fastafile2.keys():
			
			newline = item + '\n'
			min_newseqdat_len = min(fastafile2[item])
			pos3,res3 = tuple(min_newseqdat_len)
			max_newseqdat_len = max(fastafile2[item])
			pos4,res4 = tuple(max_newseqdat_len)
			
			#print min_newseqdat_len, max_newseqdat_len	
			
			min_len = pos1 - pos3
			max_len = pos4 - pos2
				
			for plr in conserved_set.keys():
					
				pre1 = range(1, min_len + 1)
				post1 = range(1, max_len +1)
				
			for num in pre1:
				
				positions = min_pdb_len[0] - num
				new_pos = (positions, '*')
				conserved_set[items].insert(0,new_pos)
			
			for num1 in post1:
				
				positions1 = max_pdb_len[0] + num1
				new_pos1 = (positions1, '*')
				conserved_set[items].append(new_pos1) 
				 
			new_conserved_set[items] = conserved_set[items]
	
	
	#print fastafile3
	#print fastafile3.keys()
	#print new_conserved_set.keys()
	#print pdb_file
	#print new_conserved_set
	
	for fast in fastafile3.keys()[:]:

		fastahits = fastafile3[fast]
		fasta_dict_file = open(curdir2 + '/%sfastafile3' % fast, 'w')
		fasta_content = ''
		
		for hit in fastahits:
			
			fastcontents = str(hit[0]) + ',' + str(hit[1]) + '\n'
			fasta_content += fastcontents
		fasta_dict_file.write(fasta_content)
		
		fasta_dict_file.close()				
	
	for consd in new_conserved_set.keys()[:]:
		
		conservedhits = new_conserved_set[consd]
		conserved_dict_file = open(curdir2 + '/%sconserved_set_file' % consd, 'w')
		conserved_content = ''
		
		for hit in conservedhits:
			
			conscontents = str(hit[0]) + ',' + str(hit[1]) + '\n'
			conserved_content += conscontents
		conserved_dict_file.write(conserved_content)
		
		conserved_dict_file.close()	
	
	plr_pos = open(curdir2 + '/plr_pos', 'w')
	
	#print fastafile3, new_conserved_set
	#print pdb_file
	for item in fastafile3.keys()[:]:
		
		line = item + '\n'
		hits = fastafile3[item]	
		residues = []
		residues1 = []
		residues2 = []
		cosposs = []
		poss = []
		for cos in new_conserved_set.keys()[:]:
			
			for costuples in conserved_set[cos]:
				
				cospos,cosres = costuples
				
				if cos in item:
					
					residues.append(cosres) 
					cosposs.append(cospos)
		##print residues	
		
		for tuples in hits:
			
			pos,residue_name = tuples
			residues1.append(residue_name)
			poss.append(pos)
				
		for res in fastafile1.keys()[:]:	
			
			if res in item:
				
				residues2.append(''.join(fastafile1[res]))
				
		#print len(residues) #,cosposs
		#print 'yo'
		#print len(residues1) #, poss
		#print 'yoho'
		#print len(residues2[0]) #, len(residues2[0])
		#print 'yoyoyo' ,pdb_file
	
		
		for hit in residues:		
			
			line += '%s' %hit
		
		line += '\n'
	
		for hits in residues1:		
			
			line += '%s' %hits
		line += '\n'
		
		for hitss in residues2:		
			
			line += '%s' %hitss
			
		plr_pos.write("%s\n" % line)	
	
	plr_pos.close()
	
print "End of Module 4"	
	
	
