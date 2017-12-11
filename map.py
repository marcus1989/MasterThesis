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
rate_path = "%s/RESULT/MODULE2-LARGER" %curdir_up

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
					
			if 'ATOMCB' in file:
				
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
	protein_name = str.upper(pdb_file.split('/')[-1].strip('ATOMCB'))
	protein_name_lower = pdb_file.split('/')[-1].strip('ATOMCB')

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
		
	print clustal_dict1['4BBJ:A|PDBID|CHAIN|SEQUENCE'][0] ==  clustal_dict1['4BBJ:A|PDBID|CHAIN|SEQUENCE'][1]
	break
