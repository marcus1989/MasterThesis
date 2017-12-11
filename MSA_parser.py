import os
import re
import fnmatch
from collections import OrderedDict

print "		step 1.4:  Preparing the ClustalO input files"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])

path = open('config.txt').read().splitlines()[2].split('=')[-1].strip()

colfiles = path
lists = os.listdir(colfiles)


seq_input = ["%s/S.AureusCOL" %curdir_up,"%s/S.AureusN315" %curdir_up]

filepaths = []
for root, dirs, files in os.walk(path): 
	  
	for file in files:
		
		if '.fasta' in file:
			
			paths = os.path.join(root, file)
			filepaths.append(paths)



for strains in seq_input[:]:
	
	def seqfilter():
			
		files = open(strains).read().split(">")
		q_ids = []
		qrys = []
		for i,ids in enumerate(files[:]):
			
			query_ids = ids.strip("\n")
			q_ids.append(query_ids)
			query = ids.strip(" ")[0:9]
			querys = query.split(" ")[0]
			
			for filepath in filepaths[:]:
				
				if querys in filepath:
					
					try:
						
						content = open(filepath).read()
					
					
						if ids not in content:
							
							new_content = ">" + ids  + "\n" + content
							newfile = open(path+'/%s.fasta' % querys, 'w')
							newfile.write(new_content)	
					except:
						
						pass
						
				else:
					
					pass
			
	#########################################################################
	
	def blancofilter():		
				
		blanco  = open('%s/blancoDBproteins' %curdir_up).read().split('>')
		blancos = []
		
		for i,pdbs in enumerate(blanco[:]):
			
			pdb_seqs = pdbs.strip("\n")
			blancos.append(pdb_seqs)
			pdbids = pdbs.strip(" ")[0:6]
			
			path3 = '%s/pdbseq' %curdir_up
			
			if not os.path.exists(path3):
				
				os.makedirs(path3)	
				
			with open(os.path.join(path3,"%s.fasta" % pdbids),'w') as outfile:
				outfile.write(pdbs)
	
	#########################################################################
	
	pdbfiles = '%s/pdbseq' %curdir_up
	
	
	def msafileprep():
		
		dirlistings = os.listdir(pdbfiles)
		editfiles =[]
		structure_names_ext = []
		
		for item in dirlistings:
			
			if ".fasta" in item:
				
				editfiles.append(item)

		new_file = ["%s/RESULT/MODULE1/P3/S.AureusCOL_Prot_Blast.txt" %curdir_up,"%s/RESULT/MODULE1/P3/S.AureusN315_Prot_Blast.txt" %curdir_up]
			
		input_fasta_files = []
		
		for input_file in lists:
			
			if input_file.endswith("fasta"):
				
				input_fasta_files.append(input_file)

		fasta_names = []
		
		for fasta_name_ext in input_fasta_files:

			fasta_name = fasta_name_ext.split(".")[0]
			fasta_names.append(fasta_name)

		prot_names = []
		new_fil_pdbs = []
		str_dict = {}
		
		
		for p3out in new_file[:]:
			
			new_files = open(p3out).readlines()
			
			for i,prots in enumerate(new_files[:]):
				
				values = prots.split('\t')
				key = values[0]
				prot_values = [item.split(',')[0].strip() + '.fasta' for item in values[1:]]
				
				if key not in str_dict:
					
					str_dict[key] = prot_values
					
				elif key in str_dict:
					
					str_dict[key].append(prot_values)
		
		set1 = set(editfiles)
		newlist = [key for key, value in str_dict.iteritems() if str_dict.values() and set1]	
		structure_name_ext = [value for key, value in str_dict.iteritems() if str_dict.values() and set1]
		
		
		for strain_name in newlist[:]:	
			
			if strain_name in fasta_names[:]:
				
				fasta_file = open(colfiles + "/" + strain_name+'.fasta').read()
				pdbs = str_dict[strain_name]
				content = fasta_file
				
				for newseqs in pdbs:
					
					sequence = open(pdbfiles + "/" + newseqs).read()
					sequence_str = '>' + sequence
					content += sequence_str
					fasta_files = open(colfiles + "/" + strain_name + '.fasta','w')
					fasta_files.write(content)
					fasta_files.close()
			
			
		for fastas in fasta_names[:]:
			
			filename = colfiles + "/" + fastas + '.fasta'
			
			if fastas not in newlist[:]:
				
				try:
					os.remove(filename)	
					
				except:
					
					pass
				
	seqfilter()	
	blancofilter()
	msafileprep()


