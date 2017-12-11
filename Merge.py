import os
import fnmatch
from collections import Counter
import glob

print "		step 1.3:  Merge the NCBI sequences to one sequence.fasta file"

path = open('config.txt').read().splitlines()[2].split('=')[-1].strip()

dirpaths = []

for root, dirs, files in os.walk(path):  
	
	for dir in dirs:
		
		path = os.path.join(root, dir)
		dirpaths.append(path)
		

for dirpath in dirpaths[:]:
	
	folders = dirpath.split("/")[-1]
	merged_content = ''
	files =  glob.glob('%s/*' % dirpath)
	newcontent = {}
	
	for file in files[:]:
		
		content = open(file).read()	
		key = ">" + content.split()[0]
		value = content.split()[-1]
		
		if key not in newcontent:
			newcontent[key] = value
		elif key in content:
			newcontent[key].append(value)

	new_ncbi_file = dirpath + "/" + "sequence.fasta"

	with open(new_ncbi_file,"w") as fasta_file:
		
		for merge in newcontent.keys()[:]:
			line = merge
			sequences = newcontent[merge]
			line += '\n'
			for merged in sequences:
				
				merged2 = ' '.join(merged)
				line += '%s' % merged2 
			
			fasta_file.write("%s\n" % line)
	fasta_file.close()
	
	
