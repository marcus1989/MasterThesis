import os
import re
import glob
import shutil

print "		step 2.5 : Create the pdb folders with opm files"
   

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
rate_path = "%s/RESULT/MODULE2" %curdir_up


for curdir, subdirs, files in os.walk(rate_path): 

	for file in files[:]:
		
		if '.aln' in file:
				
			newfile = open(os.path.join(curdir,file)).read().split('>')
			
			for ite in newfile[:]:
				
				if "PDBID" in ite:
					
					fold = ite.split('|')[0]
					new_fold = fold.lower()
					new_folds = new_fold.split(":")[0]
					
					if os.path.isdir(os.path.join(curdir,new_folds)):
						
						pass
					else:
						
						os.makedirs(os.path.join(curdir,new_folds))
					
					 
src = os.chdir(curdir_up + "/alpha")
alpha_files = glob.glob("*")

for curdir, subdirs, files in os.walk(rate_path): 

	for subs in subdirs[:]:
		
		pdbs = subs.split(":")[0] + ".pdb"
		
		if pdbs in alpha_files:
			
			src = curdir_up + "/alpha/%s" %pdbs
			dest = curdir + "/%s" %subs
			shutil.copy2(src, dest)

print "End of Module 2"			 
	
