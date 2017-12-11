import os

print "		step 3.1 : Create the coordinate file of the opm files"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
rate_path = "%s/RESULT/MODULE2-LARGER" %curdir_up

for curdir, subdirs, files in os.walk(rate_path): 

	for file in files[:]:
		
		if '.pdb' in file:
			
			atom_data = open(os.path.join(curdir, file),'r').readlines()
			outfile = open(os.path.join(curdir, file.split(".")[0]) + "ATOM", 'w')
			for atoms in atom_data:
				if atoms.startswith('ATOM'):
					outfile.write(atoms)
			outfile.close()
			
for curdir, subdirs, files in os.walk(rate_path): 

	for file in files[:]:
		
		if 'ATOM' in file:
			file = os.path.join(curdir, file)
			os.system('grep CB %s > %sCB' %(file,file))
			print file
print 'done'
			
	
