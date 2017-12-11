import os
import glob

print "		step 3.5 : Filtering out the PLRs"

#curdir = os.getcwd()
#curdir_up = '/'.join(curdir.split('/')[:-1])
#rate_path = "%s/RESULT/MODULE2-LARGER" %curdir_up

rate_path = '/home/jacob/Downloads/PROPORES/'

for curdir, subdirs, files in os.walk(rate_path): 
	
		for file in files[:]:	
			
			if 'plr_list.txt' in file:
				
				amino_name = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]
				code = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"] 
				
				amino_dict = {}
				
				destination = os.path.join(curdir,file)
					
				for i in range(len(amino_name)):
					amino_dict[amino_name[i]]=code[i]
				
				os.chdir('/'.join(destination.split('/')[:-1]))
				
				plr_files = set(glob.glob('plr_list.txt'))
				
				dictl = {}
				positions = []
				for line in plr_files:
					hold = True
					contents = open(line).readlines()
				
					for  content in contents[:]:
						con = content.split()
						key = con[1]
				 		pos = content.split()[0]
						positions.append(pos)
					
						if key not in dictl:
							dictl[key] = [pos]
						elif key in dictl:
							dictl[key].append(pos)
				
				newdict = {}
				
				for key, val in dictl.iteritems():
				
					new_key = amino_dict[key]
					new_val = val
					
					if new_key not in newdict:
						newdict[new_key] = [new_val]
					elif new_key in newdict:
						newdict[new_key].append(new_val)
						
					 	 
				plrs = open('plrs.txt','w')
				for items in newdict.keys():
					lines = items
					hit = newdict[items]
					for hits in hit:
						lines +='\t%s' % hits
					plrs.write("%s\n" %lines)

print "End of Module 3"	
