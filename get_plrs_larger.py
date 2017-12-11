import os
import glob
import re 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import box


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")

print "		step 3.4 : Filtering out the PLRs"

def sorted_nicely(l):
	
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key:[convert(c) for c in re.split('([0-9]+)',key)] 
    return sorted(l, key = alphanum_key)

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
rate_path = "%s/RESULT/MODULE2-LARGER" %curdir_up
ftp_path = "%s/RESULT/MODULE1" %curdir_up
range_paths = []
ftp_fastas = []

for curdir, subdirs, files in os.walk(ftp_path): 
	
		for file in files[:]:	
			
			if '.aln' in file:
				
				name = file.split('.')[0]
				fasta_path = os.path.join(curdir, name + ".fasta")
				ftp_fastas.append(fasta_path)	

#print ftp_fastas			

for curdir, subdirs, files in os.walk(rate_path): 

	for file in files[:]:

		if 'range' in file:
			
			filepath = os.path.join(curdir, file)
			range_paths.append(filepath)			

#print range_paths


source_files = []

for paths in range_paths:
	
	srce = paths.split('/')[:-1]
	srce_path = '/'.join(srce)
	source_files.append(srce_path)	

#print source_files

for i,srce in enumerate(source_files[:]):
	
	source_file = os.chdir(srce)	
	files = sorted_nicely(glob.glob('*.PTin'))
	#print files
	zcors = []
	xcors = []
	ycors=[]
	coins = []
	
	for file in sorted_nicely(files[:]):
		
		#print file

		cord = open(file, 'r').readlines()
		record = False
		for line in cord:
			
			if line.startswith("ATOM"):

				#print f
				x = float(line[31:38])
				y = float(line[39:46])
				z = float(line[47:54])
				plrs = [x,y,z]
				range_file = range_paths[i]
	
				with open(range_file) as box_range:
					
					cords = box_range.read().strip('[').strip(',').strip(']').split(',')
					#print range_file
					#print cords
					
					if z >= float(cords[5]) and z<= float(cords[4]) and y>= float(cords[3]) and y <= float(cords[2]) and x>= float(cords[1]) and x<= float(cords[0]):
					#if z >= float(cords[5]) and z<= float(cords[4]) and x>= float(cords[1]) and x<= float(cords[0]):	
						record = True
						xcors.append(round(x,2))
						ycors.append(round(y,2))
						zcors.append(round(z,2))
		#print xcors
						
		if record:
			
			coins.append(file)
	
		#print sorted_nicely(coins)
		listcontent = []
		listrows = []
		plr = open('plr.txt', 'w')
		
		for lis in coins:
			
			##print lis
			name = lis.split(".")[0]
			new_name = name.split("_")[0]
			plr_list_name = new_name + '.plrs'
			#print plr_list_name
			content =  open(plr_list_name,'r').readlines()
			#print content
			for i,poreplrs in enumerate(content[:]):
				
				if i == 0:
							
					True
							
				else:
							
					if not poreplrs.startswith('AA') and poreplrs != '\n':
						
						listrows.append(poreplrs.split())
						xcords_collumn = poreplrs.split()[3]
						xcords = float(xcords_collumn)
						listcontent.append(round(xcords,2))
		
			#print listcontent
		for xvals in listcontent:
			
			print xvals				
				
			#content = ''.join(content)
			#listcontent += content
		
	#plr.write(listcontent)
	#plr.close()
	
	#plr_list = open('plr_list.txt','w')
	#unipid = srce.split('/')[-2]
	#print unipid
	
	#for fas_path in ftp_fastas[:]:
		
		#if unipid in fas_path:

			#pdb_chain = open(fas_path).read().split(">")
			
			#for pdbids in pdb_chain[:]:
				
				#pds_identifier = pdbids.find('PDBID')
				##print pds_identifier
				
				#if pds_identifier >=1:
					
					#chain = pdbids.split("|")[0]
					#chain_info = chain.split(":")[-1]

			#merge = open('%s/plr.txt'% srce).readlines()
			#final = []
			
			#for info in merge[:]:
				
				#final_lists = list(info.split()[1:])
				#final.append(final_lists)
				
			#for resi in final[:]:
				
				#resi = '\t'.join([str(x) for x in resi])
				#plr_list.write("%s\n" %resi)
						
			#plr_list.close()	

			#ax.scatter([x for x in xcors],[y for y in ycors],[z for z in zcors], alpha = 0.2,color = "r")
	
			#comname = range_paths[i].split('/')[-1].split('_')[0]+"_com .txt"
			#compath = srce +'/' + comname
			#comfiles = open(compath).read().strip('[').strip(']').split(',')
		
			#points = box.box([float(comfiles[0]),float(comfiles[1]),float(comfiles[2])])
			#plt.savefig('filter.png')
		#plt.close("all")
			
	
	

	
