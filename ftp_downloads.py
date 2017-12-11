from Bio import Entrez
from urllib2 import HTTPError
import time
import os
import glob
import ntpath

print "		step 1.2: FTP Download of the NCBI sequences"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])

if os.path.isdir("%s/RESULT/MODULE1/FTP" %curdir_up):
	
	print "		FTP already exists"
	 
else:	
	
	Entrez.email ="ambat.jacob@gmail.com"
	
	path = open('config.txt').read().splitlines()[1].split('=')[-1].strip()

	for root, dirs, files in os.walk(path):

		filepaths = []
		
		for file in files:
			
			filepath = os.path.join(root, file)
			filepaths.append(filepath)
		
		for filepath in filepaths[:]:
			if 'Prot' not in filepath:
				
				p3files = open(filepath,'r').readlines()
				d = {}
				
				for hits in p3files[:]:
					
					key,value = hits.strip("\n").split("\t")[0],hits.strip("\n").split("\t")[1:]
					key = key.strip('\t')
					value = [x.split('|')[0:] for x in value]
					value2 =[x[1] for x in value]
					value3 = [ x for x in value2 if len(x) > 4]

					for items_ids in value3:
						if key not in d:
							d[key] = [items_ids]
						else:
							d[key].append(items_ids)	

			new = ['_ids.txt']
			if 'Prot' not in filepath: 
				
				ftps = ntpath.basename("%s" %filepath)
				newftps = ftps.split('.txt')[:-1]
				filtered_ids_names = ''.join(newftps + new)
				
				with open(os.path.join(path,filtered_ids_names),'w') as output_ids:	
					
					for item in d.keys()[:]:
						
						line = item
						hits = d[item]
						line += '\t%s' %hits
						output_ids.write("%s\n" % line)	
						
					output_ids.close()

				for items in d.keys()[:]:
					
					ftp_path = '%s/RESULT/MODULE1/FTP' %curdir_up
					os.system('mkdir %s' % ftp_path)
					item_path = ftp_path + '/%s' %items
					
					if item_path not in os.listdir(ftp_path):
						
						os.system('mkdir %s' % item_path)
						
					os.chdir(item_path)
					print os.getcwd()
					counter = 0
					values = d[items]
	
					for value in values:
						
						out = open('%s' % value, 'w')
						fetch_success = False
						
						while(fetch_success == False):
							
							try:
								handle = Entrez.efetch(db="protein", id=value, retmode="xml")
								records = Entrez.read(handle)
								fastaseq = value.rstrip()+" "+records[0]["GBSeq_primary-accession"]+" "+records[0]["GBSeq_definition"]+"\n"+records[0]["GBSeq_sequence"]
								out.write(fastaseq)
								out.close()
								fetch_success = True
								
							except:
								
								pass
							time.sleep(1) # to make sure not many requests go per second to ncbi	
							
					counter += 1
				
				

