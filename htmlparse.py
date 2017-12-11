from bs4 import BeautifulSoup
import urllib2
import csv
import os
import sys
import glob
import shutil

print "Module 2"

print "		step 2.1 : Creating the transporter folders and alignment files"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])

os.chdir("%s/RESULT/MODULE1/FTP" %curdir_up)
aln_files = glob.glob("*.aln")

strains = {}

for i in range(3):
	
	url=('http://www.membranetransport.org/all_type_btab.php?oOID=saur%s' % i)
	header = {'User-Agent': 'Morzilla/5.0'}
	req=urllib2.Request(url, headers = header)
	page = urllib2.urlopen(req)
	soup = BeautifulSoup(page)
	
	
	ORF = ""
	Family_ID = ""
	Family_TC = ""
	Family_Name = ""
	Transporter_Type = ""
	Substrate = ""
	Note = ""
	TC = ""
	
	table = soup.find('body', {"bgcolor" : "#FFFFFF"})
	

	for row in table.findAll("tr"):
		
		cells = row.findAll("td")
		
		if len(cells) == 8:
			
			ORF = cells[0].find(text = True)
			ORF = str(ORF)
			
			if 'ORF' not in ORF:
				
				Family_ID = cells[1].find(text = True)
				Family_ID = str(Family_ID)
				Family_TC = cells[2].find(text = True)
				Family_Name = cells[3].find(text = True)
				Transporter_Type = cells[4].find(text = True)
				Substrate = cells[5].find(text = True)
				Note = cells[6].find(text = True)
				TC = cells[7].find(text = True)
				
				if Family_ID not in strains:
					
					strains[Family_ID] = [ORF]
					
				elif Family_ID in strains:
					
					strains[Family_ID].append(ORF)
				 

file_names = []

for alns in aln_files[:]:
	names  =  alns.split('.aln')[0]
	file_names.append(names)
	

f = open('%s/RESULT/transport_families.txt' %curdir_up, 'wb')

if os.path.isdir("%s/RESULT/MODULE2" %curdir_up):
	
	print "	folder already exists"
	 
else:	
	
	for item in strains.keys()[:]:
		
		line = item
		folders = os.makedirs(curdir_up + "/RESULT/MODULE2/%s" %line)
		hits = strains[item]
		
		for news in file_names:
			
			if news in hits:
				
				news_folders = os.makedirs(curdir_up + "/RESULT/MODULE2/%s/%s" %(line,news))
				src = curdir_up + "/RESULT/MODULE1/FTP/%s.aln" %news
				src_1 = curdir_up + "/RESULT/MODULE1/FTP/%s.fasta" %news
				dest = curdir_up + "/RESULT/MODULE2/%s/%s" %(line,news)
				shutil.copy2(src, dest)
				
		for hit in hits:
			
			line += '\n\t%s\n' % hit
			
		f.write("%s\n" % line)	 


root_dir = "%s/RESULT/MODULE2" %curdir_up
newaln = glob.glob('*.aln')

empty_count = 0

for curdir, subdirs, files in os.walk(root_dir):
	
	if len(subdirs) == 0 and len(files) == 0: 
		
		empty_count += 1 
		os.rmdir(curdir)
		
		
