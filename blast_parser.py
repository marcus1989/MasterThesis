###########Seq_Profile_Blast_Parser_tool################################
import csv
import time
import re
import os
import sys
from collections import Counter
import operator
from fractions import *
import glob
import ntpath
from collections import defaultdict



path = open('config.txt').read().splitlines()[0].split('=')[-1]

rootDir = '.'		
blast_files = []
curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])


for dirName, subdirList, fileList in os.walk(rootDir, topdown = False):
	
	for fname in fileList:
		
		if fname.startswith("S.A"):
			
			fname = os.path.join(dirName, fname)
			blast_files.append(fname)

print 'Module1'
print '		step 1.1 : Parsing the input Blastp files'	


for blastfiles in blast_files[:]:

	if 'Prot' not in blastfiles:
		
		qids=[]
		query_lengths = []
		counter = 0
		seqsas = []
		file1 = open(blastfiles,'r').read()
		queries = file1.split('Query=')
		datas = queries[1:]	
		
		for item in datas[:]:
			
				lines = item.split('\n')
				qid = item.split()[0]
				qids.append(qid)
				
				for line in lines[:]:
					
					if line.startswith('Length='):
						
						query_lengths.append(int(line.split('=')[-1]))
						break
					
		for i,data in enumerate(datas[:]):
			
			lines = data.split('\n')
			record = False
			
			for line in lines[:]:
				
				if line.startswith(">") :
					
					tmp = line.split(">")
					tmp_name = tmp[1]
					tmp_name1 = tmp_name.split("[")	
					tmp_hit = ''.join(tmp_name[0:-1])
					
					if 'Staphylococcus' in line:
						
						record  = True
						
					else:
						
						record = False
						
				if line.startswith(" Score") and record:
					
					tmp = line.strip().split()
					tmp_score_s = tmp[2]
					tmp_score = float(tmp_score_s)
					tmp_evalue = float(tmp[7].replace(",",""))
					seqsas.append([qids[i],tmp_hit,tmp_score,tmp_evalue])
		
				if line.startswith(" Identities")and counter <len(seqsas) and record:
					
					tmp = line.strip().split() 
					tmp_id = tmp[3]
					tmp_ids = tmp_id.replace('(','').replace(')','').replace('%','').replace(',','')
					ids = int(tmp_ids)
					tmp_po = tmp[7]
					tmp_pos = tmp_po.replace('(','').replace(')','').replace('%','').replace(',','')
					pos = int(tmp_pos)
					tmp_gap = tmp[11]
					tmp_gaps = tmp_gap.replace('(','').replace(')','').replace('%','').replace(',','')
					gaps_percent = int(tmp_gaps)
					gap_number = int(tmp[10].split('/')[0])
					alignment_length = int(tmp[10].split('/')[-1])
					coverage_percent = round(float((alignment_length - gap_number))/query_lengths[i] * 100, 2)
					seqsas[counter].append(ids)
					seqsas[counter].append(pos)
					seqsas[counter].append(gaps_percent)
					seqsas[counter].append(gap_number)
					seqsas[counter].append(alignment_length)
					seqsas[counter].append(coverage_percent)
					counter+=1
					
		path1 = '%s/RESULT/MODULE1/P1' % curdir_up

		if not os.path.exists(path1):
			
			os.makedirs(path1)
		file_name = ntpath.basename('blast_out1%s' % blastfiles) + '.txt'
		
		with open(os.path.join(path1,file_name),'w') as out1:
			
			for item in seqsas[:]:
				
				item = '\t'.join([str(x) for x in item])
				out1.write('%s\n' %item)
				
			out1.close()
					
	else:
		
		strsas = []
		qids=[]
		query_lengths = []
		counter = 0
		file2 = open(blastfiles,'r').read()
		queries = file2.split('Query=')
		datas = queries[1:]	
		
		for item in datas[:]:
			
				lines = item.split('\n')
				qid = item.split()[0]
				qids.append(qid)
				
				for line in lines[:]:
					
					if line.startswith('Length='):
						
						query_lengths.append(int(line.split('=')[-1]))
						break
						
		for i,data in enumerate(datas[:]):
			
			lines = data.split('\n')
			record = False
			
			for line in lines[:]:
				
				if line.startswith(">") :
					
					tmp = line.split(">")
					tmp_name = tmp[1]
					tmp_hit = tmp_name.split("|")[0]
			
					
				if line.startswith(" Score") :
					
					tmp = line.strip().split()
					tmp_score_s = tmp[2]
					tmp_score = float(tmp_score_s)
					tmp_evalue = float(tmp[7].replace(",",""))
					strsas.append([qids[i],tmp_hit,tmp_score,tmp_evalue])
		
				if line.startswith(" Identities") and counter < len(strsas):
					
					tmp = line.strip().split()
					tmp_id = tmp[3]
					tmp_ids = tmp_id.replace('(','').replace(')','').replace('%','').replace(',','')
					ids = int(tmp_ids)
					tmp_po = tmp[7]
					tmp_pos = tmp_po.replace('(','').replace(')','').replace('%','').replace(',','')
					pos = int(tmp_pos)
					tmp_gap = tmp[11]
					tmp_gaps = tmp_gap.replace('(','').replace(')','').replace('%','').replace(',','')
					gaps_percent = int(tmp_gaps)
					gap_number_1 = Fraction(tmp[10])
					gap_number = int(tmp[10].split('/')[0])
					alignment_length = int(tmp[10].split('/')[-1])
					coverage_percent = round(float((alignment_length - gap_number))/query_lengths[i] * 100, 2)
					strsas[counter].append(ids)
					strsas[counter].append(pos)
					strsas[counter].append(gaps_percent)
					strsas[counter].append(gap_number)
					strsas[counter].append(alignment_length)
					strsas[counter].append(coverage_percent)
					counter +=1
					
		path1 = '%s/RESULT/MODULE1/P1' %curdir_up
		
		if not os.path.exists(path1):
			
			os.makedirs(path1)
		prot_file_name = ntpath.basename('prot_blast_out1%s' % blastfiles) + '.txt'
		
		with open(os.path.join(path1,prot_file_name),'w') as out2:
			
			for item in strsas[:]:
				
				item = '\t'.join([str(x) for x in item])
				out2.write('%s\n' %item)
				
			out2.close()
		
def parser2():
	
		os.chdir('%s/RESULT/MODULE1/P1' %curdir_up)
		
		for file1 in glob.glob('*.txt'):
			file_s = open(file1).readlines()
			prepsas = []
			
			for item in file_s[:]:
				
				item = item.strip().split('\t')
				hit = item[1]
				e = float(item[3])
				ids = int(item[4])
				cov = float(item[9])
				if e <=1e-10 and ids >= 35 and cov >= 75:
					
					prepsas.append(item)
						
				if len(item) < 10:
					
					print 'not match'
			
			prot_file_name_s = str(file1) 
			
			path2 = '%s/RESULT/MODULE1/P2' %curdir_up
			
			if not os.path.exists(path2):
				
				os.makedirs(path2)
		
			with open(os.path.join(path2,prot_file_name_s),'w') as prepsas1:
				
				for hits in prepsas[:]:
					
					hits = '\t'.join([str(x) for x in hits])
					prepsas1.write('%s\n' %hits)
					
				prepsas1.close()
		
				
def parser3():
	
	os.chdir('%s/RESULT/MODULE1/P2' %curdir_up)
	
	for file2 in glob.glob('*.txt'):
		
		file3 =open(file2).readlines()
		d = {}
		
		for filters in file3[:]:
			
			key, value = filters.strip("\n").split("\t")[0],filters.strip("\n").split("\t")[1:]
			key = key.strip('\t')
			value = [str(x)[0:]for x in value]
			
			if key not in d:
				
				d[key] = [value]
				
			elif key in d and len(d[key]) <= 250:
				
				d[key].append(value)
				
		prot_file_name_s = str(file2) 
		
		path2 = '%s/RESULT/MODULE1/P3' %curdir_up
		
		if not os.path.exists(path2):
			
			os.makedirs(path2)	
				
		with open(os.path.join(path2,prot_file_name_s),'w') as fp:
			
			for item in d.keys()[:]:
				
				line = item
				hits = d[item]

				for hit in hits:
					
					hit2 = ','.join(hit)
					line += '\t%s' % hit2
					
				fp.write("%s\n" % line)

	
parser2()
parser3()
	
	
	
	

