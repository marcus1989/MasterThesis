import subprocess
import sys
import os
import glob
from Bio.Align.Applications import ClustalOmegaCommandline

print "		step 1.5:  Multiple Sequence Alignment using ClustalO"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])

os.chdir("%s/RESULT/MODULE1/FTP" %curdir_up)

files = glob.glob("*.fasta")

def clustalo():
	
	for fasta_files in files[:]:
		
		fasta_file_s = fasta_files.split('.')[0]
		print fasta_file_s
		cline = ClustalOmegaCommandline('clustalo',infile = fasta_files,outfile = fasta_file_s + '.aln' ,verbose= False, auto=True)
		cline()
		
clustalo()

print "End of Module 1"
