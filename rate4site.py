import os
import glob

print "		step 2.2 : Rate4Site program"

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
os.chdir(curdir_up + "/RESULT/MODULE2/")

newdir = os.getcwd()

def rate4site():
	
	for curdir, subdirs, files in os.walk(newdir):
	
		for file in files:
			
			rate_file = file.split('.')[0] + '_rate'
			file_path = os.path.join(curdir, file)
			outpath = os.path.join(curdir, rate_file)
			os.system('rate4site -s %s -o %s' %(file_path,outpath))
			
rate4site()	
