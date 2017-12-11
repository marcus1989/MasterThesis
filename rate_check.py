import os
import re


print "		step 2.3 : Check if the rate files are empty"
   

curdir = os.getcwd()
curdir_up = '/'.join(curdir.split('/')[:-1])
rate_path = "%s/RESULT/MODULE2" %curdir_up

for curdir, subdirs, files in os.walk(rate_path): 
	
	for file in files[:]:
		
		newpath = os.path.join(curdir, file)
		if os.stat(newpath).st_size == 0:
			os.remove(newpath)
