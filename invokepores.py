import os
import subprocess

print "		step 3.3 : Create the coordinate file of the opm files"

def invoke():
	
	bpath = '/home/jacob/Downloads/PROPORES'
	os.environ['PERLLIB'] = bpath + '/module_lib' 
	os.environ['ROTALIB'] =  bpath + '/ROTA_lib_30'
	subprocess.check_call(['perl', bpath + '/Pore_ID.pl' ,'-f', bpath + '/3dhwATOMCB', '-r' ,'1.0','-s','1.2', '-c' ,'1.4', '-n' ,'3dhw'], env = os.environ)
