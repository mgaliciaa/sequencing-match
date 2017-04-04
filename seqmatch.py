#! /usr/bin/python

import os
import sys
import string
import signal
import fnmatch 

from datetime import datetime
from subprocess import Popen
from subprocess import PIPE

size_min_num = 17

#targetDB = '/share/biocore/internal_projects/seqmatch/genomes/59814.6/59814.6.fna'
#queryDB = '/share/biocore/internal_projects/seqmatch/03-SpadesAssemblies/26-3_S36/26-3_S36.Scaffolds.fna'
#if not os.path.isfile(queryDB):
#	print 'Query file does not exist\n'
#	exit()
#print 'Query file exists\n'

#create list of fna files
#rootdir ='/share/biocore/internal_projects/seqmatch/genomes'
#fnas = []

bevelPath = "/share/biocore/internal_projects/seqmatch/bevel/bin/bevel"


def createDBList (rd): #create list of db files
	dblist = []
	for root, dirnames, filenames in os.walk(rd):
    		for filename in fnmatch.filter(filenames, '*.fna'):
        		dblist.append(os.path.join(root, filename))

	if not dblist:
		raise
	return dblist

def wrapBev(bevelPath, targetDB, queryDB, writeDB = False, nMinimizer = 100, sizeMinimizer = size_min_num):
	
	if size_min_num >= 32:
        	print 'Minimizers must be less than or equal to 32\n'
        	exit()

	call = bevelPath

	if writeDB:
		call = call + ' -d'

	call = call + ' -w ' + str(nMinimizer) + ' -k ' + str(sizeMinimizer)

	call = call + ' ' + targetDB + ' ' + queryDB

	p = Popen(call,
				stdout=PIPE,
				stderr = None,
				bufsize=-1,
				shell=True,
				executable='/bin/bash',
				preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
	if p.returncode:
		raise
	return p.stdout 

def analyzeTarget(targetDB, queryDB, oOutput):
	try:

		i = 0
		

		# get target DB name
		targetDBName = os.path.splitext(os.path.basename(targetDB))[0]

		# get query DB name
		queryDBName = os.path.splitext(os.path.basename(queryDB))[0]

		oOutput.write('\ntarget: ' + targetDBName + ' ==> ' + 'query: ' + queryDBName + '\n')
	
		for line in wrapBev(bevelPath, targetDB, queryDB):
			i += 1
			oOutput.write(line)
	
		print ("Finished processing file, read ", i ,  " lines\n")

	except:
		sys.stderr.write('There is a problem!\n')
	
def main():

	# Create list of target DBs
	targetDBs = createDBList('/share/biocore/internal_projects/seqmatch/genomes')

	# Create list of query DBs
	queryDBs = createDBList('/share/biocore/internal_projects/seqmatch/03-SpadesAssemblies/')

	# Create output filename which includes time stamp
	outputFile = 'output_' + datetime.now().strftime("%Y%m%d-%H%M%S") + '.txt'

	# Create/open output file object
	oFile = open(outputFile, "w")

	for queryDB in queryDBs[:5]:
		for targetDB in targetDBs[:10]:
			analyzeTarget(targetDB, queryDB, oFile)

	# Done writing to output file so close it
	oFile.close()

if __name__ == '__main__':
	main()
