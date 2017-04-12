#! /usr/bin/python
import os
import sys
import string
import signal
import fnmatch
import pandas as pd

from datetime import datetime
from subprocess import Popen
from subprocess import PIPE
from operator import itemgetter
from itertools import groupby

size_min_num = 17

bevelPath = "/share/biocore/internal_projects/seqmatch/bevel/bin/bevel"


def createDBListOrig (file, rd): #create list of db files
	outfile = open(file+'.txt', 'w')

	dblist = []
	for root, dirnames, filenames in os.walk(rd):
			for filename in fnmatch.filter(filenames, '*.fna'):
				dblist.append(os.path.join(root, filename))

	for item in dblist:
		outfile.write("%s\n" % item)

	outfile.close()

def createDBList (file): #create list of db files from text file
	dblist = open(file+".txt",'r').readlines()
	dblist = [x.strip() for x in dblist]

	return dblist

def wrapBev(bevelPath, targetDB, queryDB, writeDB = False, nMinimizer = 100, sizeMinimizer = size_min_num):

	if size_min_num >= 32:
			print 'Minimizer size must be less than or equal to 32\n'
			exit()

	call = bevelPath

	if writeDB:
		call = call + ' -d'

	call = call + ' -w ' + str(nMinimizer) + ' -k ' + str(sizeMinimizer)

	call = call + ' ' + targetDB + ' ' + queryDB

	with open(os.devnull, 'w') as DEVNULL:
		p = Popen(call,
				stdout=PIPE,
				stderr = DEVNULL,
				bufsize=-1,
				shell=True,
				executable='/bin/bash',
				preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
	if p.returncode:
		raise
	return p.stdout

def analyzeTarget(targetDBs, queryDBs, oOutput):
	# try:

		total=len(targetDBs)*len(queryDBs) #how many targets processed
		# target: 149385.12 ==> query: 41-3_S39.
		#targetDB = '/share/biocore/internal_projects/seqmatch/genomes/149385.12/149385.12.fna'
		#queryDB = '/share/biocore/internal_projects/seqmatch/03-SpadesAssemblies/41-3_S39/41-3_S39.AllContigs.fna'

		count = 0
		for queryDB in queryDBs:
			for targetDB in targetDBs:
				count += 1
				# progress=count/float(total)
				update_progress(count/float(total)) #status bar

				# get target DB name
				targetDBName = os.path.splitext(os.path.basename(targetDB))[0]

				# get query DB name
				queryDBName = os.path.splitext(os.path.basename(queryDB))[0]

				lOutput=[]

				for line in wrapBev(bevelPath, targetDB, queryDB):

			#query Seqid	target Seqid	query Start	target Start	# minimizers found in target	# minimizers found in query
#use dataframe
					result = {} #create a dictionary
					line2 = line.strip().split() # split line into separate items

					result['qDB'] = queryDBName #take each line and put into a dictionary and add dictionary to list
					result['qseqID'] = line2[0]
					result['tDB'] = targetDBName
					result['tseqID'] = line2[1].replace('accn|','')
					result['tminz'] = int(line2[4])
					result['qminz'] = int(line2[5])

					lOutput.append(result.copy()) #create list of dictionaries

				if lOutput:
					sortedlist = []

					sortedlist = list(lOutput)

					sortedlist.sort(key=itemgetter('qseqID','tseqID'))

					df = pd.DataFrame() #create empty dataframe

					# Create dataframe from list of dictionaries
					df = pd.DataFrame(sortedlist)
					# print ('create DF\n')

					# Sum query minimizers for lines with same qseqID and same tseqID
					df = df.groupby(['qseqID', 'tseqID','qDB','tDB'])['qminz'].sum().reset_index()
					# df = df.groupby(['qseqID', 'tseqID','qDB','tDB'], as_index=False)['qminz'].sum()
					# print ('sum query minimizer\n')

					# Sort dataframe by qminz in descending order
					df = df.sort(['qminz'], ascending=[False])
					# print ('sort DF\n')

					df = df[['qDB','qseqID','tDB','tseqID','qminz']]

					# Write results to output file
					df.to_string(oOutput,index=False,header=False)
					oOutput.write('\n')

				else:
					oOutput.write('No Minimizers -- ' + queryDBName + ' || ' + targetDBName + '\n')
					continue

	# except:
	# 	sys.stderr.write('There is a problem! (analyzeTarget)\n')

def update_progress(progress):
	barLength = 20 # Modify this to change the length of the progress bar
	status = ""
	if isinstance(progress, int):
		progress = float(progress)
	if not isinstance(progress, float):
		progress = 0
		status = "error: progress var must be float\r\n"
	if progress < 0:
		progress = 0
		status = "Halt...\r\n"
	if progress >= 1:
		progress = 1
		status = "Done...\r\n"
	block = int(round(barLength*progress))
	text = "\rPercent: [{0}] {1}% {2}".format( "="*block + " "*(barLength-block), format(progress*100,'.2f'), status)
	sys.stdout.write(text)
	sys.stdout.flush()

def main():
	# For testing
	# targetDB = '/share/biocore/internal_projects/seqmatch/genomes/149385.12/149385.12.fna'
	# queryDB = '/share/biocore/internal_projects/seqmatch/03-SpadesAssemblies/41-3_S39/41-3_S39.AllContigs.fna'
#Create list from text files - 60000targetdbs 291 querydbs , 5000 target dbs
	# print os.path.splitext(os.path.basename(sys.argv[0]))[0]
	# Create list of target DBs
	targetDBs = createDBList('target-1000-3')
	# print len(targetDBs)

	# Create list of query DBs
	queryDBs = createDBList('query-files-all')
	# print len(queryDBs)

	# Create output filename which includes time stamp
	outputFile = 'output_' + datetime.now().strftime("%Y%m%d-%H%M%S") + '_' + os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.txt'

	# Create/open output file object
	oFile = open(outputFile, "w")

	analyzeTarget(targetDBs, queryDBs, oFile)

	# Done writing to output file so close it
	oFile.close()

if __name__ == '__main__':