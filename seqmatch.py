#! /usr/bin/python

import os
import sys
import getopt
import string
import signal
import fnmatch
import pandas as pd

from datetime import datetime
from subprocess import Popen
from subprocess import PIPE
from operator import itemgetter
from itertools import groupby

nKmer = 15
nFilter = 10
nMin = 100
# nMin = 2

bevelPath = "/share/biocore/internal_projects/seqmatch/bevel/bin/bevel"

def createDBList (file): #create list of db files from text file
	dblist = open(file+".txt",'r').readlines()
	dblist = [x.strip() for x in dblist]

	return dblist

def sortFile(outFile):
	with open(outFile) as f:
	    sorted_file = sorted(f)

	sorted_file.close()


def wrapBev(bevelPath, targetDB, queryDB, writeDB = False, nMinimizer = nMin, sizeMinimizer = nKmer, filter = nFilter):

		if nKmer >= 32:
				print 'Minimizer size must be less than or equal to 32\n'
				exit()

		s = ' '
		strList = (bevelPath,'-w',str(nMinimizer),'-k',str(sizeMinimizer),'-n',str(nFilter),targetDB,queryDB)
		call = s.join(strList )

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

def analyzeTarget(targetDB, queryDB, oOutput):
	# try:

		# columns=['qDB','qseqID','tDB','qminz']
		#columns=['qseqID','tseqID','qminz']
		columns=['qDB','qseqID','tDB','tseqID','qminz']
		dfall = pd.DataFrame(columns=columns)
		dftseq = pd.DataFrame(['qseqID','tseqID'])
		dftotal = pd.DataFrame(columns=columns)

		# get target DB name
		targetDBName = os.path.splitext(os.path.basename(targetDB))[0]

		# get query DB name
		queryDBName = os.path.splitext(os.path.basename(queryDB))[0]

		lOutput=[]

		for line in wrapBev(bevelPath, targetDB, queryDB):

			result = {}
			line2 = line.strip().split()

			result['qDB'] = queryDBName
			result['qseqID'] = line2[0]
			result['tDB'] = targetDBName
			result['tseqID'] = line2[1].replace('accn|','')
			# result['tminz'] = int(line2[4])
			result['qminz'] = int(line2[5])

			lOutput.append(result.copy())

		if lOutput:
			sortedlist = []

			sortedlist = list(lOutput)

			sortedlist.sort(key=itemgetter('qseqID','tseqID'))

			df = pd.DataFrame()

			# Create dataframe from list of dictionaries
			df = pd.DataFrame(sortedlist)
			dftseq = df.copy(deep=True)

			# delete qminz from dftseq
			dftseq = dftseq.drop('qminz', 1)

			dftseq = dftseq.groupby(['qseqID'],as_index=False).first()

			# Sum query minimizers for lines with same qseqID and same tDB
			df = df.groupby(['qseqID'])['qminz'].sum().reset_index()

			# Sort dataframe by qminz in descending order
			df = df.sort(['qminz'], ascending=[False]).reset_index()

			df = df.groupby(["qseqID"]).apply(lambda x: x[x["qminz"] == x["qminz"].max()])

			df = df.reindex(columns=['qseqID','tseqID','qminz'])

			dfall=dfall.append(df, ignore_index = True)

			dfall = dfall.drop('tseqID', 1)

			dftotal = pd.merge(dfall, dftseq, how='outer', on='qseqID')

		else:
			# oOutput.write('No Minimizers -- ' + queryDBName + ' || ' + targetDBName + '\n')
			return

		# Write results to output file
		dftotal.to_string(oOutput,index=False,header=False)
		oOutput.write('\n')

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

def main(argv):
	targetfile = ''
	queryfile = ''
	try:
		opts, args = getopt.getopt(argv,"ht:q:",["tfile=","qfile="])
	except getopt.GetoptError:
		print 'seqmatch.py - t <targetfile> -q <queryfile>'
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			print 'seqmatch.py -t <targetfile> -q <queryfile>'
			sys.exit()
		elif opt in ("-t", "--tfile"):
				targetfile = arg
		elif opt in ("-q", "--qfile"):
			queryfile = arg

	# Create list of target DBs
	targetDBs = createDBList(targetfile)

	# Create list of query DBs
	queryDBs = createDBList(queryfile)

	total = len(targetDBs)*len(queryDBs)

    # Create output filename which includes time stamp
	outputFile = 'output_' + datetime.now().strftime("%Y%m%d-%H%M%S") + '_' + os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.txt'

	# Create/open output file object
	oFile = open(outputFile, "w")

	count=0
	#Looping through each db list and running bevel against each db
	for queryDB in queryDBs:
		for targetDB in targetDBs:
			count += 1
			# progress=count/float(total)
			update_progress(count/float(total))
			analyzeTarget(targetDB, queryDB, oFile)

	# sortFile(oFile)

	# Done writing to output file so close it
	oFile.close()

	outList = pd.read_csv(outputFile, sep=" ", header = None)
	print(outList)

if __name__ == '__main__': main(sys.argv[1:])