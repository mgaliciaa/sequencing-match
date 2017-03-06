#! /usr/bin/python

import os
import sys
import string
import signal

from datetime import datetime
from subprocess import Popen
from subprocess import PIPE

bevelPath = "/share/biocore/internal_projects/seqmatch/bevel/bin/bevel"

def wrapBev(bevelPath, targetDB, queryDB, writeDB = False, nMinimizer = 100, sizeMinimizer = 17):
	call = bevelPath
	if writeDB:
		call = call + " -d"
	call = call + " -w " + nMinimizer + " -k " + sizeMinimizer
	call = call + " " + targetDB + " " + queryDB
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

def analyzeTarget(targetDB, queryDB, outputFile):
#do stuff
	try:
		oFile = open(outputFile, "w")
		i = 0
		for line in wrapBev(bevelPath, targetDB, queryDB):
			i += 1
			oFile.write(line)
			line2 = line.strip().split()
	except:
		print ("Finished processing file, read ", i ,  " lines\n")
		oFile.close()

def main():
	targetDB = "/share/biocore/internal_projects/seqmatch/genomes/59814.6/59814.6.fna"
	queryDB = "/share/biocore/internal_projects/seqmatch/03-SpadesAssemblies/26-3_S36/26-3_S36.Scaffolds.fna"
	outputFile = "output_" + datetime.now().strftime("%Y%m%d-%H%M%S") + ".txt"
	analyzeTarget(targetDB, queryDB, outputFile)



if __name__ == '__main__':
	main()