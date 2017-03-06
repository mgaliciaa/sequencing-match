#! /usr/bin/python

import os
import sys
import string

from datetime import datetime
from subprocess import Popen

bevelPath = "/share/biocore/internal_projects/seqmatch/bevel/bin/bevel"
targetDB = "/share/biocore/internal_projects/seqmatch/genomes/59814.6/59814.6.fna"
queryDB = "/share/biocore/internal_projects/seqmatch/03-SpadesAssemblies/26-3_S36/26-3_S36.Scaffolds.fna"
outputFile = "output_" + datetime.now().strftime("%Y%m%d-%H%M%S") + ".txt"

p = Popen([bevelPath, targetDB, queryDB])

#p = Popen([bevelPath, targetDB, queryDB, '> ' + outputFile])

p.communicate()
if p.returncode != 0:
	raise OSError('bevel error')

