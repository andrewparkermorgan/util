#! /usr/bin/env python

import os
import sys
import csv
import argparse

parser = argparse.ArgumentParser(description = "Simple parser for *.psl output from blat.")
parser.add_argument(	"-w","--minsize", type = int,
			default = 0,
			help = "minimum hit size (match+mismatch+repmatch) [default: %(default)d]" )
args = parser.parse_args()

## skip first 5 rows
for i in range(0, 5):
	sys.stdin.next()

## open tab-delimted from stdin
psl = csv.reader(sys.stdin, delimiter = "\t")
bed = csv.writer(sys.stdout, delimiter = "\t")
for line in psl:
	if sum([int(x) for x in line[0:3]]) < args.minsize:
		continue
	qname = line[9]
	tname = line[13]
	strand = line[8]
	qsizes = line[18].strip(",").split(",")
	qstarts = line[19].strip(",").split(",")
	tstarts = line[20].strip(",").split(",")
	qends = [ int(s)+int(l) for (s,l) in zip(qstarts, qsizes) ]
	tends = [ int(s)+int(l) for (s,l) in zip(tstarts, qsizes) ]
	for i in range(0, len(qstarts)):
		bed.writerow([ qname, qstarts[i], qends[i], "block"+str(i), ".", strand, tname, tstarts[i], tends[i] ])
