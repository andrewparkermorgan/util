#! /usr/bin/env python

import os
import sys
import csv
import argparse
import re

parser = argparse.ArgumentParser(description = "Simple parser for *.out files from RepeatMasker.")
parser.add_argument(	"-x","--exclude", nargs = "*",
			default = None,
			help = "repeat classes to exclude; regex permitted [default: None]" )
args = parser.parse_args()

## skip first 3 rows
for i in range(0, 4):
	sys.stdin.next()

## open tab-delimted from stdin
#rmsk = csv.reader(sys.stdin, delimiter = " ")
bed = csv.writer(sys.stdout, delimiter = "\t")

for ll in sys.stdin:
	line = ll.strip().split()
	#print "|".join(line)
	chrom = line[4]
	start = int(line[5])
	end = int(line[6])
	strand = line[8]
	#score = int(line[0])
	score = float(line[1])/100
	reptype = line[10]
	repfam = line[9]
	if strand == "C":
		strand = "-"

	keep = True
	if args.exclude is not None:
		for x in args.exclude:
			if re.search(str(x), reptype):
				keep = False

	if keep:
		bed.writerow([ chrom, start, end, reptype, score, strand, repfam ])
