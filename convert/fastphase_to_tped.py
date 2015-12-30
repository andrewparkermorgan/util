#! /usr/bin/env python

## --- fastphase_to_tped.py --- ##
##	Date: 2 March 2015
##	Purpose: transpose fastPHASE output and attach map to make a tped file

import os
import sys
import csv
import argparse

parser = argparse.ArgumentParser(description = "Converter for fastPHASE output. NB: no error checking is performed; user beware.")
parser.add_argument("-m", "--map", type = argparse.FileType("rU"),
					required = True,
					help = "marker map file (plink *.map or *.bim format)" )
parser.add_argument("phased", type = argparse.FileType("rU"),
					help = "output file from fastPHASE" )
args = parser.parse_args()

def plinkify(x):
	if x == "?":
		x = 0
	return str(int(x))

## read phased genotypes
cols = []
for line in args.phased:

	if not line.startswith("BEGIN GENOTYPES"):
		continue
	else:
		for chrom in args.phased:
			if chrom.startswith("#") or chrom.startswith("END GENOTYPES"):
				continue
			else:
				this_chr = chrom.split()
				cols.append(this_chr)

## read marker map
snf = csv.Sniffer().sniff(args.map.read(1024))
args.map.seek(0)
markers = csv.reader(args.map, snf)
i = 0
for line in markers:
	print " ".join(line[0:4]), " ".join([ plinkify(g[i]) for g in cols ])
	i += 1
