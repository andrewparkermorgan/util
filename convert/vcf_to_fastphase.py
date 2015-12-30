#! /usr/bin/env python

from __future__ import print_function # for python3 compatibility

import os
import sys
import collections
import numpy as np
import argparse

from slimvcf import SlimVCF

parser = argparse.ArgumentParser(description = "Convert phased VCF file to fastPHASE-like output format.")
parser.add_argument(	"-m","--map", type = argparse.FileType("w"),
						default = "./map.inp",
						help = "file path for SNP map [default: %(default)s]" )
args = parser.parse_args()

## initialize VCF reader
vcf = SlimVCF(sys.stdin)

## initialize container for phased haplotypes
haps = collections.OrderedDict()
sites = []
for s in vcf.samples:
	haps.update({ s: ([],[]) })

sys.stderr.write("Reading phased genotypes for {} samples from VCF file...\n".format(len(vcf.samples)))

## loop on sites in VCF file
j = 0
for site in vcf:

	sites.append( ("{}.{}".format(site.chr, site.pos), site.chr, site.pos, 1, 2) )
	for i in range(0, len(vcf.samples)):
		
		a1 = int(site.geno_matrix[i,0])
		a2 = int(site.geno_matrix[i,1])
		
		if a1 == -1:
			a1 = "?"
		else:
			a1 += 1
		
		if a2 == -1:
			a2 = "?"
		else:
			a2 += 1

		haps[ vcf.samples[i] ][0].append(str(a1))
		haps[ vcf.samples[i] ][1].append(str(a2))

	j += 1
	if not j % 10000:
		sys.stderr.write("\t{} sites processed\n".format(j))

sys.stderr.write("\t{} total sites.\n".format(j))
sys.stderr.write("Map contains {} sites.\n".format(len(sites)))
sys.stderr.write("\nWriting fastPHASE-format haplotypes...\n")
print("********************************************")
print("*                                          *")
print("*      Output from fastPHASE 1.3.0c        *")
print("*      Code by P Scheet                    *")
print("*                                          *")
print("********************************************")
print("BEGIN COMMAND_LINE")
print("./fastPHASE --bogus")
print("END COMMAND_LINE")
print("")
print("BEGIN GENOTYPES")
for s in vcf.samples:
	sys.stderr.write("\t{}\n".format(s))
	print("{}".format(s))
	print(" ".join(haps[s][0]))
	print(" ".join(haps[s][1]))
print("END GENOTYPES")

sys.stderr.write("\nWriting map file to {}...\n".format(args.map))
for s in sites:
	print(s[0], s[1], s[2], s[3], s[4], file = args.map)

sys.stderr.write("Done.\n\n")