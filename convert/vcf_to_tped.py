#! /usr/bin/env python

## --- vcf_to_tped.py -- ##
## 	Date: 9 March 2015
##	Purpose: convert vcf to tped format, respecting phase and emitting either character (A/C/T/G) or numeric (0/1) alleles

import os
import sys
import argparse
import gzip

## parse command-line arguments
parser = argparse.ArgumentParser(description = "Lightweight and dumb utility to convert (phased) VCF to plink TPED format.")
parser.add_argument("--recode", action = "store_true",
					help = "recode alleles to nucleotide (A/C/T/G) [default: use 1/2 coding]" )
parser.add_argument("--map", type = argparse.FileType("rU"),
					required = False,
					help = "map file in plink format (chr marker cM pos)" )
parser.add_argument("--missing-char",
					default = 0,
					help = "missing-genotype character [default: %(default)s]" )
parser.add_argument("--status-every", type = int,
					default = 0,
					help = "report status ever X markers (0 for no status updates) [default: %(default)d]" )
parser.add_argument("vcf",
					help = "path to VCF file (can be gzipped)" )
args = parser.parse_args()

def parse_geno(g, calls = None):
	geno = g.split(":").pop(0)
	if "|" in geno:
		alleles = geno.split("|")
	elif "/" in geno:
		alleles = geno.split("/")
	else:
		raise Exception("Allele delimiter not recognized: reported genotype was '{}'.".format(geno))

	if calls:
		alleles = [ calls[ int(a) ] for a in alleles ]
	else:
		alleles = [ str(int(a)) for a in alleles ]

	return alleles

## check if input is gzipped
if args.vcf.endswith(".gz") or args.vcf.endswith(".Z"):
	infile = gzip.open(args.vcf, "r")
elif args.vcf == "-":
	infile = sys.stdin
else:
	infile = open(args.vcf, "r")

gmap = None
cm_rescale = 1
if args.map:
	cms = []
	gmap = {}
	for line in args.map:
		(chrom, marker, cm, pos) = line.split()[0:4]
		cm = float(cm)
		cms.append(cm)
		gmap[marker] = (chrom, cm, pos)
	sys.stderr.write("Loaded {} markers from genetic map.\n".format(len(gmap.keys())))
	if max(cms) < 20:
		cm_rescale = 100
		sys.stderr.write("Rescaling genetic distances by {}.\n".format(cm_rescale))

## loop on lines in vcf
i = 0
for line in infile:

	## skip header lines
	if line.startswith("##"):
		continue
	
	if line.startswith("#CHROM"):
		## this is final header line; parse samples
		samples = line.split()[9:]

	else:
		
		## this is a variant entry
		pieces = line.split()
		(chrom, pos, marker, ref, alt, qual, filtr, info, fmt) = pieces[:9]
		geno = pieces[9:]

		## check biallelicity
		nalt = len(alt.split(","))
		if nalt > 1:
			sys.stderr.write("Warning: {} ({}:{}) is multiallelic and will be skipped.\n".format(marker, chrom, pos))
			continue

		## parse genotypes
		calls = None
		if args.recode:
			calls = [ref].extend(alt.split(","))

		alleles = [ parse_geno(g, calls) for g in geno ]
		alleles_txt = " ".join([ " ".join(a) for a in alleles ])
		alleles_txt = alleles_txt.replace(".", str(args.missing_char))

		cm = 0
		if gmap is not None and marker in gmap:
			cm = gmap[marker][1]

		print chrom, marker, cm_rescale*cm, pos, alleles_txt

	i += 1
	if args.status_every > 0:
		if not (i % args.status_every):
			sys.stderr.write("\t... {} markers processed\n".format(i))

