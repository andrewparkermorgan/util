#! /usr/bin/env python
"""
vcf_to_tped.py
Conert VCF file to PLINK *.tped format, preserving phase.
"""

from __future__ import print_function

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
parser.add_argument("--keep", type = argparse.FileType("rU"),
					required = False,
					help = "list of individuals to keep [default: all]" )
parser.add_argument("--exclude", type = argparse.FileType("rU"),
					required = False,
					help = "list of individuals to exclude; applied after '--keep' [default: all]" )
parser.add_argument("--missing-char",
					default = 0,
					help = "missing-genotype character [default: %(default)s]" )
parser.add_argument("--status-every", type = int,
					default = 0,
					help = "report status ever X markers (0 for no status updates) [default: %(default)d]" )
parser.add_argument("--header", action = "store_true",
					default = False,
					help = "include a header line with column names [default: no header]" )
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

def parse_ids(infile):
	ids = []
	for line in infile:
		if line.startswith("#"):
			continue
		pieces = line.split()
		if len(pieces) > 1:
			fid = pieces.pop(0)
			iid = pieces.pop(0)
			ids.append(iid)
		else:
			ids.append(pieces.pop())
	return set(ids)


## check if input is gzipped
if args.vcf.endswith(".gz") or args.vcf.endswith(".Z"):
	infile = gzip.open(args.vcf, "r")
elif args.vcf == "-":
	infile = sys.stdin
else:
	infile = open(args.vcf, "r")

## read keep/exclude lists
keepers = set()
if args.keep:
	keepers |= parse_ids(args.keep)
if args.exclude:
	keepers -= parse_ids(args.exclude)

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
		if not len(keepers):
			keepers |= set(samples)
		if args.header:
			print("#chr","marker","cM","pos"," ".join([ "{}.1 {}.2".format(s,s) for s in samples if s in keepers ]), file = sys.stdout)

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
		alleles_txt = " ".join([ " ".join(a) for i,a in enumerate(alleles) if samples[i] in keepers ])
		alleles_txt = alleles_txt.replace(".", str(args.missing_char))

		cm = 0
		if gmap is not None and marker in gmap:
			cm = gmap[marker][1]

		print(chrom, marker, cm_rescale*cm, pos, alleles_txt, file = sys.stdout)

	i += 1
	if args.status_every > 0:
		if not (i % args.status_every):
			sys.stderr.write("\t... {} markers processed\n".format(i))
