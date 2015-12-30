#! /usr/bin/env python
"""
split_vcf.py
Split a VCF file into smaller files, each with a fixed number of variants.
"""

from __future__ import print_function # for py3 compatibility

import os
import sys
import gzip
import argparse

parser = argparse.ArgumentParser(description = "Split a VCF file into smaller files, each with a fixed number of variants.")
parser.add_argument(	"-i","--infile", type = argparse.FileType("rU"),
						default = sys.stdin,
						help = "VCF file/stream to process, uncompressed [default: stdin]" )
parser.add_argument(	"-w","--window", type = int,
						default = 1000,
						help = "window size, in number of sites [default: %(default)d]" )
parser.add_argument(	"-s","--step", type = int,
						default = 1000,
						help = "step size, in number of sites [default: %(default)d]" )
parser.add_argument(	"-o","--out",
						default = "out",
						help = "prefix for output files [default: %(default)s]" )
parser.add_argument(	"-d","--digits", type = int,
						default = 4,
						help = "log10(estimated max number of windows) [default: %(default)d]" )
parser.add_argument(	"-z","--gzip", action = "store_true",
						default = False,
						help = "output gzipped files [default: %(default)s]" )
args = parser.parse_args()

if args.step > args.window:
	sys.exit("Step size ({}) must be less than window size ({}).".format(args.step, args.window))

def block_writer(ff):
	if args.gzip:
		return gzip.open(ff + ".gz", "wb")
	else:
		return open(ff, "w")

def write_block(sites, block):
	ff = (args.out + "{:0" + str(args.digits) + "}.vcf").format(block)
	lastpos = ": ".join(sites[-1:].pop().split()[0:2])
	with block_writer(ff) as outfile:
		print("\t... {} sites read: writing block {} to <{}>...\t@ {}".format(i, block, outfile.name, lastpos), file = sys.stderr)
		for h in header:
			print(h, file = outfile)
		for s in sites:
			print(s, file = outfile)

i = 0
block = 0
header_done = False
header = []
sites = []
print("Reading VCF header...", file = sys.stderr)
for line in args.infile:
	if line.startswith("#"):
		if not header_done:
			header.append(line.strip())
		else:
			continue
	else:
		if not header_done:
			print("Reading sites...", file = sys.stderr)
		header_done = True
		i += 1
		sites.append(line.strip())
		if not len(sites) % args.window:
			write_block(sites, block)
			block += 1
			if args.window == args.step:
				sites = []
			else:
				sites = sites[-(args.window-args.step):]

if len(sites):
	write_block(sites, block)
	block += 1

print("Done: created {} files.\n".format(block), file = sys.stderr)