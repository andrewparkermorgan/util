#! /usr/bin/env python
"""
nocpg.py
Remove SNV sites from a VCF file if they fall at a CpG site in the reference.
Ignores effect of adjacent SNP which changes the C to something else, and ignores indels, since the CpG distinction
is only meaningful in the context of a SNV.
"""

from __future__ import print_function # for py3 compatibility

import os
import sys
import argparse

import pyfasta

parser = argparse.ArgumentParser(description = "Remove CpG sites from a VCF file streamed to stdin.")
parser.add_argument(	"-f","--fasta",
						default = "$GENOMES/mm10/mm10_all.fa",
						help = "path to reference genome compatible with the one used to make this VCF [default: %(default)s]" )
args = parser.parse_args()

args.fasta = os.path.expandvars(os.path.expanduser(args.fasta))
fa = pyfasta.Fasta(args.fasta)

cpgs = 0
for line in sys.stdin:
	do_print = True
	if not line.startswith("#"):
		chrom, pos, var_id, ref, alt, notused = line.strip().split(None, 5)
		pos = int(pos)
		alts = alt.split(",")
		if len(ref) == 1 and all(len(a) == 1 for a in alts):
			dinuc = fa[chrom][ (pos-1):(pos+1) ].upper()
			if dinuc == "CG":
				cpgs += 1
				do_print = False
	
	if do_print:
		print(line.strip())

print("Filtered {} CpG sites.\n".format(cpgs), file = sys.stderr)