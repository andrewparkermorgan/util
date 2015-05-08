#! /usr/bin/env python

import os
import sys
import argparse

parser = argparse.ArgumentParser(description = "Command-line translation of cDNA to protein sequence.")
parser.add_argument(	"-f", "--frame", type = int,
						default = 1,
						help = "translation frame (+1,+2,+3,-1,-2,-3) [default: %(default)d]" )
args = parser.parse_args()

from Bio import SeqIO

fa = SeqIO.parse(sys.stdin, "fasta")
for seq in fa:
	codons = seq.seq.ungap("-")
	print ">" + seq.id
	if args.frame > 0:
		print codons[ (args.frame-1): ].translate()
	else:
		print codons.reverse_complement()[ (abs(args.frame)-1): ].translate()
