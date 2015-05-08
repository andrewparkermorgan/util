#! /usr/bin/env python

import os
import sys
import argparse

from myio import common

from Bio import SeqIO
from cogent.align.algorithm import nw_align, sw_align

parser = argparse.ArgumentParser(description = "Simple pairwise alignment of two DNA or protein sequences. Local alignment by default.")
parser.add_argument(	"-f","--fasta", type = common.readable_or_stdin_handle,
						default = "-",
						help = "fasta file containing sequences; first 2 will be aligned [default: stdin]" )
parser.add_argument(	"-a","--algorithm", choices = ["sw","nw"],
						default = "sw",
						help = "algorithm to use (nw = Needleman-Wunsch/global, sw = Smith-Waterman/local) [default: %(default)s]" )
args = parser.parse_args()

seqs = []
fa = SeqIO.parse(args.fasta, "fasta")
for next_seq in fa:
		seqs.append( (next_seq.id, str(next_seq.seq.upper())) )

for j in range(0, len(seqs)/2):
	if args.algorithm == "nw":
		aln_fn = nw_align
	else:
		aln_fn = sw_align
	aligned = aln_fn(seqs[2*j][1], seqs[2*j+1][1])
	for i in range(0, len(aligned)):
		print ">{}".format(seqs[2*j+i][0])
		print aligned[i]
