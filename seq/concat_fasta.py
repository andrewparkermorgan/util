#! /usr/bin/env python

import os
import sys
import argparse

from Bio import SeqIO
from myio import common

parser = argparse.ArgumentParser(description = "Concatenate two or more fastas and deal with missing data.")
parser.add_argument(	"-u","--union", action = "store_true",
 						help = "Take the 'union' of sequences across all input files [default: do intersection]" )
parser.add_argument(	"-w","--width", default = 50,
 						help = "wrap sequence output to this many characters per line [default %(default)d]" )
parser.add_argument(	"fastas", nargs = "+", metavar = "FASTA",
						type = common.readable_or_stdin_handle,
						help = "one or more fastas to concatenate; NB: order of sequences wont' be preserved" )
args = parser.parse_args()

def longest_sequence(seqdict):
	return max([ len(seq) for name,seq in seqdict.iteritems() ])

def fold(x, width = 50):
	for i in range(0, len(x), width):
		yield x[i:(i+50)]

seqnames = set()
all_seqs = []
for fa in args.fastas:
	these_seqs = SeqIO.to_dict( SeqIO.parse(fa, "fasta") )
	all_seqs.append(these_seqs)
	if args.union:
		seqnames = seqnames.union( these_seqs.keys() )
	else:
		seqnames = seqnames.intersection( these_seqs.keys() )

for s in sorted(seqnames):
	final_seq = ""
	for fa in all_seqs:
		longest = longest_sequence(fa)
		if s in fa:
			this_seq = str(fa[s].seq).strip()
		else:
			this_seq = "-"*longest
		final_seq += this_seq
	print ">{}".format(s)
	for line in fold(final_seq):
		print line
