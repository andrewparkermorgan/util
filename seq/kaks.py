#! /usr/bin/env python

## TODO: move code for stripping gaps, cutting at stop codon to separate module

import os
import sys
import argparse
import tempfile
import itertools # for combinations
import collections # for Counter
import re

from myio import common

from Bio import AlignIO
from Bio.Alphabet import Gapped, IUPAC
from Bio.Data.CodonTable import standard_dna_table
from Bio.Phylo.PAML import yn00

parser = argparse.ArgumentParser(description = "Simple wrapper for PAML::yn00, to compute Ka/Ks from codon alignment.")
parser.add_argument( "-s","--sequences", type = common.readable_file,
					help = "file of sequences (fasta format)" )
parser.add_argument( "-w", "--working-dir", type = common.writeable_dir,
 					required = False,
					help = "working directory for RAxML output [default: an auto-generated temp directory]" )
parser.add_argument( "-m", "--method", choices = ["YN00","NG86"],
 					default = "YN00",
					help = "algorithm for computing Ka/Ks, either 'NG86' (Nei-Gojobori) or 'YN00' (Yang-Nielsen) [default: %(default)s]" )
parser.add_argument( "-g", "--trim-gaps", type = float,
 					default = 1.0,
					help = "trim triplets with more than this fraction gap characters [default: %(default)f]" )
args = parser.parse_args()

## set working directory
wd = args.working_dir
if not wd:
	wd = tempfile.mkdtemp()

## make sure this is a valid codon alignment (by force if needed) and convert to phylip-sequential
prefix = re.sub(r"\.\w+$", "", os.path.basename(args.sequences))
seq_converted = os.path.join(wd, "{}.phylip".format(prefix))
fa = AlignIO.read(args.sequences, "fasta")
untrimmed_len = len(fa[0])
i = 0
while i <= len(fa[0]) - 3:
	counts = collections.Counter(fa[ :,i ])
	#print i/3, float(counts["-"])/len(fa)
	if float(counts["-"])/len(fa) >= args.trim_gaps:
		fa = fa[ :,:i, ] + fa[ :,(i+3): ]
		#print "just trimmed a gappy codon:", len(fa[0])
	else:
		i += 3
first_stop = len(fa[0]) - (len(fa[0]) % 3)
for s in fa:
	s.id = re.sub(r"\/\d+\-\d+$","", s.id)[ 0:8 ]
	for stop in standard_dna_table.stop_codons:
		has_stop = [ x.start() for x in re.finditer(stop, str(s.seq)) if not (x.start() % 3) ]
		if len(has_stop):
			#sys.stderr.write("\nfound stop codon at {} in sequence '{}'\n".format(has_stop[0], s.id))
			first_stop = has_stop[0]
fa = fa[ :, 0:first_stop ]
trimmed_len = len(fa[0])

with open(seq_converted, "w") as ph:
	AlignIO.write(fa, ph, "phylip-sequential")

## print some diagnostics to stderr
sys.stderr.write("\nWorking directory: {}\n".format(wd))
sys.stderr.write("Input file: {} (length {})\n".format(args.sequences, untrimmed_len))
sys.stderr.write("Converted to phylip: {} (length {})\n\n".format(seq_converted, trimmed_len))

## run yn00
yn = yn00.Yn00(alignment = seq_converted, working_dir = wd, out_file = "results.out")
rez = yn.run(verbose = False)

## get results for all sequence pairs
ph = AlignIO.read(seq_converted, "phylip-relaxed")
ids = [ s.id for s in ph ]
for s1,s2 in itertools.combinations(ids, 2):
	omega = rez[s1][s2][args.method]["omega"]
	if omega < 0 or omega > 10:
		omega = "NA"
	print s1, s2, omega, abs(rez[s1][s2][args.method]["dN"]), abs(rez[s1][s2][args.method]["dS"])
