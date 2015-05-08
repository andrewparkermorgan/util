#! /usr/bin/env python

import os
import sys

from Bio import SeqIO

fa = SeqIO.parse(sys.stdin, "fasta")
for seq in fa:
	newseq = seq.reverse_complement()
	newseq.id = seq.id
	newseq.description = ""
	sys.stdout.write(newseq.format("fasta"))
