#! /usr/bin/env python

##	--- dna.py --- ##
##	Date: 18 June 2014
##	Updated: 14 Oct 2014
##	Purpose: utility functions for DNA sequences

import re

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

DNACODE = { "A":"T", "C":"G", "G":"C", "T":"A", "N":"N" }

## undo lowercase repeat-masking on DNA sequence
def unmask(dna):
	dnastr = "".join(dna)
	return dnastr.upper()

## mask lowercase repeats to Ns (or whatever)
def mask(dna, repeatchar = "N"):
	dnastr = "".join(dna)
	return re.sub(r"[acgtn]", "N", dnastr)

def _as_seq_object(dna, alphabet = IUPAC.ambiguous_dna):

	if not isinstance(dna, Seq):
		dna = Seq(dna, alphabet)

	return dna

## reverse complement, including on-the-fly conversion to Bio.Seq
def revcomp(dna, alphabet = IUPAC.unambiguous_dna):

	dna = _as_seq_object(dna, alphabet)
	return str(dna.reverse_complement())

def has_ambiguities(dna, alphabet = IUPAC.ambiguous_dna):

	dna = _as_seq_object(dna, alphabet)
	flag = any([ nt in dna for nt in alphabet.letters[4:] ])
	return flag


## DEPRECATED VERSION
# def revcomp(dna, complementarity = DNACODE):
#
# 	dna = unmask(dna)
# 	newdna = []
# 	for bp in dna:
# 		if bp not in complementarity:
# 			newdna.append("N")
# 		else:
# 			newdna.append( complementarity[bp] )
#
# 	return "".join( newdna[::-1] )
