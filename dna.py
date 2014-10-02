#! /usr/bin/env python

##	--- dna.py --- ##
##	Date: 18 June 2014
##	Updated: NA
##	Purpose: utility functions for DNA sequences

import re

DNACODE = { "A":"T", "C":"G", "G":"C", "T":"A", "N":"N" }

## undo lowercase repeat-masking on DNA sequence
def unmask(dna):
	dnastr = "".join(dna)
	return dnastr.upper()

## mask lowercase repeats to Ns (or whatever)
def mask(dna, repeatchar = "N"):
	dnastr = "".join(dna)
	return re.sub(r"[acgtn]", "N", dnastr)

def revcomp(dna, complementarity = DNACODE):

	dna = unmask(dna)
	newdna = []
	for bp in dna:
		if bp not in complementarity:
			newdna.append("N")
		else:
			newdna.append( complementarity[bp] )

	return "".join( newdna[::-1] )