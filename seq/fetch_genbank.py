#! /usr/bin/env python

import os
import sys
import argparse

from Bio import Entrez
from Bio import SeqIO

parser = argparse.ArgumentParser(description  = "Query a bunch of NCBI accession numbers and dump result as fasta.")
parser.add_argument(	"-e","--email",
						default = "apm@email.unc.edu",
						help = "declare your identity to NCBI using this email [default: %(default)s]" )
parser.add_argument(	"-i","--ids", type = argparse.FileType("rU"),
						help = "file containing accession numbers (one per line)" )
parser.add_argument(	"-d","--db",
						default = "nucleotide",
						help = "NCBI database name to query [default: %(default)s]" )
args = parser.parse_args()

idlist = []
for line in args.ids:
	ll = line.strip()
	if ll.startswith("#") or not len(ll):
		continue
	idlist.append(ll)

if not len(idlist):
	sys.sdterr.write("Oops: list of accession numbers is empty.\n")
	sys.exit(1)

Entrez.email = args.email
for i in idlist:
	sys.stderr.write("--- {} ---\n".format(i))
	try:
		handle = Entrez.efetch(db = args.db, id = str(i), rettype = "gb", retmode = "text")
		record = SeqIO.read(handle, "genbank")
		print record.format("fasta")
	except Exception as e:
		sys.stderr.write(str(e) + "\n")
