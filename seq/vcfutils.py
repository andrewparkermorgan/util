import os
import sys

import pysam

def parse_region(reg):
	pieces = reg.split(":")
	if len(pieces) == 2:
		(start, end) = pieces.pop().split("-")
		return (pieces[0], int(start), int(end))
	else:
		return (pieces[0],)

def geno_to_str(geno, short = False):
	## geno is a dictionary of lists:
	## { 'sample1': ['1','/','1'], 'sample2': ['0','/','0'], ... }
	geno_str = []
	if not short:
		geno_str = [ "".join([ str(i) for i in a ]) for s,a in geno.iteritems() ]
	else:
		geno_str = [ a[0] for s,a in geno.iteritems() ]

	return geno_str

def is_het(geno):
	return geno[0] != geno[2]

def is_polymorphic_site(geno):
	geno_str = geno_to_str(geno)
	return len(set(geno_str)) > 1

def is_variant_site(geno):
	geno.update( { ".": [0,"/",0] })
	return is_polymorphic_site(geno)

def is_private_variant(samples, geno):
	include = { s: g for s,g in geno.iteritems() if s in samples }
	exclude = { s: g for s,g in geno.iteritems() if s not in samples }
	geno_str_in = set(geno_to_str(include))
	geno_str_ex = set(geno_to_str(exclude))
	return ( len(geno_str_in) == 1 and len(geno_str_in & geno_str_ex) == 0 )

def get_geno(site, samples = None):
	geno = { s: site[s]["GT"][0] for s in samples }
	return geno

def get_alleles(site):
	alleles = [ site.ref ] + site.alt
	return alleles

def parse_ensembl_csq(site):
	if "CSQ" in site.info:
		csq = []
		csq_pieces = ",".join(site.info["CSQ"])
		csq_alleles = csq_pieces.split("+")
		for a in csq_alleles:
			pieces = a.split(":")
			csq.append({ "chr": site.contig, "pos": site.pos, "transcript": pieces[0], "gene": pieces[1], "consequence": pieces[2] })
		return csq
	else:
		return None

def has_csq(site, so_terms = []):

	csq = parse_ensembl_csq(site)
	if (not csq) or (not len(so_terms)):
		return False

	flags = []
	for s in so_terms:
		flags.append( reduce(	lambda x,y: x or y,
					[ c["consequence"].upper() == s.upper() for c in csq ],
					False ) )

	return reduce(lambda x,y: x or y, flags, False)
