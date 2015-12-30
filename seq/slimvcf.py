#! /usr/bin/env python
## --- slimvcf.py -- ##
##	Date:		28 July 2015
##	Purpose:	A simple parser for VCF files.  Little or no validation, but no problems with VCF format spec either.

from __future__ import print_function # python3 compatibility

import os
import sys
import re
import numpy as np
from collections import Iterable, OrderedDict
from copy import deepcopy

class VariantSite:

	MISSING_ALLELES = [ "X", "." ]

	def _float_or_missing(self, x):
		try:
			return float(x)
		except:
			return -1
	def _numeric_to_missing(self, x):
		if x < 0:
			return "."
		else:
			return(str(x))

	def __init__(self, kwargs):

		## basic properties
		self.chr = None
		self.pos = None
		self.id = None
		self.ref = None
		self.alt = []
		self.qual = None
		self.filter = []
		self.info = None
		self.format = None
		self.geno = None
		self.phased = None
		self.geno_matrix = None

		## a few useful derived properties
		self.nalleles = 1
		self.is_indel = False

		## get properties from function args
		if "chr" in kwargs:
			self.chr = kwargs["chr"]
		if "pos" in kwargs:
			self.pos = max(int(kwargs["pos"]),1)
		if "id" in kwargs:
			self.id = kwargs["id"]
		if "ref" in kwargs:
			self.ref = str(kwargs["ref"])
		if "alt" in kwargs and isinstance(kwargs["alt"], Iterable):
			self.alt = [ str(a) for a in kwargs["alt"] ]
		if "qual" in kwargs:
			self.qual = self._float_or_missing(kwargs["qual"])
		if "filter" in kwargs:
			self.filter = [ str(f) for f in kwargs["filter"] ]
		if "info" in kwargs:
			self.info = kwargs["info"]
		if "format" in kwargs:
			self.format = kwargs["format"]
		if "geno" in kwargs:
			self.geno = kwargs["geno"]
		if "phased" in kwargs:
			self.phased = kwargs["phased"]

		## set the derived properties
		if self.alt is not None:
			nonmissing = [ a for a in self.alt if a not in self.MISSING_ALLELES ]
			self.nalleles = len(nonmissing) + 1
		if len(self.alt) and len(self.ref):
			nonmissing = [ a for a in self.alt if a not in self.MISSING_ALLELES ]
			if max(len(a) for a in [self.ref] + nonmissing) > 1:
				self.is_indel = True

	def _dict_to_keyval(self, x):
		out = []
		for k,v in x.iteritems():
			if v is not None:
				out.append("{}={}".format(k,v))
			else:
				out.append(k)
		return ";".join(out)

	def __str__(self, samples = None):

		out = [	str(self.chr), self._numeric_to_missing(self.pos), str(self.id), str(self.ref), str(",".join(self.alt)),
				self._numeric_to_missing(self.qual), str(";".join(self.filter)), self._dict_to_keyval(self.info), ":".join(self.format) ]
		if samples is None:
			samples = self.geno.keys()
		for s in samples:
			if s in self.geno:
				out = out + [ ":".join(self.geno[s].values()) ]

		return "\t".join(out)

	def _copy_geno_mat(self, sample_idx = None):

		if self.geno_matrix is None:
			raise ValueError("No genotypes matrix present.")
		g = self.geno_matrix.copy()
		if sample_idx is not None:
			g = g[ sample_idx,: ]
		return(g)

	def is_polymorphic(self, sample_idx = None):

		## make ourselves a local copy
		g = self._copy_geno_mat(sample_idx)
		## scrub missing values
		g[ g == -1 ] = 0
		## monomorphic site iff max=min
		return np.min(g) != np.max(g)

	def is_filtered(self, filters = [".","PASS"]):
		return len(self.filter) and len([ f for f in self.filter if f not in filters ])

	def mac(self, sample_idx = None):
		## make ourselves a local copy
		g = self._copy_geno_mat(sample_idx)
		## ignore samples with missing data...
		g = g[ np.min(g, 1) >= 0,: ]
		if g.shape[0] > 0:
			return float(np.sum(g[ np.min(g, 1) > 0,: ] == 0))
		else:
			return 0

	def maf(self, sample_idx = None):
		## make ourselves a local copy
		g = self._copy_geno_mat(sample_idx)
		## ignore samples with missing data...
		g = g[ np.min(g, 1) >= 0,: ]
		if g.shape[0] > 0:
			return float(np.sum(g[ np.min(g, 1) > 0,: ] == 0))/(2.0*g.shape[0])
		else:
			return 0.0

	def slice(self, sample_idx = None, het_missing = False):

		new_site = deepcopy(self)
		if sample_idx is not None:
			keep = [ self.geno.keys()[i] for i in sample_idx ]
		else:
			keep = self.geno.keys()
		new_geno = OrderedDict()
		for s in keep:
			new_geno.update({ s: self.geno[s] })
		new_site.geno = new_geno
		new_site.geno_matrix = self._copy_geno_mat(sample_idx)
		if het_missing:
			g = new_site.geno_matrix
			hets = g[:,0] != g[:,1]
			g[hets,:] = -1
			new_site.geno_matrix = g
		return new_site

	def is_biallelic(self, sample_idx = None):
		#g = self._copy_geno_mat(sample_idx)
		## make ourselves a local copy
		g = self._copy_geno_mat(sample_idx)
		## scrub missing values
		g[ g == -1 ] = 0
		return np.max(g) <= 1

	def dosages(self, sample_idx = None):
		if not self.is_biallelic(sample_idx):
			raise ValueErorr("This site not biallelic; dosages are bogus.")
		else:
			g = self._copy_geno_mat(sample_idx)
			## make ourselves a local copy
			g = self._copy_geno_mat(sample_idx)
			## scrub missing values
			g[ g == -1 ] = 0
			return np.sum(g, 1)

	def allele_counts(self, samples, ploidy = None, fold = False):
		if ploidy is None:
			ploidy = np.array([ 2 for s in samples ])
			ploidy = np.reshape(ploidy, (len(samples), 1))
		else:
			if len(ploidy) != len(samples):
				raise ValueError("If ploidy is specified, it should be a list with same length as samples.")
			else:
				ploidy = np.array(ploidy)
				ploidy = np.reshape(ploidy, (len(samples), 1))
		avail_samples = self.geno.keys()
		sample_idx = [ avail_samples.index(s) for s in samples ]
		g = self._copy_geno_mat(sample_idx)
		nonmissing = np.sum((g != -1) * ploidy/2.0)
		g[ g == -1 ] = 0
		count = np.sum(g*ploidy/2.0)
		if fold:
			count = min(count, (np.sum(ploidy)-count))
		return nonmissing, count

	def apply_expr(self, expr):
		ops = re.findall(r"\s+AND\s+|\s+OR\s+", str(expr))
		exprs = re.split(r"\s+AND\s+|\s+OR\s+", str(expr))
		flags = []
		for e in exprs:
			m = re.search(r"[=><\!]+", e)
			if m is not None:
				op = m.group(0)
				tag, val = re.split(r"[=><\!]+", e)
				if not tag in self.info:
					raise ValueError("Can't find tag '{}'' in INFO field.".format(tag))
				if op == "=":
					flags.append( str(self.info[tag]) == str(val) )
				elif op == ">":
					flags.append( float(self.info[tag]) > float(val) )
				elif op == "<":
					flags.append( float(self.info[tag]) < float(val) )
				elif op == ">=":
					flags.append( float(self.info[tag]) >= float(val) )
				elif op == "<=":
					flags.append( float(self.info[tag]) <= float(val) )
				else:
					raise ValueError("Comparison operator in filtering expression not recognized.")
			else:
				raise ValueError("Filter expression lacks valid comparison operator.")
		if len(flags) == len(exprs):
			if len(flags) == 0:
				return True
			elif len(flags) == 1:
				return flags[0]
			else:
				rez = flags[0]
				for i in range(1, len(flags)):
					if ops[i-1].strip() == "AND":
						rez = rez and flags[i]
					elif ops[i-1].strip() == "OR":
						rez = rez or flags[i]
					else:
						raise ValueError("Unsupported logical operator in compound filtering expression.")
				return rez
		else:
			raise ValueError("Expression list was malformed.")


class SlimVCF():

	def _split_header_piece(self, x):
		pieces = x.split("=", 1)
		if (len(pieces) < 2):
			raise ValueError("That wasn't a valid VCF header element.")
		else:
			keep = pieces.pop()
			if not keep.startswith("<") and keep.endswith(">"):
				raise ValueError("That wasn't a valid VCF header element.")
			else:
				keep = keep[1:][:-1]
				kv = keep.split(",")
				vals = {}
				for i in kv:
					y = i.split("=")
					if len(y) != 2:
						continue
					else:
						vals.update({ y[0].upper(): y[1] })
				return vals

	def _parse_sample_header(self, x):
		pieces = x.split()
		return pieces[9:]

	def _pairlist_to_dict(self, x):
		d = OrderedDict()
		if x == ".":
			d = {}
		else:
			for e in x.split(";"):
				try:
					k,v = e.split("=")
					d.update({k:v})
				except:
					d.update({e: None})
		return d

	def _parse_variant_entry(self, x):

		ns = len(self.samples)
		pieces = x.split()
		fields = ["chr","pos","id","ref","alt","qual","filter","info","format"] + self.samples

		if len(pieces) != len(fields):
			raise ValueError("Malformatted line:\n" + str(pieces))
		else:
			## parse required fields
			out = {}
			for i in range(0,9):
				out.update({ fields[i]: pieces[i] })
			out["alt"] = out["alt"].split(",")
			out["filter"] = list(out["filter"].split(";"))
			out["format"] = list(out["format"].split(":"))
			out["info"] = self._pairlist_to_dict(out["info"])
			## fix SVTYPE field for files generated by GenomeSTRiP
			if "SVTYPE" in out["info"] and "GSCNCATEGORY" in out["info"]:
				out["info"]["SVTYPE"] = out["info"]["GSCNCATEGORY"]
			## now parse genotypes
			gt_raw = OrderedDict()
			phased = OrderedDict()
			gt = np.zeros((ns, 2), dtype = np.int)
			gt[:,:] = -1
			for i in range(9, len(pieces)):
				gt_pieces = pieces[i].split(":")
				gt_raw.update({ self.samples[i-9]: OrderedDict(zip(out["format"], gt_pieces)) })
				p0 = gt_pieces[0]
				if p0 == ".":
					p0 = "./."
				phased.update({ self.samples[i-9]: ("|" in p0) })
				alleles = re.split(r"[\|\/]", p0)
				for j in range(0, len(alleles)):
					try:
						alleles[j] = int(alleles[j])
					except:
						alleles[j] = -1
				gt[i-9,:] = alleles
			out["geno"] = gt_raw
			out["phased"] = phased
			site = VariantSite(out)
			site.geno_matrix = gt
			return(site)


	def __init__(self, stream = None):

		self.header = []
		self.samples = []
		self.contigs = []
		self.ready = False
		self._stream = stream

		if stream is not None:
			for line in stream:
				line = line.strip()
				if line.startswith("##"):
					if line.upper().startswith("##CONTIG"):
						self.contigs.append( self._split_header_piece(line)["ID"] )
					self.header.append(line)
				elif line.upper().startswith("#CHROM"):
					self.samples = self._parse_sample_header(line)
					self.header.append(line)
					break

			if len(self.samples):
				self.ready = True
			else:
				raise Exception("Couldn't read VCF header.")

	def __iter__(self):

		if not self.ready:
			raise Exception("Haven't read the VCF header yet.")

		while True:
			try:
				line = next(self._stream)
				if not line.startswith("#"):
					yield self._parse_variant_entry(line.strip())
			except StopIteration as e:
				raise e


	def print_header(self, sample_idx = None):
		for i in range(0, len(self.header)-1):
			print(self.header[i], file = sys.stdout)
		last_line = self.header.pop().strip().split()
		if sample_idx is None:
			sample_idx = range(9, len(last_line))
		else:
			sample_idx = [ i+9 for i in sample_idx ]
		to_write = last_line[0:9] + [ last_line[i] for i in sample_idx ]
		print("\t".join(to_write))

	## get positional index of sample column in a VCF, failing safely if sample not found
	def index_of(self, query):
		out = OrderedDict()
		for s in query:
			if s in self.samples:
				out.update({ s: self.samples.index(s) })
		return out
