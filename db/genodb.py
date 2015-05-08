#! /usr/bin/env python -u

## --- genodb.py --- ##
##	Date: 12 Jan 2015
##	Purpose: query genotype databases and dump result to file for reading in from R, because RSQLite interface is really slow

import os
import sys
import sqlite3
import csv
import argparse

parser = argparse.ArgumentParser(description = "Python interface to MegaMUGA and GigaMUGA genotype databases, because RSQLite is so slow.")
parser.add_argument(	"-d", "--db",
						default = "/Volumes/apm_passport/db/arrays/megamuga/megamuga.db",
 						help = "path to genotype database [default: %(default)s]" )
parser.add_argument(	"-o", "--operation", choices = ["samples","genotypes","intensities"],
 						default = "samples",
						help = "operation to perform (ie. what to select)" )
parser.add_argument(	"-s", "--samples", type = argparse.FileType("r"),
						required = False,
						help = "file containing sample IDs" )
parser.add_argument(	"-m", "--markers", type = argparse.FileType("r"),
						required = False,
						help = "file containing marker IDs" )
parser.add_argument(	"--out",
						help = "where to dump output" )
parser.add_argument(	"--by", choices = ["id","name"],
						required = False,
						help = "which chromosome" )
parser.add_argument(	"--exact", action = "store_true",
						default = False,
						help = "perform exact matching of ids" )
parser.add_argument(	"--mda", action = "store_true",
						default = False,
						help = "assume the schema of the MDA database from Jax" )
parser.add_argument(	"--chr", type = str,
						required = False,
						help = "which chromosome" )
parser.add_argument(	"--from-bp", type = int,
						required = False,
						help = "start position" )
parser.add_argument(	"--to-bp", type = int,
						required = False,
						help = "end position" )
parser.add_argument(	"--pos",
						default = "position",
						help = "name of column with chromosomal position [default: %(default)s]" )
args = parser.parse_args()

def read_ids(ff):

	ids = []
	ifile = csv.reader(ff)
	for row in ifile:
		if not row[0].startswith("#"):
			ids.append(row[0])
	return ids

def insert_samples(conn, ids, coltype = None, colname = None, quote_char = "'", table = "_mysamples"):

	## attempt to guess at data type of IDs if none provided
	if coltype is None and colname is None:
		quote_char = ""
		colname = "id"
		coltype = "int"
		try:
			int(ids[0])
		except ValueError as e:
			quote_char = "'"
			colname = "name"
			coltype = "text"

	#sys.stderr.write("{} {} {}\n".format(ids[0], colname, coltype))

	sql = "CREATE TEMP TABLE {} ({} {});".format(table, colname, coltype)
	sys.stderr.write(sql + "\n")
	c = conn.cursor()
	c.execute(sql)

	for i in ids:
		sql = "INSERT INTO {} VALUES ({}{}{});".format(table, quote_char, i, quote_char)
		#sys.stderr.write(sql + "\n")
		c.execute(sql)

	conn.commit()
	return colname


try:
	conn = sqlite3.connect(args.db)
except Exception as e:
	print e
	sys.exit(1)

if args.operation == "samples":

	ids = read_ids(args.samples)
	if args.exact:
		cols = [ "s."+x for x in ["id", "name", "well", "batch", "sex", "flags", "timestamp"] ]
		insert_samples(conn, ids)
		sql = "SELECT {} FROM samples s INNER JOIN _mysamples as sg ON s.{} = sg.{};".format(cols, args.by, args.by)
	else:
		sql = "SELECT {} FROM samples WHERE s.name LIKE '{}';".format(ids[0])

	sys.stderr.write(sql + "\n")
	c = conn.cursor()
	c.execute(sql)

	sys.stderr.write("dumping output to: {}\n".format(args.out))

	with open(args.out, "w") as ff:
		outfile = csv.writer(ff)
		for row in c.fetch():
			outfile.writerow(row)

elif args.operation == "intensities":

	ids = read_ids(args.samples)
	insert_samples(conn, ids)

	fields = []
	sql = ""

	if not args.mda:

		fields = ["marker","sid","id","chr","pos","x","y","call"]
		sql =	"SELECT m.name as marker, s.id as sid, s.name as id, m.chromosome as chr, \
				m.{} as pos, g.x, g.y, g.call \
				FROM samples as s \
				INNER JOIN _mysamples as sg ON s.id = sg.id \
				INNER JOIN genotypes as g ON g.sampleID = s.id \
				INNER JOIN snps as m on g.snpID = m.id ".format(args.pos)

		if args.markers:
			markers = read_ids(args.markers)
			colname = insert_samples(conn, markers, table = "_mymarkers")
			sql = sql + "INNER JOIN _mymarkers as mym ON m.{} = mym.{} ".format(colname, colname)

		if args.chr:
			sql = sql + "WHERE chr = '{}' ".format(args.chr)

		if args.from_bp:
			sql = sql + "AND pos >= {} ".format(args.from_bp)

		if args.to_bp:
			sql = sql + "AND pos <= {} ".format(args.to_bp)

	else:

		fields = ["sid","name","marker","chr","pos","cM","mm10pos","flag"]
		sql =	"SELECT s.id as sid, s.name as name, m.snpID as marker, m.chrID as chr, \
				m.positionBp as pos, m.positioncM as cM, m.mm10positionBp as mm10pos, \
				m.flagged as flag, {} \
				FROM samples as s \
				INNER JOIN _mysamples as sg ON s.id = sg.id \
				{} "

		if args.operation == "intensities":
			fields = fields + ["avg","contrast"]
			sql = sql.format("i.average as avg, i.contrast as contr", "INNER JOIN intensity as i on i.sampleID = s.id \
																	   INNER JOIN snpInfo as m on i.snpID = m.snpID")
		else:
			fields = fields + ["alleleA","alleleB","call"]
			sql = sql.format("m.alleleA, m.alleleB, g.call", "INNER JOIN genotypes as g ON g.sampleID = s.id \
															  INNER JOIN snpInfo as m on g.snpID = m.snpID ")

		if args.markers:
			markers = read_ids(args.markers)
			insert_samples(conn, markers, "_mymarkers")
			sql = sql + "INNER JOIN _mymarkers as mym ON m.{} = mym.{} ".format(args.by, args.by)

		if args.chr:
			sql = sql + "WHERE m.chrID = '{}' ".format(args.chr)

		if args.from_bp:
			sql = sql + "AND m.positionBp >= {} ".format(args.from_bp)

		if args.to_bp:
			sql = sql + "AND m.positionBp <= {} ".format(args.to_bp)

	sys.stderr.write(sql + "\n")
	c = conn.cursor()

	sys.stderr.write("dumping output to: {}\n".format(args.out))

	i = 0
	with open(args.out, "w") as ff:
		outfile = csv.writer(ff)
		outfile.writerow(fields)
		for row in c.execute(sql):
			i += 1
			if i > 1 and not i % 1000:
				sys.stderr.write("... {} records ...\n".format(i))
			outfile.writerow(row)

elif args.operation == "genotypes":

	ids = read_ids(args.samples)
	insert_samples(conn, ids, coltype = "varchar", colname = "id")

	fields = []
	sql = ""

	if args.mda:

		fields = ["sid","name","marker","chr","pos","flag","call"]
		#sql =	"SELECT c.celid as sid, c.name as name, m.jaxid as marker, p.chromosome as chr, p.offset as pos, m.active as flag, g.value as call \n" +
		sql = 	["SELECT c.celid as sid, c.name as name, m.jaxid as marker, m.active as flag, g.value as call ",
				"FROM samples as s ",
				"INNER JOIN samplesinbatch as sb on s.id = sb.sampleid ",
				"INNER JOIN sample_info as c ON s.celid = c.celid ",
				"INNER JOIN callvalues as g ON (g.snpid = m.jaxid AND g.sbid = sb.id) ",
				"INNER JOIN snps as m ON m.jaxid = g.snpid",
				#"INNER JOIN positions as p ON p.snpid = m.jaxid",
				"INNER JOIN _mysamples AS sg ON c.celid = sg.id",
				"LIMIT 10;"]
		sql = "\n".join(sql)
				#"WHERE p.build = {} ".format(args.pos)

		# if args.chr:
		# 	sql = sql + "AND\np.chromosome = '{}' ".format(args.chr)

		# if args.from_bp:
		# 	sql = sql + "AND\np.offset >= {} ".format(args.from_bp)

		# if args.to_bp:
		# 	sql = sql + "AND\np.offset <= {} ".format(args.to_bp)

		sys.stderr.write(sql + "\n")
		c = conn.cursor()

		sys.stderr.write("dumping output to: {}\n".format(args.out))

		i = 0
		with open(args.out, "w") as ff:
			outfile = csv.writer(ff)
			outfile.writerow(fields)
			for row in c.execute(sql):
				i += 1
				if i > 1 and not i % 1000:
					sys.stderr.write("... {} records ...\n".format(i))
				outfile.writerow(row)

	else:
		raise NotImplementedError("To get genotypes from MegaMUGA or GigaMUGA databases, use operation = 'intensities' instead.")

conn.close()
