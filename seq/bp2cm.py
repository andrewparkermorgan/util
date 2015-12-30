#! /usr/bin/env python

## --- bp2cm.py --- ##
##	Date: 17 Nov 2014
##	Last: NA
##	Purpose: Assign genetic positions to loci for which only physical position is avaialable, using linear interpolation on some reference map.

import os
import sys
import bisect
import collections
import csv
import argparse
import re

#sys.path.insert(0, os.path.expanduser("~/lib/util/"))
#from myio import common

parser = argparse.ArgumentParser(description = "Convert phsyical to genetic positions using a reference map using linear interpolation.")
parser.add_argument(	"-m","--map", type = argparse.FileType("rU"),
						help = "reference genetic-to-physical (bp-to-cM) map" )
parser.add_argument(	"--delim",
						default = ",",
						help = "delimiter character in the reference map" )
parser.add_argument(	"-o", "--offset", type = int,
						default = 3000000,
						help = "bp position of origin on genetic map [default: %(default)d]" )
parser.add_argument(	"positions", nargs = "*", type = argparse.FileType("rU"),
						default = sys.stdin,
						help = "list of positions to interpolate" )
args = parser.parse_args()

def interpolate_between(a, x1, x2, y1, y2):
	m = float(y2-y1)/(x2-x1)
	dx = a - x1
	return float(dx)*m + y1

rawmap = collections.defaultdict(list)
bp = collections.defaultdict(list)
cm = collections.defaultdict(list)

sys.stderr.write("--- Reading reference map ... ---\n")
mapfile = csv.reader(args.map, delimiter = args.delim)
for row in mapfile:
	if not row[0].startswith("#"):
		try:
			_chrom, _bp, _cm = row[0], int(row[1]), float(row[2])
			if not _chrom.startswith("chr"):
				_chrom = "chr" + str(_chrom)
				if _cm > 0.0:
					rawmap[_chrom].append( (_cm, _bp) )
		except Exception as e:
			pass
			#sys.stderr.write("SKIPPING: {} {} {}\n".format(row[0], row[1], row[2]))
for (k,v) in rawmap.iteritems():
	(cm_pos, bp_pos) = ( list(x) for x in zip(*v) )
	cm_pos.append(0.0)
	bp_pos.append(args.offset)
	cm_pos.sort()
	bp_pos.sort()
	bp[k] = bp_pos
	cm[k] = cm_pos

chroms = sorted(bp.keys())
sys.stderr.write("\nReference map summary:\n")
sys.stderr.write("\tchr\tmin bp\tmax bp\tmin cM\tmax cM\n")
sys.stderr.write("\t---\t---\t---\t---\t---\n")
for c in chroms:
	sys.stderr.write( "\t{}\t{}\t{}\t{}\t{}\n".format(c, bp[c][0], bp[c][-1], cm[c][0], cm[c][-1]) )

#sys.exit(0)
i = 0
sys.stderr.write("\n--- Interpolating positions ... ---\n")
posfile = csv.reader(args.positions, delimiter = args.delim)
for row in posfile:

	chrm = row[0].strip()

	if not chrm in bp:
		#chrm = re.sub("^chr","", chrm)
		if not chrm in bp:
			sys.stderr.write("Warning: {} is not in the map\n".format(row[0]))
			continue

	pos = int(row[3])
	j = bisect.bisect(bp[chrm], pos)
	interp_cm = 0
	if j == 0:
		#if cm[ row[0] ][0] > 0:
		#	interp_cm = interpolate_between(pos, 0, bp[ row[0] ][0], 0.0, cm[ row[0] ][0])
		sys.stderr.write("Warning: bp position {} not interpolatable.\n".format(int(pos)))
	elif j == len(bp[chrm]):
		interp_cm = cm[chrm][-1]
	else:
		(x1, x2) = bp[chrm][ (j-1):(j+1) ]
		(y1, y2) = cm[chrm][ (j-1):(j+1) ]
		interp_cm = interpolate_between(pos, x1, x2, y1, y2)

	print row[0], row[1], interp_cm, int(pos)

	i += 1
	if not i % 10000:
		sys.stderr.write("\t{} positions read\n".format(i))
