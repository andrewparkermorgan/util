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

sys.path.insert(0, os.path.expanduser("~/lib/util/"))
from io import common

parser = argparse.ArgumentParser(description = "Convert phsyical to genetic positions using a reference map using linear interpolation.")
parser.add_argument(	"-m","--map", type = argparse.FileType("rU"),
			help = "reference genetic-to-physical (bp-to-cM) map" )
parser.add_argument(	"--delim",
			default = ",",
			help = "delimiter character in the reference map" )
parser.add_argument(	"positions", type = common.readable_or_stdin_handle,
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
		rawmap[row[0]].append( (float(row[2]), float(row[1])) )
for (k,v) in rawmap.iteritems():
	(cm_pos, bp_pos) = ( list(x) for x in zip(*v) )
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

	if not row[0] in bp:
		sys.stderr.write("Warning: {} is not in the map\n".format(row[0]))
	else:
		pos = float(row[3])
		j = bisect.bisect(bp[ row[0] ], pos)
		interp_cm = 0
		if j == 0:
			if cm[ row[0] ][0] > 0:
				interp_cm = interpolate_between(pos, 0, bp[ row[0] ][0], 0.0, cm[ row[0] ][0])
		elif j == len(bp[ row[0] ]):
			interp_cm = cm[ row[0] ][-1]
		else:
			(x1, x2) = bp[ row[0] ][ (j-1):(j+1) ]
			(y1, y2) = cm[ row[0] ][ (j-1):(j+1) ]
			interp_cm = interpolate_between(pos, x1, x2, y1, y2)
		print row[0], row[1], interp_cm, int(pos)

	i += 1
	if not i % 10000:
		sys.stderr.write("\t{} positions read\n".format(i))
