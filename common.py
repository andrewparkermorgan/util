#! /usr/bin/env python

##	--- common.py --- ##
##	Date: 21 Feb 2014
##	Updated: 7 Aug 2014
##	Purpose: miscellaneous common utility functions

import os
import sys
import argparse
import csv
import subprocess

## functions for command-line argument validation; intended to work with argparse module
def expand_all(path):
	return os.path.expanduser( os.path.expandvars(path) )

def readable_dir(indir):
	indir = expand_all(indir)
	if not os.path.isdir(indir):
		raise argparse.ArgumentError("readable_dir:{0} is not a valid path".format(indir))
	if os.access(indir, os.R_OK):
		return indir
	else:
		raise argparse.ArgumentError("readable_dir:{0} is not a readable dir".format(indir))

def writeable_dir(indir):
	indir = expand_all(indir)
	if not os.path.isdir(indir):
		raise argparse.ArgumentError("writeable_dir:{0} is not a valid path".format(indir))
	if os.access(indir, os.W_OK):
		return indir
	else:
		raise argparse.ArgumentError("writeable_dir:{0} is not a writeable dir".format(indir))

def readable_file(infile):
	infile = expand_all(infile)
	if not os.path.isfile(infile):
		raise argparse.ArgumentError("readable_file:{0} is not a valid file path".format(infile))
	if os.access(infile, os.R_OK):
		return infile
	else:
		raise argparse.ArgumentError("readable_file:{0} is not a readable file".format(infile))

def writeable_file(infile):
	infile = expand_all(infile)
	if not os.path.isfile(infile):
		raise argparse.ArgumentError("writeable_file:{0} is not a valid file path".format(infile))
	if os.access(infile, os.W_OK):
		return infile
	else:
		raise argparse.ArgumentError("writeable_file:{0} is not a writeable file".format(infile))

def readable_or_stdin(infile):
	if not infile == "-":
		return readable_file(infile)
	else:
		return infile

def readable_or_stdin_handle(infile):
	if not infile == "-":
		return argparse.FileType("rU")(infile)
	else:
		return sys.stdin

def writeable_or_stdout_handle(infile):
	if not infile == "-":
		return argparse.FileType("w")(infile)
	else:
		return sys.stdout

def comma_list(value):
	return value.split(",")

def list_from_file(infile):
	if not (os.path.isfile(infile) and os.access(infile, os.R_OK)):
		raise argparse.ArgumentError("list_from_file:{0} is not a readable file".format(infile))
	else:
		ll = []
		with open(infile, "rU") as ff:
			samples = csv.reader(ff, delimiter = ",")
			for line in samples:
				ll.append(line[0])
		return(ll)

def count_lines(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])
