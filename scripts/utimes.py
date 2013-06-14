#!/usr/bin/env python

usage = \
"""
Computes and and stores stats on all ultracam files it can locate from the directory tree
starting with the current directory. For each run it locates,  it prints out the times of all frames.

The script specifically searches for all directories of the form YYYY-MM-DD, but you can specify
a regular expression to narrow the search. e.g. '2010' will find all runs in directories containing
2010, as well as having the YYYY-MM-DD format.
"""

# builtins
import argparse, os, re, sys, traceback

# thirdparty
import numpy as np

# mine
from trm import ultracam

parser = argparse.ArgumentParser(description=usage)

# optional
parser.add_argument('--regex', '-r', help='regular expression for matching specific directories')
parser.add_argument('--suppress', '-s', action='store_true', help='suppress all but bad output')

# OK, done with arguments.
args = parser.parse_args()

if args.regex is not None:
    rext = re.compile(args.regex)
else:
    rext = None

rdir = re.compile('\d\d\d\d-\d\d-\d\d$')
rmat = re.compile('^run\d\d\d\.xml$')

for rpath, rnames, fnames in os.walk('.'):

    # search only YYYY-MM-DD directories
    if rdir.search(rpath) and (rext is None or rext.search(rpath)):
        
        fnames.sort()
        runs  = [os.path.join(rpath, fname[:-4]) for fname in fnames if rmat.match(fname)]

        print '\n\nFound',len(runs),'runs in directory = ',rpath,'\n'

        for run in runs:

            if not os.path.exists(run + '.dat'):
                print run + '.dat does not exist.'
                continue

            try:
                print 'Processing',run
                tdat = ultracam.Rtime(run)

                first = True
                for time in tdat:
                    if not args.suppress or (first and not time[0].good):
                        print time
                        first = False
                    
            except Exception, err:
                print 'Encountered problem on',run
                traceback.print_exc(file=sys.stdout)
