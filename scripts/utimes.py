#!/usr/bin/env python

usage = \
"""
Prints out times of an ULTRACAM run.  The script specifically searches for all
directories of the form YYYY-MM-DD, but you can specify a regular expression
to narrow the search. e.g. '2010' will find all runs in directories containing
2010, as well as having the YYYY-MM-DD format.
"""

# builtins
import argparse, os

# mine
from trm import ultracam

parser = argparse.ArgumentParser(description=usage)

# positional
parser.add_argument('run', help='run to list times of')

# optional
parser.add_argument('-s', dest='suppress', action='store_true', help='suppress all but bad output')

# OK, done with arguments.
args = parser.parse_args()

run = args.run
if not os.path.isfile(run + '.xml') or not os.path.isfile(run + '.dat'):
    print 'One or both of',run+'.xml','and',run+'.dat','does not exist.'
    exit(1)

tdat = ultracam.Rtime(run)

first = True
for nf, time in enumerate(tdat):
    if not args.suppress or (first and not time[0].good):
        print 'Frame %d, mid-time = %s, GPS = %s, exposure = %7.4f, status = %s' % \
            (nf+1,ultracam.mjd2str(time[0].mjd,True),ultracam.mjd2str(time[3]['gps'],True),\
                 time[0].expose, 'T' if time[0].good else 'F'),
        if not time[0].good:
            print ', reason =',time[0].reason
        else:
            print
        first = False
