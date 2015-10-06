from __future__ import absolute_import
from __future__ import print_function
#!/usr/bin/env python

usage = \
"""
Checks the present working directory for runs containing the
ULTRASPEC timing problem spotted in January 2014.
"""

import argparse, os, re, sys
from trm import ultracam

parser = argparse.ArgumentParser(description=usage)
parser.add_argument('-s','--server',action="store_true",
                    help="get runs from server")
args = parser.parse_args()

server = args.server

rmat = re.compile('^run\d\d\d\.xml$')

for rpath, rnames, fnames in os.walk('.'):
    runs  = [fname[:-4] for fname in fnames if rmat.match(fname)]
    runs.sort()
    for run in runs:
        try:
            tdat  = ultracam.Rtime(run,server=server)
            badFrames = []
            for nf, time in enumerate(tdat):
                if nf:
                    if (not time[0].good) and \
                       (time[0].reason == 'timestamp too early'):
                        badFrames.append(nf+1)

            if len(badFrames) == 1:
                print(run,'has timing bug in frame',badFrames[0])
            elif len(badFrames) > 1:
                print(run,'has timing bug in multiple frames:',badFrames)
                print('       please e-mail Tom Marsh at Warwick about this!!')
            else:
                print(run,'is OK')
        except Exception as err:
            print(run,'could not be read (probably a power on)')
            print(err)
