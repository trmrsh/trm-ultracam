#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function


usage = \
"""
Checks all runs in the present working directory for timing problems. It does this by stripping
out GPS time stamps, fitting a straight line to them and reporting if there any time not at the
extreme ends that are further than expected off the fit.
"""

import argparse, os, re, sys
import numpy as np
from scipy.optimize import leastsq
from trm import ultracam

parser = argparse.ArgumentParser(description=usage)
parser.add_argument('-s','--server',action="store_true",
                    help="get runs from server")
parser.add_argument('-v','--verbose',action="store_true",
                    help="print more information")
parser.add_argument('--recursive','-r',action="store_true",
                    help="run recursively")
parser.add_argument('-b','--bad',action="store_true",
                    help="report only identifiably bad frames")
parser.add_argument('--absolute', '-a', type=float, default=0.0001,
                    help='absolute rejection threshold to eliminate bad times, units of days')
parser.add_argument('--rms', type=float, default=8.,
                    help='RMS rejection threshold to eliminate bad times')
args = parser.parse_args()

if not args.verbose:
    import warnings
    warnings.filterwarnings("ignore")

server = args.server

rmat = re.compile('^run\d\d\d\.xml$')

def func(x, n, gps):
    t0  = x[0]
    cad = x[1]
    return gps-t0-cad*n

# Compile a list of runs

if args.recursive:
    runs = []
    for rpath, rnames, fnames in os.walk('.'):
        runs.extend([os.path.join(rpath,fname[:-4]) for fname in fnames \
                 if rmat.match(fname)] )
else:
    runs = [fname[:-4] for fname in os.listdir('.') \
            if rmat.match(fname)]
runs.sort()

for run in runs:
    try:
        tdat  = ultracam.Rtime(run,server=server)
        if args.verbose: print('Starting on run',run)
        n   = []
        gps = []
        for nf, time in enumerate(tdat):
            n.append(nf+1)
            gps.append(time[1]['gps'])

        n   = np.array(n)
        gps = np.array(gps)

        if len(gps) > 2:
            rej  = True
            x    = [gps.mean(),0.0001]
            ok   = n > -1
            ind  = np.arange(len(ok),dtype=int)
            while rej:
                rej = False
                res   = leastsq(func, x, (n[ok],gps[ok]))
                x     = res[0]
                d     = func(x,n[ok],gps[ok])
                std   = d.std()
                d     = np.abs(d)
                if d.max() > args.absolute:
                    if args.verbose: print('Rejected time for frame',n[ok][d.argmax()])
                    ok[ind[ok][d.argmax()]] = False
                    rej = True
                elif d.max() > args.rms*std:
                    if args.verbose: print('Rejected time for frame',n[ok][d.argmax()])
                    ok[ind[ok][d.argmax()]] = False
                    rej = True

            bstart = 0
            for val in ok:
                if not val:
                    bstart += 1
                else:
                    break

            bend = len(ok)
            for val in ok[::-1]:
                if not val:
                    bend -= 1
                else:
                    break

            if bstart < bend:
                nbad = len(ok[bstart:bend][~ok[bstart:bend]])
                if nbad:
                    print(run,'has',nbad,'bad times')
                elif not args.bad:
                    print(run,'has',nbad,'bad times')
            elif not args.bad:
                print(run,'has 0 bad times')

        elif not args.bad:
            print(run,'has',len(gps),'times total (< 3)')

    except Exception as err:
        if args.verbose:
            print(run,'could not be read (probably a power on)')
            print(err)
