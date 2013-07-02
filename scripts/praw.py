#!/usr/bin/env python

usage = \
"""
Plots frames from a raw ULTRACAM file sequentially.
"""

# just import these for speed. After arguments are OK-ed, some more imports
import argparse, os

parser = argparse.ArgumentParser(description=usage)

# positional
parser.add_argument('run', help='run to plot, e.g. "run045"')

# optional
parser.add_argument('--nccd', '-n', type=int, default=0, help='CCD to plot (0 for all)')
parser.add_argument('-plo', type=float, default=2., help='Lower percentile for intensity display')
parser.add_argument('-phi', type=float, default=98., help='Upper percentile for intensity display')
parser.add_argument('--first', '-f', type=int, default=1, help='first frame to plot (default = 1)')
parser.add_argument('--last', '-l', type=int, default=0, help='last frame to plot (0 to go up to last one)')
parser.add_argument('--sleep', '-s', type=float, default=0, help='number of seconds to pause between frames')
parser.add_argument('--back', '-b', action='store_true', help='subtract median background from each window')
parser.add_argument('-x1', type=float, help='left-hand X-limit')
parser.add_argument('-x2', type=float, help='right-hand X-limit')
parser.add_argument('-y1', type=float, help='lower Y-limit')
parser.add_argument('-y2', type=float, help='upper Y-limit')

# OK, done with arguments.
args = parser.parse_args()

# Check arguments
run   = args.run
if not os.path.exists(run + '.xml'):
    print 'ERROR: could not find',run+'.xml'
    exit(1)
if not os.path.exists(run + '.dat'):
    print 'ERROR: could not find',run+'.dat'
    exit(1)

first = args.first
if first < 0:
    print 'ERROR: first frame must be >= 0'
    exit(1)

nccd = args.nccd
if nccd < 0:
    print 'ERROR: nccd must be >= 0'
    exit(1)
nccd -= 1

plo = args.plo
if plo < 0. or plo > 100.:
    print 'ERROR: plo must lie from 0 to 100'
    exit(1)

phi = args.phi
if phi < 0. or phi > 100.:
    print 'ERROR: phi must lie from 0 to 100'
    exit(1)

# more imports
import time
from trm import ultracam
from ppgplot import *

# Now do something
pgopen('/xs')

fnum = args.first
for mccd in ultracam.Rdata(run,args.first):

    if args.back:
        mccd.rback(nccd)

    prange = mccd.plot(plo,phi,nccd,close=False,x1=args.x1,x2=args.x2,y1=args.y1,y2=args.y2)
    print 'Plotted',mccd.head.value('Data.run'),'frame',mccd.head.value('Data.frame'),'plot range(s):',prange
    
    fnum += 1
    if args.last > 0 and fnum > args.last:
        break
    time.sleep(args.sleep)

pgclos()

