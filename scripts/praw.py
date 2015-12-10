#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function

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
parser.add_argument('-n', dest='nccd', type=int, default=0, help='CCD to plot (0 for all)')
parser.add_argument('-plo', type=float, default=2., help='Lower percentile for intensity display')
parser.add_argument('-phi', type=float, default=98., help='Upper percentile for intensity display')
parser.add_argument('-f', dest='first', type=int, default=1, help='first frame to plot (default = 1)')
parser.add_argument('-l', dest='last', type=int, default=0, help='last frame to plot (0 to go up to last one)')
parser.add_argument('-s', dest='sleep', type=float, default=0, help='number of seconds to pause between frames')
parser.add_argument('-r', dest='back', action='store_true', help='remove median background from each window')
parser.add_argument('-b', dest='bias', help='bias frame to subtract (ucm file)')
parser.add_argument('-u', dest='ucam', action='store_true', help='Get data via the ULTRACAM FileServer')
parser.add_argument('-e', dest='every', action='store_true', help='Show every blue frame, not just good ones')
parser.add_argument('-x1', type=float, help='left-hand X-limit')
parser.add_argument('-x2', type=float, help='right-hand X-limit')
parser.add_argument('-y1', type=float, help='lower Y-limit')
parser.add_argument('-y2', type=float, help='upper Y-limit')

# OK, done with arguments.
args = parser.parse_args()

# Check arguments
run   = args.run
if not args.ucam and not os.path.exists(run + '.xml'):
    print('ERROR: could not find',run+'.xml')
    exit(1)
if not args.ucam and not os.path.exists(run + '.dat'):
    print('ERROR: could not find',run+'.dat')
    exit(1)

first = args.first
if first < 0:
    print('ERROR: first frame must be >= 0')
    exit(1)

nccd = args.nccd
if nccd < 0:
    print('ERROR: nccd must be >= 0')
    exit(1)
nccd -= 1

plo = args.plo
if plo < 0. or plo > 100.:
    print('ERROR: plo must lie from 0 to 100')
    exit(1)

phi = args.phi
if phi < 0. or phi > 100.:
    print('ERROR: phi must lie from 0 to 100')
    exit(1)

# more imports
import time, copy
from trm import ultracam
from ppgplot import *

if args.bias:
    bias = ultracam.MCCD.rucm(args.bias)

# Now do something
pgopen('/xs')

fnum  = args.first
first = True
rdat  = ultracam.Rdata(run,args.first,server=args.ucam)

if rdat.instrument == 'ULTRACAM':
    saveBlue = None
    for mccd in rdat:

        if rdat.nblue > 1 and not args.every:
            # save the last read blue frame to avoid
            # plotting the intermediate bias frames
            # when "nblue" is in operation
            if fnum % rdat.nblue == 0:
                saveBlue = mccd[2]
            elif saveBlue:
                if args.bias or args.back:
                    mccd[2] = copy.deepcopy(saveBlue)
                else:
                    mccd[2] = saveBlue

        if args.bias:
            if first and bias != mccd:
                try:
                    bias = bias.cropTo(mccd)
                except ultracam.UltracamError as err:
                    print('UltracamError:',err)
                    print('Bias format:\n',bias.format())
                    print('Data format (should be subset of the bias):\n',mccd.format())
                    exit(1)
            first = False
            mccd -= bias

        if args.back:
            mccd.rback(nccd)

        prange = mccd.plot(plo,phi,nccd,close=False,x1=args.x1,x2=args.x2,y1=args.y1,y2=args.y2)
        print('Plotted',mccd.head.value('Run.run'),'frame',mccd.head.value('Frame.frame'),'plot range(s):',prange)
    
        fnum += 1
        if args.last > 0 and fnum > args.last:
            break
        time.sleep(args.sleep)

elif rdat.instrument == 'ULTRASPEC':

    x1    = 0.5 if args.x1 is None else args.x1
    x2    = rdat.nxmax+0.5 if args.x2 is None else args.x2
    y1    = 0.5 if args.y1 is None else args.y1
    y2    = rdat.nymax+0.5 if args.y2 is None else args.y2
    
    pgenv(x1,x2,y1,y2,1,0)

    for ccd in rdat:

        if args.bias:
            if first and bias != ccd:
                try:
                    bias = bias.cropTo(ccd)
                except ultracam.UltracamError as err:
                    print('UltracamError:',err)
                    print('Bias format:\n',bias.format())
                    print('Data format (should be subset of the bias):\n',ccd.format())
                    exit(1)
            first = False
            ccd -= bias

        if args.back:
            ccd.rback(nccd)

        vmin, vmax = ccd.centile((plo,phi))

        ccd.plot(vmin,vmax)
        print('Plotted',ccd.head.value('Run.run'),'frame',ccd.head.value('Frame.frame'),'plot range:',vmin,'to',vmax)
    
        fnum += 1
        if args.last > 0 and fnum > args.last:
            break
        time.sleep(args.sleep)

pgclos()

