#!/usr/bin/env python

usage = \
"""
Plots times versus frame number to look for jumps. Rejection is
applied by absolute value first, then by RMS scatter. You can
optionally subtract off the final best-fit line.
"""

# builtins
import argparse, os, re, sys, traceback, subprocess, datetime

# thirdparty
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

# mine
from trm import ultracam

parser = argparse.ArgumentParser(description=usage,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('run',  help='run to analyse')

# optional
parser.add_argument('--absolute', '-a', type=float, default=0.001,
                    help='absolute rejection threshold to eliminate bad times, units of days')
parser.add_argument('--rms', '-r', type=float, default=4., help='RMS rejection threshold to eliminate bad times')
parser.add_argument('--linear', '-l', action='store_true', help='remove linear trend in times or not')
parser.add_argument('--server', '-s', action="store_true",help="get runs from server")

# OK, done with arguments.
args = parser.parse_args()

n   = []
gps = []
tdat  = ultracam.Rtime(args.run,server=args.server)
for nf, time in enumerate(tdat):
    n.append(nf+1)
    gps.append(time[1]['gps'])

n   = np.array(n)
gps = np.array(gps)

print 'Loaded',len(n),'timestamps.'

def func(x, n, gps):
    t0  = x[0]
    cad = x[1]
    return gps-t0-cad*n

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
        print 'Rejected time for frame',n[ok][d.argmax()]
        ok[ind[ok][d.argmax()]] = False
        rej = True
    elif d.max() > args.rms*std:
        print 'Rejected time for frame',n[ok][d.argmax()]
        ok[ind[ok][d.argmax()]] = False
        rej = True

if args.linear:
    gps[ok] = func(x, n[ok], gps[ok])
    plt.ylabel('Time - linear fit (days)')
else:
    plt.ylabel('Time (days)')
plt.plot(n[ok],gps[ok])
plt.plot(n[ok],gps[ok],'g.')
plt.xlabel('Frame number')
plt.xlim(0,len(n)+1)
plt.show()

