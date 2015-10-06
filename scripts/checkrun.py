#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function


usage = \
"""
Checks the data from a raw ULTRACAM file sequentially. This is really to test the
method 'checkData' more than anything, but could be useful otherwise.
"""

# just import these for speed. After arguments are OK-ed, some more imports
import argparse, os

parser = argparse.ArgumentParser(description=usage)

# positional
parser.add_argument('run', help='run to plot, e.g. "run045"')

# optional
parser.add_argument('-f', dest='first', type=int, default=1, help='first frame to check (default = 1)')
parser.add_argument('-l', dest='last', type=int, default=0, help='last frame to check (0 to go up to last one)')

# OK, done with arguments.
args = parser.parse_args()

# Check arguments
run   = args.run
if not os.path.exists(run + '.xml'):
    print('ERROR: could not find',run+'.xml')
    exit(1)
if not os.path.exists(run + '.dat'):
    print('ERROR: could not find',run+'.dat')
    exit(1)

first = args.first
if first < 0:
    print('ERROR: first frame must be >= 0')
    exit(1)

# more imports
from trm import ultracam

fnum = args.first
rdat = ultracam.Rdata(run,args.first,False)

for mccd in rdat:
    red,green,blue = mccd.checkData()
    mess = ''
    if red[0]:   mess += ' r: ' + red[1] + '.'
    if green[0]: mess += ' g: ' + green[1] + '.'
    if blue[0]:  mess += ' b: ' + blue[1] + '.'
    if red[0] or green[0] or blue[0]:
        print('Frame',fnum,'has possible problems.'+mess)
    fnum += 1

