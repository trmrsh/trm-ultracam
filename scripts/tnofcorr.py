#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function


usage = \
"""
TNO flats suffer a pronounced brightening in the middle of the field (~10%).
On the basis that this is caused by scattered light and is not real this
script attempts to correct for this by dividing out a blurred version
(2D gaussian convolution) of the image. The idea is that this is larger
scale than the features you want to remove (specks of dust, etc), but
smaller scale than the scattered light
"""

import argparse, os
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from trm import ultracam

parser = argparse.ArgumentParser(description=usage)

# positional
parser.add_argument('iflat', help='name of input flat field (full frame)')
parser.add_argument('oflat', help='name of output flat field')

# optional
parser.add_argument('-f', dest='fwhm', type=float, default=40.,
                    help='FWHM of smoothing gaussian')


# OK, done with arguments.
args = parser.parse_args()

# Check arguments
if not os.path.exists(args.iflat + '.ucm'):
    print('ERROR: could not find',args.iflat+'.ucm')
    exit(1)

# Load the flat
iflat = ultracam.MCCD.rucm(args.iflat)

if len(iflat) != 1:
    print('Expected only 1 CCD for ULTRASPEC')
    exit(1)

ccd = iflat[0]
if len(ccd) != 1:
    print('Expected only 1 Window for full frame flat')
    exit(1)

win   = ccd[0]
arr   = win.data
ny,nx = arr.shape

if win.llx != 1 or win.lly != 1 or win.xbin != 1 or win.ybin != 1 or \
        nx != 1056 or ny != 1072:
    print('Only full frame, unbinned TNO data accepted.')
    exit(1)

# good data region X: 18 to 1039, Y: 3 to 1026.
# we replace data outside of this region with the nearest pixel values, first
# horizontally, then vertically.

# left side

arr[:,:17] = arr[:,17:18]

# right side
arr[:,1039:] = arr[:,1038:1039]

# bottom
arr[:2,:] = arr[2:3,:]

# top
arr[1026:,:] = arr[1025:1026,:]

smoothed = gaussian_filter(arr, args.fwhm/2.3548)

arr /= smoothed

iflat.wucm(args.oflat)

print('Corrected file written to',args.oflat)

