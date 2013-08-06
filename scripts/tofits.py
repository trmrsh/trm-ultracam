#!/usr/bin/env python

usage = \
"""
Reads in ULTRACAM raw data files and writes out to FITS files.
"""

# just import these for speed. After arguments are OK-ed, some more imports
import argparse, os

parser = argparse.ArgumentParser(description=usage)

# positional
parser.add_argument('run', help='run to plot, e.g. "run045"')

# optional
parser.add_argument('-n', dest='nccd', type=int, default=0, help='which CCD to process, 0 for the lot')
parser.add_argument('-f', dest='first', type=int, default=1, help='first frame to read (default = 1)')
parser.add_argument('-l', dest='last', type=int, default=0, help='last frame to read (0 to go up to last one)')
parser.add_argument('-r', dest='back', action='store_true', help='remove median background from each window')
parser.add_argument('-b', dest='bias', help='bias frame to subtract (ucm file)')
parser.add_argument('-u', dest='ucam', action='store_true', help='get data via the ULTRACAM FileServer')
parser.add_argument('-c', dest='clobber', action='store_true', help='clobber existing files')

# OK, done with arguments.
args = parser.parse_args()

# Check arguments
run   = args.run
if not args.ucam and not os.path.exists(run + '.xml'):
    print 'ERROR: could not find',run+'.xml'
    exit(1)
if not args.ucam and not os.path.exists(run + '.dat'):
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

# more imports
import pyfits as fits
from trm import ultracam

if args.bias:
    bias = ultracam.MCCD.rucm(args.bias)

# Now do something

fnum  = args.first
first = True
rdat  = ultracam.Rdata(run,args.first,server=args.ucam)

for mccd in rdat:

    if args.bias:
        if first and bias != mccd:
            try:
                bias = bias.cropTo(mccd)
            except ultracam.UltracamError, err:
                print 'UltracamError:',err
                print 'Bias format:\n',bias.format()
                print 'Data format (should be subset of the bias):\n',mccd.format()
                exit(1)
        first = False
        mccd -= bias

    if args.back:
        mccd.rback(nccd)

    # write out the data in a FITS file

    # first cook up the primary HDU and its header
    header = fits.Header()
    if args.bias:
        header['BSUB'] = (True,'Was a bias subtracted or not?')
    else:
        header['BSUB'] = (False,'Was a bias subtracted or not?')
    header.add_comment('File created by tofits.py')
    phdu = fits.PrimaryHDU(header=header)
    hdus = [phdu,]

    # now tack on an image one for each window of each CCD
    for nc, ccd in enumerate(mccd):
        for nw, win in enumerate(ccd):
            wheader = fits.Header()
            wheader['NCCD'] = nc+1
            wheader['NWIN'] = nw+1
            ihdu = fits.ImageHDU(win.data, wheader)
            hdus.append(ihdu)
    hdul = fits.HDUList(hdus)

    # create filename and write
    fname = os.path.basename(run) + '_' + str(fnum) + '.fits'
    hdul.writeto(fname, clobber=args.clobber)

    print 'Saved frame',fnum,'to',fname

    fnum += 1
    if args.last > 0 and fnum > args.last:
        break
    
