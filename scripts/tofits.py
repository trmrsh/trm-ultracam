from __future__ import absolute_import
from __future__ import print_function
from six.moves import range
#!/usr/bin/env python

usage = \
""" 
Reads in ULTRACAM/SPEC raw data files and writes out to FITS files. The
output files can be one per exposure, or in the case of multi-CCD formats, one
per CCD per exposure.  """


# just import these for speed. After arguments are OK-ed, some more imports
import argparse, os

parser = argparse.ArgumentParser(description=usage)

# positional
parser.add_argument('run', help='run to plot, e.g. "run045"')

# optional
parser.add_argument('-n', dest='nccd', type=int, default=0,
                    help='which CCD to process, 0 for the lot')
parser.add_argument('-f', dest='first', type=int, default=1,
                    help='first frame to read (default = 1)')
parser.add_argument('-l', dest='last', type=int, default=0,
                    help='last frame to read (0 to go up to last one)')
parser.add_argument('-r', dest='back', action='store_true',
                    help='remove median background from each window')
parser.add_argument('-b', dest='bias',
                    help='bias frame to subtract (ucm file)')
parser.add_argument('-u', dest='ucam', action='store_true',
                    help='get data via the ULTRACAM FileServer')
parser.add_argument('-c', dest='clobber', action='store_true',
                    help='clobber existing files')
parser.add_argument('-i', dest='interval', type=int, default=100,
                    help='interval for reporting progress')
parser.add_argument('-s', dest='split', action='store_true',
                    help='split files by CCD')

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

# more imports
import os
import numpy as np
try:
    import astropy.io.fits as fits
except:
    import pyfits as fits
from trm import ultracam

if args.bias:
    bias = ultracam.MCCD.rucm(args.bias)

# Now do something
fnum  = args.first
first = True
rdat  = ultracam.Rdata(run,args.first,server=args.ucam)

nccd = args.nccd
if nccd < 0:
    print('ERROR: nccd must be >= 0')
    exit(1)
elif nccd > rdat.nccd:
    print('ERROR: the data contains only',rdat.nccd,'CCDs.')
    exit(1)
nccd -= 1

if nccd == -1:
    ccds = list(range(rdat.nccd))
else:
    ccds = list(range(0,nccd+1))

# Now the data reading step.

def add_head(header, fkey, mccd, hkey, comment=None):
    """
    Checked addition of parameter into a FITS header.

    header  : the FITS header
    fkey    : FITS keyword
    mccd    : the MCCD object
    hkey    : the header key within the MCCD
    comment : comment. If None, the value from the MCCD will be used
    """
    if hkey in mccd.head:
        if comment is None:
            header[fkey] = (mccd.head.value(hkey),mccd.head.comment(hkey))
        else:
            header[fkey] = (mccd.head.value(hkey),comment)

for mccd in rdat:

    if args.bias:
        if first and bias != mccd:
            try:
                bias = bias.cropTo(mccd)
            except ultracam.UltracamError as err:
                print('UltracamError:',err)
                print('Bias format:\n',bias.format())
                print('Data format (should be subset of the bias):\n',\
                    mccd.format())
                exit(1)
        first = False
        mccd -= bias

    if args.back:
        mccd.rback(nccd)

    # write out the data in a FITS file

    if args.split:
        # now tack on an image one for each window of each CCD
        for nc in ccds:
            ccd = mccd[nc]

            # first cook up the primary HDU and its header
            header = fits.Header()
            if args.bias:
                header['BSUB'] = (True,'Was a bias subtracted or not?')
            else:
                header['BSUB'] = (False,'Was a bias subtracted or not?')
            header['NCCD']    = (nc+1,'CCD number')
            add_head(header, 'OBJECT', mccd, 'User.target')
            add_head(header, 'RUNNUM', mccd, 'Run.run')
            add_head(header, 'FILTER', mccd, 'Run.filters')
            add_head(header, 'PI', mccd, 'User.pi')
            add_head(header, 'ID', mccd, 'Programme ID')
            add_head(header, 'OBSRVRS', mccd, 'User.observers')
            add_head(header, 'OUTPUT', mccd, 'Run.output')
            add_head(header, 'SPEED', mccd, 'Run.speed')
            add_head(header, 'DTYPE', mccd, 'User.dtype')
            add_head(header, 'SLIDE', mccd, 'Run.slidepos')
            add_head(header, 'FOCUS', mccd, 'Run.focus')
            add_head(header, 'CCDTEMP', mccd, 'Run.ccdtemp')
            add_head(header, 'FINGTEMP', mccd, 'Run.fingtemp')
            add_head(header, 'FINGPCEN', mccd, 'Run.fingpcen')
            add_head(header, 'RA', mccd, 'Run.RA')
            add_head(header, 'DEC', mccd, 'Run.Dec')
            add_head(header, 'PA', mccd, 'Run.PA')
            add_head(header, 'ENGPA', mccd, 'Run.EngPA')
            add_head(header, 'TRACK', mccd, 'Run.track')
            add_head(header, 'TTFLAG', mccd, 'Run.ttflag')

            header['MJDUTC']  = (mccd[nc].time.mjd,
                                 'MJD(UTC) at centre of exposure')
            header['EXPOSE']  = (mccd[nc].time.expose,'Exposure time, secs')
            header['TIMEOK']  = (mccd[nc].time.good,'Is time reliable?')
            header.add_comment('File created by tofits.py')

            phdu = fits.PrimaryHDU(header=header)
            hdus = [phdu,]

            for nw, win in enumerate(ccd):
                wheader = fits.Header()
                wheader['NWIN']   = (nw+1,'Window number')

                wheader['CTYPE1'] = ('LINEAR', 'Transformation of X scale')
                wheader['CTYPE2'] = ('LINEAR', 'Transformation of Y scale')
                wheader['CUNIT1'] = ('pixels', 'Units of transformed X scale')
                wheader['CUNIT2'] = ('pixels', 'Units of transformed Y scale')

                fnumber = 1. - float(win.llx - 1)/win.xbin
                wheader['CRPIX1'] = (fnumber,
                                     'Pixel equivalent in X of reference point')
                fnumber = 1. - float(win.lly - 1)/win.ybin
                wheader['CRPIX2'] = (fnumber,
                                     'Pixel equivalent in Y of reference point')
                fnumber = 1.
                wheader['CRVAL1'] = (1., 'X value of reference point')
                wheader['CRVAL2'] = (1., 'Y value of reference point')

                wheader['CD1_1'] = (float(win.xbin),'Binning factor in X')
                wheader['CD1_2'] = 0.
                wheader['CD2_1'] = 0.
                wheader['CD2_2'] = (float(win.ybin),'Binning factor in Y')

                ihdu = fits.ImageHDU(win.data, wheader)
                hdus.append(ihdu)
            hdul = fits.HDUList(hdus)

            # create filename and write out data
            fname = os.path.basename(run) + '_' + str(nc+1) + '_' + \
                    str(fnum) + '.fits'
            hdul.writeto(fname, clobber=args.clobber)

    else:

        # first cook up the primary HDU and its header
        header = fits.Header()
        if args.bias:
            header['BSUB'] = (True,'Was a bias subtracted or not?')
        else:
            header['BSUB'] = (False,'Was a bias subtracted or not?')
        if rdat.nccd > 1:
            header['NCCD']   = (nc+1,'CCD number')

        add_head(header, 'OBJECT', mccd, 'User.target')
        add_head(header, 'RUNNUM', mccd, 'Run.run')
        add_head(header, 'FILTER', mccd, 'Run.filters')
        add_head(header, 'OUTPUT', mccd, 'Run.output')
        add_head(header, 'SPEED', mccd, 'Run.speed')
        add_head(header, 'PI', mccd, 'User.pi')
        add_head(header, 'ID', mccd, 'Programme ID')
        add_head(header, 'OBSRVRS', mccd, 'User.observers')
        add_head(header, 'DTYPE', mccd, 'User.dtype')
        add_head(header, 'SLIDE', mccd, 'Run.slidepos')
        add_head(header, 'FOCUS', mccd, 'Run.focus')
        add_head(header, 'CCDTEMP', mccd, 'Run.ccdtemp')
        add_head(header, 'FINGTEMP', mccd, 'Run.fingtemp')
        add_head(header, 'FINGPCEN', mccd, 'Run.fingpcen')
        add_head(header, 'RA', mccd, 'Run.RA')
        add_head(header, 'DEC', mccd, 'Run.Dec')
        add_head(header, 'PA', mccd, 'Run.PA')
        add_head(header, 'ENGPA', mccd, 'Run.EngPA')
        add_head(header, 'TRACK', mccd, 'Run.track')
        add_head(header, 'TTFLAG', mccd, 'Run.ttflag')

        if rdat.nccd == 1:
            header['MJDUTC'] = (mccd[0].time.mjd,
                                'MJD(UTC) at centre of exposure')
            header['EXPOSE'] = (mccd[0].time.expose,'Exposure time, secs')
            header['TIMEOK'] = (mccd[0].time.good,'Is time reliable?')
        header.add_comment('File created by tofits.py')
        phdu = fits.PrimaryHDU(header=header)
        hdus = [phdu,]

        # now tack on an image one for each window of each CCD
        for nc in ccds:
            ccd = mccd[nc]
            for nw, win in enumerate(ccd):
                wheader = fits.Header()
                if rdat.nccd > 1:
                    header['MJDUTC'] = (mccd[nc].time.mjd,
                                        'MJD(UTC) at centre of exposure')
                    header['EXPOSE'] = (mccd[nc].time.expose,
                                        'Exposure time, secs')
                    header['TIMEOK'] = (mccd[nc].time.good,'Is time reliable?')
                wheader['NCCD']   = nc+1
                wheader['NWIN']   = nw+1

                wheader['CTYPE1'] = ('LINEAR', 'Transformation of X scale')
                wheader['CTYPE2'] = ('LINEAR', 'Transformation of Y scale')
                wheader['CUNIT1'] = ('pixels', 'Units of transformed X scale')
                wheader['CUNIT2'] = ('pixels', 'Units of transformed Y scale')

                fnumber = 1. - float(win.llx - 1)/win.xbin
                wheader['CRPIX1'] = (fnumber,
                                     'Pixel equivalent in X of reference point')
                fnumber = 1. - float(win.lly - 1)/win.ybin
                wheader['CRPIX2'] = (fnumber,
                                     'Pixel equivalent in Y of reference point')
                fnumber = 1.
                wheader['CRVAL1'] = (1., 'X value of reference point')
                wheader['CRVAL2'] = (1., 'Y value of reference point')

                wheader['CD1_1'] = (float(win.xbin),'Binning factor in X')
                wheader['CD1_2'] = 0.
                wheader['CD2_1'] = 0.
                wheader['CD2_2'] = (float(win.ybin),'Binning factor in Y')

                ihdu = fits.ImageHDU(win.data, wheader)
                hdus.append(ihdu)
        hdul = fits.HDUList(hdus)

        # create filename and write out data
        fname = os.path.basename(run) + '_' + str(fnum) + '.fits'
        hdul.writeto(fname, clobber=args.clobber)

    print('Saved frame',fnum,'to',fname)

    fnum += 1
    if args.last > 0 and fnum > args.last:
        break
