#!/usr/bin/env python

usage = \
"""
Reads in ULTRACAM raw data files and writes out to 3D FITS files together with
an ASCII file of mid-exposure MJDs and exposure times in seconds. A run called 
'run123' with one CCD and two windows will generate three files called 
'run123_1.fits' and 'run123_2.fits' and 'run123.times'. If there is more than one
CCD, then an extra '_1', '_2' will appear to identify the CCDs.
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
parser.add_argument('-i', dest='interval', type=int, default=100, help='interval for reporting progress')
parser.add_argument('-s', dest='split', action='store_true', help='split files by CCD')

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

# more imports
import numpy as np
import astropy.io.fits as fits
from trm import ultracam

if args.bias:
    bias = ultracam.MCCD.rucm(args.bias)

# Now do something
fnum  = args.first
first = True
flt   = args.bias or args.back
dtype = np.float32 if flt else np.uint16
rdat  = ultracam.Rdata(run,args.first,flt=flt,server=args.ucam)

nccd = args.nccd
if nccd < 0:
    print 'ERROR: nccd must be >= 0'
    exit(1)
elif nccd > rdat.nccd:
    print 'ERROR: the data contains only',rdat.nccd,'CCDs.'
    exit(1)
nccd -= 1

if nccd == -1:
    ccds = range(rdat.nccd)
else:
    ccds = range(0,nccd+1)

# For the 3D option we need to prepare the files in advance
# first compute number of slices
if args.last:
    nz = min(rdat.ntotal(),args.last) - fnum + 1
else:
    nz = rdat.ntotal() - fnum + 1

# now for each window create a FITS file. 
nbytes = 4 if flt else 2
for nc in ccds:
    for nw,rwin in enumerate(rdat.win):
        # a 3D array of the right type
        # nx needs fiddling because of the pixel shift bug
        nx = rwin.nx-1 if rdat.version  == -1 else rwin.nx
        tdim = np.zeros((1,rwin.ny,nx),dtype=dtype)

        # the primary HDU
        phdu  = fits.PrimaryHDU(data=tdim)

        # allocate NBLOCK block of 2880 bytes for the header
        # Each allows 36 cards; one is used up for the END
        head   = phdu.header
        NBLOCK = 1
        while len(head) < NBLOCK*36 - 1:
            head.append()

        # manipulate the dimensions ...
        head['NAXIS1']   = nx 
        head['NAXIS2']   = rwin.ny
        head['NAXIS3']   = nz
        npix = nz*rwin.ny*nx

        head['OBJECT']   = (rdat.target,'Object name')
        head['PI']       = (rdat.pi,'Principal investigator')
        head['ID']       = (rdat.id,'Program ID')
        head['OBSRVRS']  = (rdat.observers,'Observers')
        head['FILTER']   = (rdat.filters,'Filter')
        head['FILTER']   = (rdat.filters,'Filter')
        head['SPEED']    = (rdat.speed,'Readout speed')
        head['DTYPE']    = (rdat.dtype,'Data type')
        head['SLIDE']  = (rdat.slidePos,'Slide position, pixels')

        # generate the file name. Include the CCD number if the instrument
        # has more than 1:
        if rdat.nccd == 1:
            fname = os.path.basename(run) + '_' + str(nw+1) + '.fits'
        else:
            fname = os.path.basename(run) + '_' + str(nc+1) + '_' + str(nw+1) + '.fits'

        # write the header
        head.tofile(fname, clobber=args.clobber)

        # now extend the size of the file to accomodate the data which will be 
        # written later. Do this with a trick from the PyFITS manual: seek
        # past the end of the file and write a 0 in the last byte. We also
        # ensure that we write an exact number of blocks to avoid some warning
        # message from pyfits about truncation
        with open(fname,'rb+') as fobj:
            tbytes = len(head.tostring()) + nbytes*npix
            abytes = 2880 * (tbytes / 2880) + (2880 if tbytes % 2880 else 0)
            fobj.seek(abytes - 1)
            fobj.write('\0')

        print 'Created container file =',fname
        print 'Will now read the data'

        # at this stage we have generated a FITS container file for 
        # each window ready to receive the data

# Now the data reading step. 

for mccd in rdat:

    # kludge to make a CCD generated when there is only one
    # behave like an MCCD
    if isinstance(mccd, ultracam.CCD):
        mccd = [mccd,]

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

    if fnum == args.first:
        # on first pass, open all the files:
        # one FITS per window
        # one ASCII file per CCD.
        hduls = []
        times = []
        for nc in ccds:
            ccd = mccd[nc]
            if len(mccd) == 1:
                fname = os.path.basename(run) + '.times'
            else:
                fname = os.path.basename(run) + '_' + str(nc+1) + '.times'
            times.append(open(fname,'w'))
            times[-1].write(
"""#
# Data written by to3dfits.py
#
# The times in this file represent the MJD (UTC) at the
# centre of each exposure, followed by the exposure time
# in seconds, followed by a flag indicating reliability.
# Typically the first frame or two, or in drift mode the
# first several frames will not have reliable times.
#
""")
            for nw, win in enumerate(ccd):
                if len(mccd) == 1:
                    fname = os.path.basename(run) + '_' + str(nw+1) + '.fits'
                else:
                    fname = os.path.basename(run) + '_' + str(nc+1) + '_' + str(nw+1) + '.fits'
                hduls.append(fits.open(fname,mode='update',memmap=True))

    # on all passes, store data slice by slice
    nhdul = 0
    for nc in ccds:
        ccd = mccd[nc]
        # write times
        times[nc].write('{0:15.9f} {1:11.6f} {2:1d}\n'.format(ccd.time.mjd,ccd.time.expose,ccd.time.good))
        for nw, win in enumerate(ccd):
            hduls[nhdul][0].data[fnum-args.first,::] = win.data
            nhdul += 1

    if fnum % args.interval == 0:
        print 'Stored data for frame',fnum
        
    fnum += 1
    if args.last > 0 and fnum > args.last:
        break

# shut down the files
for hdul in hduls:
    hdul.close()

for time in times:
    time.close()
