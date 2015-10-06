from __future__ import absolute_import
from __future__ import print_function
from six.moves import range
#!/usr/bin/env python

usage = \
"""
Displays ULTRACAM files with ds9. It leaves ds9 set for use of a ruler
'region' to allow position angle measurement. Multiple Windows of a given CCD
are displayed in a single mosaic. It can display both raw data and ucm files.
If this is any ambiguity (i.e. run.dat, run.xml and run.ucm files are all
present), it will default to the raw date; append '.ucm' to the run to force
display of the ucm file.

Note that the first time you run this, it can be slow to set up the xpa
nameserver the first time around, but should be OK on subsequent calls,
provide you leave ds9 running.

If you don't specify a CCD, all will be displayed and you can click through
with the ds9 tile options, including switching between one frame and all
three.
"""

if __name__ == '__main__':

    # just import these for speed. After arguments are OK-ed, some more imports
    import argparse, os

    parser = argparse.ArgumentParser(description=usage)

    # positional
    parser.add_argument('run', help='run to plot, e.g. "run045"')

    # optional
    parser.add_argument('-n', dest='nccd', type=int, help='CCD to plot')
    parser.add_argument('-f', dest='frame', type=int, default=1, help='frame number to plot, 0 for the last')
    parser.add_argument('-r', dest='back', action='store_true', help='remove median background from each window')
    parser.add_argument('-b', dest='bias', help='bias frame to subtract (ucm file)')
    parser.add_argument('-u', dest='ucam', action='store_true', help='Get data via the ULTRACAM FileServer')
    parser.add_argument('-v', dest='height', type=int, help='vertical height of display in pixels')
    parser.add_argument('-w', dest='width', type=int, help='width of display in pixels')

    # OK, done with arguments.
    args = parser.parse_args()

    # Check arguments
    run   = args.run

    if args.ucam:
        ucm = False
        nframe = args.frame
        if nframe < 0:
            print('ERROR: first frame must be >= 0')
            exit(1)

    elif os.path.exists(run + '.xml') and os.path.exists(run + '.dat'):
        nframe = args.frame
        if nframe < 0:
            print('ERROR: first frame must be >= 0')
            exit(1)
        ucm = False

    elif run.endswith('.ucm') and os.path.exists(run):
        ucm = True

    elif os.path.exists(run + '.ucm'):
        ucm = True
        run += '.ucm'

    else:
        print('ERROR: could not find one or both of',run+'.xml','and',run+'.dat','or',run+'.ucm')
        exit(1)

    if args.nccd:
        nccd = args.nccd
        if nccd < 1:
            print('ERROR: nccd must be > 0')
            exit(1)
        nccd -= 1

    # more imports
    import time, copy, tempfile
    import numpy as np
    import pyfits as fits
    from trm import ultracam

    # Load bias
    if args.bias:
        # load bias frame, if any
        bias = ultracam.MCCD.rucm(args.bias)

    # Load data
    if ucm:
        mccd = ultracam.MCCD.rucm(run)
    else:
        rdat  = ultracam.Rdata(run,nframe,server=args.ucam)
        mccd  = rdat()

    # Subtract bias
    if args.bias:
        mccd = mccd - bias

    # Subtract background
    if args.back:
        ccd.rback()

    # Determine indices of CCDs to plot
    if mccd.head.value('Instrument.instrument') == 'ULTRACAM':
        if args.nccd:
            nccds = [nccd,]
        else:
            nccds = list(range(len(mccd)))
    else:
        nccds = [0,]

    # Set display dimensions, pixels
    if args.height:
        height = args.height
    else:
        height = 500

    if args.width:
        width = args.width
    else:
        width = len(nccds)*height

    # Fire up ds9
    import ds9
    ds = ds9.ds9()
    ds.set('frame delete all')
    ds.set('width %d' % (width))
    ds.set('height %d' % (height))
    for i, nccd in enumerate(nccds):

        if mccd.head.value('Instrument.instrument') == 'ULTRACAM':
            ccd = mccd[nccd]
        else:
            ccd = mccd

        if nccd == 0:
            ccdname = 'Red CCD'
        elif nccd == 1:
            ccdname = 'Green CCD'
        elif nccd == 2:
            ccdname = 'Blue CCD'

        # goto and create frame if needed
        ds.set('frame %d' % (i+1))

        # write out to an iraf mosaic fits file
        fobj = ultracam.ccd2fits(ccd, ccdname)

        # the intensity scale
        ds.set('scale zscale')

        # set the colour map
        ds.set('cmap heat')

        # make the same for each window
        ds.set('scale scope global')

        # load as a mosaic image
        ds.set('mosaicimage iraf ' + fobj.name)

        # draw a box around edge of potential imaging area
        comm = 'box(' + str((1+ccd.nxmax)/2.) + ',' + \
            str((1+ccd.nymax)/2.) + ',' + str(ccd.nxmax) + ',' + \
            str(ccd.nymax) + ') # color=black'

        ds.set('regions', comm)

        # initially set to the centre of the CCD
        ds.set('pan to %d %d physical' % (ccd.nxmax/2,ccd.nymax/2))

        # zoomed out
        zh = ccd[0].xbin*(0.9*width)/ccd.nymax/len(nccds)
        zv = ccd[0].ybin*(0.9*height)/ccd.nymax
        ds.set('zoom %4.2f' % (min(zh,zv)))


    if len(nccds) > 1:
        ds.set('tile frames')
        ds.set('tile column')
        ds.set('frame lock physical')

    # set a ruler
    ds.set('regions shape ruler')



