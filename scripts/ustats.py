#!/usr/bin/env python

usage = \
"""
Computes and and stores stats on all ultracam files it can locate. For each run it locates, it computes
the min, max, median and mean of each CCD of each exposure. The results are stored in a file called 
run###_stats.fits

The script can only be run from the top-level directory containing the raw_data and meta_data directories
because it will look through the raw directories but store the results in the corresponding meta_data
directories. It specifically searches for all directories of the form YYYY-MM-DD, but you can specify
a regular expression to narrow the search. e.g. '2010' will find all runs in directories containing
2010, as well as having the YYYY-MM-DD format.

This script takes a very long time to complete on the full ULTRACAM archive.
"""

# builtins
import argparse, os, re, sys, traceback

# thirdparty
import pyfits, numpy as np

# mine
from trm import ultracam

parser = argparse.ArgumentParser(description=usage)

# optional
parser.add_argument('--regex', '-r', help='regular expression for matching specific directories')
parser.add_argument('--overwrite', '-o', action='store_true', help='overwrite existing statistics files or not')

# OK, done with arguments.
args = parser.parse_args()

if args.regex is not None:
    rext = re.compile(args.regex)
else:
    rext = None

rdir = re.compile('\d\d\d\d-\d\d-\d\d$')
rmat = re.compile('^run\d\d\d\.xml$')

raw  = 'raw_data'
meta = 'meta_data'
if not os.path.isdir(raw) or not os.path.isdir(meta):
    print 'One or both of',raw,'and',meta,'does not exist or is not a directory.'
    print 'Are you running this script from the right directory?'
    exit(1)

for rpath, rnames, fnames in os.walk(raw):

    # search only YYYY-MM-DD directories
    if rdir.search(rpath) and (rext is None or rext.search(rpath)):
        # check for equivalent directory in derived_data
        dpath = meta + rpath[len(raw):]
        if not os.path.exists(dpath):
            print 'Directory',rpath,'has no corresponding',dpath,'and will be skipped.'
            break
        
        fnames.sort()
        runs  = [os.path.join(rpath, fname[:-4]) for fname in fnames if rmat.match(fname)]
        stats = [os.path.join(dpath, fname[:-4] + '_stats.fits') for fname in fnames if rmat.match(fname)]

        print '\n\nFound',len(runs),'runs in directory = ',rpath,'\n'

        for run, stat in zip(runs, stats):

            if not os.path.exists(run + '.dat'):
                print run + '.dat does not exist.'
                continue

            if not args.overwrite and os.path.exists(stat):
                print stat,'exists and will not be overwritten.'
                continue

            try:
                print 'Processing',run
                rdat = ultracam.Rdata(run,flt=False)

                mins,maxs,means,medians = [],[],[], []
                for mccd in rdat:
                    mins.append(mccd.min())
                    maxs.append(mccd.max())
                    means.append(mccd.mean())
                    medians.append(mccd.median())

                if len(mins):
                    # only write out a file if some values were found
                    mins    = np.array(mins)
                    maxs    = np.array(maxs)
                    means   = np.array(means)
                    medians = np.array(medians)

                    # Create FITS table columns. Group stats on each CCD together
                    c1  = pyfits.Column(name='min1',   format='E', unit='DN', array=mins[:,0])
                    c2  = pyfits.Column(name='max1',   format='E', unit='DN', array=maxs[:,0])
                    c3  = pyfits.Column(name='mean1',  format='E', unit='DN', array=means[:,0])
                    c4  = pyfits.Column(name='median1',format='E', unit='DN', array=medians[:,0])
                    
                    c5  = pyfits.Column(name='min2',   format='E', unit='DN', array=mins[:,1])
                    c6  = pyfits.Column(name='max2',   format='E', unit='DN', array=maxs[:,1])
                    c7  = pyfits.Column(name='mean2',  format='E', unit='DN', array=means[:,1])
                    c8  = pyfits.Column(name='median2',format='E', unit='DN', array=medians[:,1])
                    
                    c9  = pyfits.Column(name='min3',   format='E', unit='DN', array=mins[:,2])
                    c10 = pyfits.Column(name='max3',   format='E', unit='DN', array=maxs[:,2])
                    c11 = pyfits.Column(name='mean3',  format='E', unit='DN', array=means[:,2])
                    c12 = pyfits.Column(name='median3',format='E', unit='DN', array=medians[:,2])
                    
                    phdu    = pyfits.PrimaryHDU()
                    head    = phdu.header
                    head['SPEED'] = rdat.speed
                    head['XBIN']  = rdat.xbin
                    head['YBIN']  = rdat.ybin
                    tbhdu   = pyfits.new_table([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12])
                    hdulist = pyfits.HDUList([phdu,tbhdu])
                    hdulist.writeto(stat, clobber=True)
                else:
                    print '*** No frames found in run',run
                    print '*** should not happen'

            except Exception, err:
                print 'Encountered problem on',run
                traceback.print_exc(file=sys.stdout)
