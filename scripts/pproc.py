#!/usr/bin/env python

usage = \
"""
Looks for files of the form YYYY-MM-DD/runXXX_stats.fits and prints
out one liner summary info. The aim is to use this to auto ID data
types to at least some extent.
"""

# builtins
import argparse, os, re, sys, traceback

# thirdparty
import pyfits, numpy as np

# mine
from trm import ultracam

def stats(hdu, dtype):
    data  = hdu.data
    lm    = data.field('medl').mean()
    rm    = data.field('medr').mean()
    
    # detect an asymmetry in 5 and 95%-iles relative to the median
    # indicative of the presence of real counts.
    asymm = data.field('p5l').mean()+data.field('p95l').mean()+\
        data.field('p5r').mean()+data.field('p95r').mean()-\
        2*lm-2*rm
    
    # next should be roughly 8 sigma
    spread = data.field('p95l').mean()-data.field('p5l').mean()+\
        data.field('p95r').mean()-data.field('p5r').mean()
    mean = (lm+rm)/2.
    diff = rm-lm
    print '[%5d %7.1f %7.1f %7.1f]' % (mean, diff, asymm, spread),

    if dtype['JUNK'] and spread > 1: dtype['JUNK'] = False
    if dtype['BIAS'] and (spread < 2 or spread > 100 or asymm > 10): dtype['BIAS'] = False

def rfits(fname, mtype=None, mcomm=None):
    """"
    Runs over a whole file, tries to collate results.
    """
    hdul  = pyfits.open(fname)
    head  = hdul[0].header

    dtype = {'JUNK' : True, 'BIAS' : True, 'SCIENCE' : True, 'TECHNICAL' : True, 
             'SKY' : True, 'ACQUISITION' : True}

    print fname,
    stats(hdul[1], dtype)
    stats(hdul[2], dtype)
    stats(hdul[3], dtype)

    if mtype is not None: print mtype,

    autoid = ''
    for key, value in dtype.iteritems():
        if value: autoid += '[' + key + ']' 
    print autoid,
    if mcomm is not None: print mcomm,
#    print
    hdul.close()


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description=usage)

    # optional
    parser.add_argument('--regex', '-r', \
                            help='regular expression for matching specific directories')
    parser.add_argument('--file', '-f', type=argparse.FileType('r'), \
                            help='file of manual IDs')

    # OK, done with arguments.
    args = parser.parse_args()

    if args.regex is not None:
        rext = re.compile(args.regex)
    else:
        rext = None

    rdir = re.compile('\d\d\d\d-\d\d-\d\d$')
    rmat = re.compile('^run\d\d\d_stats\.fits$')

    if args.file is None:
        for rpath, rnames, fnames in os.walk('.'):

            # search only YYYY-MM-DD directories
            if rdir.search(rpath) and (rext is None or rext.search(rpath)):
        
                fnames.sort()
                stats = [os.path.join(rpath, fname) for fname in fnames if rmat.match(fname)]

                if len(stats):
                    print '\n\nFound',len(stats),'statistics files in directory = ',rpath,'\n'

                    for stat in stats:
                        rfits(stat)

    else:
        for line in args.file:
            if not line.startswith('#'):
                arr   = line.split()
                fname = arr[0]+'_stats.fits'
                mtype = arr[1]
                mcomm = line[line.find(mtype)+len(mtype):]
                rfits(fname,mtype,mcomm)
    
