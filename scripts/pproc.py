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

parser = argparse.ArgumentParser(description=usage)

# optional
parser.add_argument('--regex', '-r', help='regular expression for matching specific directories')
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

def rfits(fname, mtype=None, mcomm=None):
    hdul  = pyfits.open(fname)
    head  = hdul[0].header
                
    red   = hdul[1].data
    rlm   = red.field('medl').mean()
    rld   = red.field('p999l').mean()-rlm
    rrm   = red.field('medr').mean()
    rrd   = red.field('p999r').mean()-rrm
    
    green = hdul[2].data
    glm   = green.field('medl').mean()
    gld   = green.field('p999l').mean()-glm
    grm   = green.field('medr').mean()
    grd   = green.field('p999r').mean()-grm
    
    blue  = hdul[3].data
    blm   = blue.field('medl').mean()
    bld   = blue.field('p999l').mean()-blm
    brm   = blue.field('medr').mean()
    brd   = blue.field('p999r').mean()-brm
    
    
    print fname,int(round(rlm)),rld,int(round(rrm)),rrd,int(round(glm)),\
        gld,int(round(grm)),grd,int(round(blm)),bld,int(round(brm)),brd,
    if mtype is not None: print mtype,
    if mcomm is not None: print mcomm,
    print '\n'
    hdul.close()

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
    

