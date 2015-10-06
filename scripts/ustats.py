from __future__ import absolute_import
from __future__ import print_function
from six.moves import range
from six.moves import zip
#!/usr/bin/env python

usage = \
"""
Computes and and stores statistics on all ultracam files it can locate. For
each run it locates, it computes the minimum, median, mean, mode and maximum,
the following percentiles in addition to the median 0.1, 1, 5, 95, 99, 99.9,
the number of occurrences of the mode, the time at mid-exposure and the
exposure time.  It does this for the left- and right-hand sides of each CCD of
each exposure, subject to a user-defined maximum. The results are stored in a
file called runXYZ_stats.fits. The script can only be run from the top-level
directory containing the raw_data and meta_data directories because it will
look through the raw directories but store the results in the corresponding
meta_data directories. It specifically searches for all directories of the
form YYYY-MM-DD, but you can specify a regular expression to narrow the
search. e.g. '2010' will find all runs in directories containing 2010, as well
as having the YYYY-MM-DD format. In all cases it starts from frame 1.
"""

# builtins
import argparse, os, re, sys, traceback

# thirdparty
import pyfits, numpy as np

# mine
from trm import ultracam

parser = argparse.ArgumentParser(description=usage,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# optional
parser.add_argument('--regex', '-r', help='regular expression for matching specific directories')
parser.add_argument('--overwrite', '-o', action='store_true', help='overwrite existing statistics files or not')
parser.add_argument('--max', '-m', type=int, default=10000, help='maximum number of frames to analyse')

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
    print('One or both of',raw,'and',meta,'does not exist or is not a directory.')
    print('Are you running this script from the right directory?')
    exit(1)

for rpath, rnames, fnames in os.walk(raw):

    # search only YYYY-MM-DD directories
    if rdir.search(rpath) and (rext is None or rext.search(rpath)):
        # check for equivalent directory in derived_data
        dpath = meta + rpath[len(raw):]
        if not os.path.exists(dpath):
            print('Directory',rpath,'has no corresponding',dpath,'and will be skipped.')
            break
        
        fnames.sort()
        runs  = [os.path.join(rpath, fname[:-4]) for fname in fnames if rmat.match(fname)]
        stats = [os.path.join(dpath, fname[:-4] + '_stats.fits') for fname in fnames if rmat.match(fname)]

        print('\n\nFound',len(runs),'runs in directory = ',rpath,'\n')

        for run, stat in zip(runs, stats):

            if not os.path.exists(run + '.dat'):
                print(run + '.dat does not exist.')
                continue

            if not args.overwrite and os.path.exists(stat):
                print(stat,'exists and will not be overwritten.')
                continue

            print('Processing',run)

            # initialise outer lists (LH, RH = left-hand, right-hand)
            minls   = [] # LH minimum values, each of which is a list over each CCD
            minrs   = [] # RH minima
            maxls   = [] # LH maxima
            maxrs   = [] # RH maxima
            meanls  = [] # LH means
            meanrs  = [] # RH means
            p01ls   = [] # LH 0.1%-iles
            p01rs   = [] # RH 0.1%-iles
            p1ls    = [] # LH 1%-iles
            p1rs    = [] # RH 1%-iles
            p5ls    = [] # LH 5%-iles
            p5rs    = [] # RH 5%-iles
            medls   = [] # LH medians
            medrs   = [] # RH medians
            p95ls   = [] # LH 95%-iles
            p95rs   = [] # RH 95%-iles
            p99ls   = [] # LH 99%-iles
            p99rs   = [] # RH 99%-iles
            p999ls  = [] # LH 99.9%-iles
            p999rs  = [] # RH 99.9%-iles
            models  = [] # LH modal values
            nmodels = [] # LH numbers of occurrences of modes
            moders  = [] # RH modal values
            nmoders = [] # RH numbers of occurrences of modes
            tims    = [] # MJD times
            exps    = [] # Exposure times
            flags   = [] # timing flags

            try:
                rdat   = ultracam.Rdata(run,flt=False)
                nframe = rdat.ntotal()
                nskip  = nframe // args.max + 1
                numMidnight    = 0
                numFrameErrors = 0

                # we go through every frame, but if nskip > 1, we only read the data
                # 1 in every nskip files. We read all times to ensure that we get good 
                # times in drift and other modes that require them.
                for nf in range(1,nframe):
                    if (nf-1) % nskip == 0:
                        mccd = rdat(nf)
                    else:
                        rdat.time(nf)
                        continue

                    if mccd.head.value('Data.midnight'):
                        numMidnight += 1

                    if mccd.head.value('Data.ferror'):
                        numFrameErrors += 1

                    # same as above, but these contain the values for one frame
                    minl   = []
                    minr   = []
                    maxl   = []
                    maxr   = []
                    meanl  = []
                    meanr  = []
                    p01l   = []
                    p01r   = []
                    p1l    = []
                    p1r    = []
                    p5l    = []
                    p5r    = []
                    medl   = []
                    medr   = []
                    p95l   = []
                    p95r   = []
                    p99l   = []
                    p99r   = []
                    p999l  = []
                    p999r  = []
                    model  = []
                    nmodel = []
                    moder  = []
                    nmoder = []
                    tim    = []
                    exp    = []
                    flag   = []

                    for ccd in mccd:
                        npix = ccd.npix()
                        # treat left and right sides separately
                        larr, rarr = [], []
                        for winl, winr in zip(ccd[::2],ccd[1::2]):
                            larr.append(winl.flatten())
                            rarr.append(winr.flatten())
                        larr = np.concatenate(larr)
                        rarr = np.concatenate(rarr)

                        minl.append(larr.min())
                        maxl.append(larr.max())
                        minr.append(rarr.min())
                        maxr.append(rarr.max())
                        meanl.append(larr.mean())
                        meanr.append(rarr.mean())
                        p01,p1,p5,med,p95,p99,p999 = np.percentile(larr,(0.1,1.,5.,50.,95.,99.,99.9))
                        p01l.append(p01)
                        p1l.append(p1)
                        p5l.append(p5)
                        medl.append(med)
                        p95l.append(p95)
                        p99l.append(p99)                        
                        p999l.append(p999)                        
                        p01,p1,p5,med,p95,p99,p999 = np.percentile(rarr,(0.1,1.,5.,50.,95.,99.,99.9))
                        p01r.append(p01)
                        p1r.append(p1)
                        p5r.append(p5)
                        medr.append(med)
                        p95r.append(p95)
                        p99r.append(p99)                        
                        p999r.append(p999)                        
                        tim.append(ccd.time.mjd)
                        exp.append(ccd.time.expose)
                        flag.append(ccd.time.good)
                        
                        hist = np.bincount(larr)
                        imax = np.argmax(hist)
                        nmax = hist[imax]
                        model.append(imax)
                        nmodel.append(nmax)

                        hist = np.bincount(rarr)
                        imax = np.argmax(hist)
                        nmax = hist[imax]
                        moder.append(imax)
                        nmoder.append(nmax)

                    minls.append(minl)
                    maxls.append(maxl)
                    minrs.append(minr)
                    maxrs.append(maxr)
                    meanls.append(meanl)
                    meanrs.append(meanr)
                    p01ls.append(p01l)
                    p01rs.append(p01r)
                    p1ls.append(p1l)
                    p1rs.append(p1r)
                    p5ls.append(p5l)
                    p5rs.append(p5r)
                    medls.append(medl)
                    medrs.append(medr)
                    p95ls.append(p95l)
                    p95rs.append(p95r)
                    p99ls.append(p99l)
                    p99rs.append(p99r)
                    p999ls.append(p999l)
                    p999rs.append(p999r)
                    tims.append(tim)
                    exps.append(exp)
                    flags.append(flag)
                    models.append(model)
                    nmodels.append(nmodel)
                    moders.append(moder)
                    nmoders.append(nmoder)

                dataProblem = False

            except ultracam.PowerOnOffError as err:
                # silently pass these ones
                pass
            except Exception as err:
                # Can get here after reading some frames if last is
                # partial.
                print('Encountered problem on',run)
                traceback.print_exc(file=sys.stdout)
                dataProblem = True

            if len(minls):
                # only write out a file if some values were found
                tims    = np.array(tims)
                exps    = np.array(exps)
                flags   = np.array(flags)
                minls   = np.array(minls)
                minrs   = np.array(minrs)
                p01ls   = np.array(p01ls)
                p01rs   = np.array(p01rs)
                p1ls    = np.array(p1ls)
                p1rs    = np.array(p1rs)
                p5ls    = np.array(p5ls)
                p5rs    = np.array(p5rs)
                medls   = np.array(medls)
                medrs   = np.array(medrs)
                meanls  = np.array(meanls)
                meanrs  = np.array(meanrs)
                p95ls   = np.array(p95ls)
                p95rs   = np.array(p95rs)
                p99ls   = np.array(p99ls)
                p99rs   = np.array(p99rs)
                p999ls  = np.array(p999ls)
                p999rs  = np.array(p999rs)
                maxls   = np.array(maxls)
                maxrs   = np.array(maxrs)
                models    = np.array(models)
                nmodels   = np.array(nmodels)
                moders    = np.array(moders)
                nmoders   = np.array(nmoders)
                
                # start FITS file construction, primary header
                # followed by tables for each CCD
                phdu             = pyfits.PrimaryHDU()
                head             = phdu.header
                head['NFRAME']   = (nframe,'Total number of frames in file')
                head['NSKIP']    = (nskip,'ustats skip distance(1=all read)')
                head['NMAX']     = (args.max,'maximum number of frames to analyse')
                head['NANAL']    = (len(tims),'actual number of frames analysed')
                head['SPEED']    = (rdat.gainSpeed,'Readout speed hex code')
                head['XBIN']     = (rdat.xbin,'X pixel binning factor')
                head['YBIN']     = (rdat.ybin,'Y pixel binning factor')
                head['NPIX']     = (rdat.npix(),'Total number of binned pixels')
                head['MODE']     = (rdat.mode,'Window readout mode')
                head['MIDNIGHT'] = (numMidnight,'Number of midnight-bug corrections')
                head['FERROR']   = (numFrameErrors,'Number of frame number clashes')
                head['DERROR']   = (dataProblem,'Flags possible problems with raw data file')
                head['']         = '-------'
                head['COMMENT'] = 'File of statistsics with one table per CCD'
                head['COMMENT'] = 'For each frame of a run this records the following:'
                head['COMMENT'] = 'MJD, Expose, Good -- exposure mid-time & length, and timing status flag'
                head['COMMENT'] = 'minl, minr    -- minima of left- and right-hand windows'
                head['COMMENT'] = 'maxl,maxr     -- maxima of left- and right-hand windows'
                head['COMMENT'] = 'meanl,meanr   -- means of left- and right-hand windows'
                head['COMMENT'] = 'medl,medr     -- medians of left- and right-hand windows'
                head['COMMENT'] = 'model,moder   -- modes of left- and right-hand windows'
                head['COMMENT'] = 'nmodel,nmoder -- modes of left- and right-hand windows'
                head['COMMENT'] = 'p1l,p1r       -- 1-percentiles of left- and right-hand windows'
                head['COMMENT'] = 'p5l,p5r       -- 5-percentiles of left- and right-hand windows'
                head['COMMENT'] = 'p95l,p95r     -- 95-percentiles of left- and right-hand windows'
                head['COMMENT'] = 'p99l,p99r     -- 99-percentiles of left- and right-hand windows'
                head['COMMENT'] = 'p999l,p999r   -- 99.9-percentiles of left- and right-hand windows'
                
                nccd = 3 if rdat.instrument == 'ULTRACAM' else 1
                hdul = [phdu,]
                timingError = False
                for nc in range(nccd):
                    c = []
                    snc = str(nc+1)
                    # Create FITS table columns. Group stats on each CCD together
                    c.append(pyfits.Column(name='MJD', format='D', unit='days', array=tims[:,nc]))
                    c.append(pyfits.Column(name='Expose', format='E', unit='secs', array=exps[:,nc]))
                    c.append(pyfits.Column(name='Good', format='L', array=flags[:,nc]))
                    c.append(pyfits.Column(name='minl', format='E', unit='DN', array=minls[:,nc]))
                    c.append(pyfits.Column(name='minr', format='E', unit='DN', array=minrs[:,nc]))
                    c.append(pyfits.Column(name='p01l', format='E', unit='DN', array=p01ls[:,nc]))
                    c.append(pyfits.Column(name='p01r', format='E', unit='DN', array=p01rs[:,nc]))
                    c.append(pyfits.Column(name='p1l', format='E', unit='DN', array=p1ls[:,nc]))
                    c.append(pyfits.Column(name='p1r', format='E', unit='DN', array=p1rs[:,nc]))
                    c.append(pyfits.Column(name='p5l', format='E', unit='DN', array=p5ls[:,nc]))
                    c.append(pyfits.Column(name='p5r', format='E', unit='DN', array=p5rs[:,nc]))
                    c.append(pyfits.Column(name='meanl', format='E', unit='DN', array=meanls[:,nc]))
                    c.append(pyfits.Column(name='meanr', format='E', unit='DN', array=meanrs[:,nc]))
                    c.append(pyfits.Column(name='medl', format='E', unit='DN', array=medls[:,nc]))
                    c.append(pyfits.Column(name='medr', format='E', unit='DN', array=medrs[:,nc]))
                    c.append(pyfits.Column(name='p95l', format='E', unit='DN', array=p95ls[:,nc]))
                    c.append(pyfits.Column(name='p95r', format='E', unit='DN', array=p95rs[:,nc]))
                    c.append(pyfits.Column(name='p99l', format='E', unit='DN', array=p99ls[:,nc]))
                    c.append(pyfits.Column(name='p99r', format='E', unit='DN', array=p99rs[:,nc]))
                    c.append(pyfits.Column(name='p999l', format='E', unit='DN', array=p999ls[:,nc]))
                    c.append(pyfits.Column(name='p999r', format='E', unit='DN', array=p999rs[:,nc]))
                    c.append(pyfits.Column(name='maxl', format='E', unit='DN', array=maxls[:,nc]))
                    c.append(pyfits.Column(name='maxr', format='E', unit='DN', array=maxrs[:,nc]))
                    c.append(pyfits.Column(name='model', format='J', unit='DN', array=models[:,nc]))
                    c.append(pyfits.Column(name='nmodel', format='J', array=nmodels[:,nc]))
                    c.append(pyfits.Column(name='moder', format='J', unit='DN', array=moders[:,nc]))
                    c.append(pyfits.Column(name='nmoder', format='J', array=nmoders[:,nc]))
                    tbhdu = pyfits.new_table(c)
                    hdul.append(tbhdu)

                    # some checks
                    if exps[:,nc].min() < 0 or exps[:,nc].max() > 1000. or \
                            tims[-1,nc]-tims[0,nc] > 0.5 or tims[:,nc].min() < ultracam.FIRST:
                        timingError = True

                head['TERROR'] = (timingError, 'Flags possible problems with the times')
                hdulist = pyfits.HDUList(hdul)
                hdulist.writeto(stat, clobber=True)
            else:
                print('*** No data found in run',run)
