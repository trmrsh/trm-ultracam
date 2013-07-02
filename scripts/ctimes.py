#!/usr/bin/env python

usage = \
"""
Compares times extracted with pure Python vs pipeline routine gettimes,
stopping when the difference exceeds a preset amount. This is a debug routine.
It runs on every suitable run it can find and stops when it finds an error.
"""

# builtins
import argparse, os, re, sys, traceback, subprocess, datetime

# thirdparty
import numpy as np

# mine
from trm import ultracam

parser = argparse.ArgumentParser(description=usage)

# optional
parser.add_argument('--diff', '-d', type=float, default=10., help='minimum number of milliseconds difference to stop on')
parser.add_argument('--max', '-m', type=int, default=20, help='maximum number of frames to read; 0 for the lot')
parser.add_argument('--number', '-n', type=int, default=0, help='specific run number')
parser.add_argument('--regex', '-r', help='regular expression for matching specific directories')
parser.add_argument('--suppress', '-s', action='store_true', help='suppress all but bad output')

# OK, done with arguments.
args = parser.parse_args()

if args.regex is not None:
    rext = re.compile(args.regex)
else:
    rext = None

rdir = re.compile('\d\d\d\d-\d\d-\d\d$')
rmat = re.compile('^run\d\d\d\.xml$')

# the ULTRACAM times command path
times = os.path.join(os.environ['TRM_SOFTWARE'], 'bin', 'ultracam', 'times')

# Special cases, not worth fixing, so just skip. The
# first two start with no satellites which causes difficulties.
SPECIALS = ['./2003-11-12/run001','./2004-05-18/run008']

for rpath, rnames, fnames in os.walk('.'):

    # search only YYYY-MM-DD directories
    if rdir.search(rpath) and (rext is None or rext.search(rpath)):
        
        fnames.sort()
        runs  = [os.path.join(rpath, fname[:-4]) for fname in fnames if rmat.match(fname)]

        print '\n\nFound',len(runs),'runs in directory = ',rpath,'\n'

        for run in runs:

            if args.number:
                if int(run[-3:]) != args.number:
                    continue
            
            print 'Processing',run
            if not os.path.exists(run + '.dat'):
                print '....... dat file does not exist.'
                continue

            if run in SPECIALS:
                print 'Special case; will skip it'
                continue

            try:
                # first run the Python version
                gps, mjd, expose, res = [], [], [], []
                rtime = ultracam.Rtime(run)

                for time in rtime:
                    tim,blueTime,badBlue,info = time
                    gps.append(info['gps'])
                    mjd.append(tim.mjd)
                    expose.append(tim.expose)
                    res.append(tim.reason)
                    if args.max and len(mjd) == args.max:
                        break
                
                # build the pipeline equivalent command, run it, read the results.
                comm = (times, 'source=l', 'file=' + run, 'first=1', 'last='+str(args.max), 
                        'clock=yes', 'twait=1', 'tmax=0', 'nodefs')
                sout, serr = subprocess.Popen(comm, stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE).communicate()
                output = sout.split('\n')
                nline = 0
                for line in output:
                    if not line.startswith('#') and nline < len(mjd) and not line.startswith('U'):
                        nf,gp,mj,ok,expos,date = line.split('|')
                        gp, mj, expos = float(gp), float(mj), float(expos)
                        dgp = 1000*ultracam.DSEC*(gps[nline]-gp)
                        dmj = 1000*ultracam.DSEC*(mjd[nline]-mj)
                        dex = 1000*(expose[nline]-expos)
                        if abs(expos) < 10000. and (abs(dgp) > args.diff or abs(dmj) > args.diff or abs(dex) > args.diff):
#                            raise Exception('%-4d GPS=%12.6f %8.2f, MJD=%12.6f %8.2f, Exp=%6.3f %8.2f, GPS,MJD=%s, %s' % \
#                                                (nline,gp,dgp,mj,dmj,expos,dex,dstr(gp),dstr(mj)))
                            print 'format,which,mode,instrument,defTstamp,vclock_frame =',\
                                frmt,whichRun,rtime.mode,rtime.instrument,defTstamp,vclock_frame,res[nline]
                            print '%-4d GPS=%12.6f %8.2f, MJD=%12.6f %8.2f, Exp=%6.3f %8.2f, GPS,MJD=%s, %s' % \
                                (nline+1,gps[nline],dgp,mjd[nline],dmj,expose[nline],dex,
                                 ultracam.mjd2str(gp,True),ultracam.mjd2str(mj,True))
                            exit(1)
                        nline += 1
                        if args.max and nline == args.max:
                            break

            except ultracam.PowerOnOffError, err:
                print '....... power on/off; ignoring'

            except Exception, err:
                print 'Encountered problem on',run
                traceback.print_exc(file=sys.stdout)
