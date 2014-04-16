#!/usr/bin/env python

usage = \
"""
Prints warnings to the terminal if it detects problems with the timing of the
current ULTRACAM run during observing. It polls the FileServer for the times
of every frame and tries to spot any that come out of sequence. Initially just
try with default arguments:

talert.py

"""

# builtins
import argparse, time, urllib2
import numpy as np

# mine
from trm import ultracam

# argument parsing
parser = argparse.ArgumentParser(description=usage,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# optional
parser.add_argument('-t', dest='tmax', type=int, default=30,
                    help='maximum time before warning of no new frame (secs)')
parser.add_argument('-w', dest='wait', type=int, default=10,
                    help='number of seconds wait between updates')
parser.add_argument('-f', dest='fmax', type=float, default=0.1,
                    help='Maximum fractional deviation in interval')
parser.add_argument('-n', dest='nmin', type=int, default=10,
                    help='Minimum number of frames to accumulate before reporting problems')
parser.add_argument('-b', dest='bmax', type=int, default=0,
                    help='Maximum number of times to read before pausing, 0 to ignore')

# OK, done with arguments.
args = parser.parse_args()

# these indicate the most recently read run, frame, and the last time
# anything changed
lastRun    = None
lastNframe = None
lastNew    = time.time()

print \
"""
Time alert script started. Will poll once every {0:d}
seconds and will warn of no change after {1:d} seconds.
It reads every frame of the most recent run looking for
times that do not occur at close to the median interval
after the previous one.
""".format(args.wait,args.tmax)

while True:

    # get a couple of times just once
    uttime  = time.asctime(time.gmtime())
    tstamp  = time.time()

    # get list of runs from the FileServer
    try:
        runs = ultracam.get_runs_from_server()

        if len(runs):

            # Check whether the final run has changed
            currentRun = runs[-1]

            # skip poweron/offs
            rhead = ultracam.Rhead(currentRun, server=True)
            if rhead.isPonoff(): continue

            newrun     = currentRun != lastRun
            lastRun    = currentRun

            # Check whether the number of frames has changed
            nframe      = ultracam.get_nframe_from_server(currentRun)
            newframe    = newrun or nframe != lastNframe
            lastNframe  = nframe

            # Update last time anything was new
            lastNew = tstamp if newframe else lastNew

            # Issue alert if nothing seems to be happening
            if tstamp - lastNew >= args.tmax:
                print uttime + ': >>>>>>> WARNING: Nothing has changed for ' + \
                    str(int(tstamp-lastNew)) + ' seconds! <<<<<<<<'

            # update the Rdata object if the run has changed
            if newrun:
                tdat  = ultracam.Rtime(currentRun,server=True)
                frame = 0
                # use this to accumulate time intervals
                diffs = []

            if args.bmax:
                nmax = min(nframe, frame+args.bmax)
            else:
                nmax = nframe

           # Read the times of the frame up to the current frame
            if nmax > frame:
                for nf in xrange(frame+1,nmax+1):
                    mjd = tdat(nf)[1]['gps']
                    if nf > 1:
                        diffs.append(86400.*(mjd-mjdold))
                    mjdold = mjd

                if len(diffs) > args.nmin:
                    # only start checking when we have a few in the bag
                    adiffs = np.array(diffs)
                    mdiff  = np.median(adiffs)

                    devs   = np.abs(adiffs[frame:] - mdiff)
                    bad    = devs > args.fmax*mdiff
                    if len(devs[bad]):
                        for nb, ng in enumerate(bad):
                            if ng:
                                print 'WARNING: run ' + str(currentRun) + ', frame',\
                                    frame+nb+2,'occurred',adiffs[frame+nb],\
                                    'secs after previous cf median =',mdiff
                    else:
                        print 'Run ' + str(currentRun) + ', frames',frame+1,'to',\
                            nmax,'have OK times.'

                frame = nmax

        else:
            if tstamp - lastNew > args.tmax:
                print uttime + ': >>>>>>> WARNING: Nothing has changed for ' + \
                    str(int(tstamp-lastNew)) + ' seconds! <<<<<<<<'

    except urllib2.URLError, err:
        print uttime + ': ' + str(err) + '; have you started the ATC FileServer?'
    except ultracam.UltracamError, err:
        print uttime + ': ' + str(err)

    # now wait before polling the server again
    time.sleep(args.wait)
