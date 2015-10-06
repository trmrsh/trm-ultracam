#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
from six.moves import urllib


usage = \
"""
Produces alerts if it detects problems with the current ULTRACAM run 
during observing. It works by polling the FileServer every so often for
the latest frame. It produces alerts if too long a period of time has 
elapsed since the frame changed or if there are problems with the data.
"""

# builtins
import argparse, time

# mine
from trm import ultracam

# argument parsing
parser = argparse.ArgumentParser(description=usage,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# optional
parser.add_argument('-t', dest='tmax', type=int, default=10, help='maximum time before warning of no new frame')
parser.add_argument('-w', dest='wait', type=int, default=10, help='number of seconds wait between updates')

# OK, done with arguments.
args = parser.parse_args()

# these indicate the most recently read run, frame, and the last time
# anything changed
lastRun    = None
lastNframe = None
lastNew    = time.time()

print('\nAlerter script started. Will poll once every',args.wait,\
    'seconds and will warn of no change after',args.tmax,'minutes.')
print("""
It reads the final frame of the most recent run and makes some fairly crude
checks for problems (too many pixels of the same value, too high, too low, too
large a difference between right and left CCDs). Please let trm know if you
think it is being either too sensitive or not sensitive enough, and how it is
going wrong. For options to the script, invoke it with -h.

""")
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
            newrun     = currentRun != lastRun
            lastRun    = currentRun

            # Check whether the number of frames has changed
            nframe      = ultracam.get_nframe_from_server(currentRun)
            newframe    = newrun or nframe != lastNframe
            lastNframe  = nframe

            # Update last time anything was new
            lastNew = tstamp if newframe else lastNew

            # Issue alert if nothing seems to be happening
            if tstamp - lastNew >= 60*args.tmax:
                print(uttime + ': >>>>>>> WARNING: Nothing has changed for ' + \
                    str(int((tstamp-lastNew)/ 60.)) + ' minutes! <<<<<<<<')

            # update the Rdata object if the run has changed
            if newrun:
                rdat = ultracam.Rdata(currentRun,flt=False,server=True)

            # Read the last frame if it has changed and there is one
            if newframe and nframe:
                mccd = rdat(nframe)
            
                # check it
                r,g,b = mccd.checkData()
                if r[0] or g[0] or b[0]:
                    print(uttime,'>>>>>>> WARNING:', end=' ')
                    if r[0]: print(' red: ' + r[1] + '.', end=' ')
                    if g[0]: print(' green: ' + g[1] + '.', end=' ')
                    if b[0]: print(' blue: ' + b[1] + '.', end=' ')
                    print(' <<<<<<<<')
                else:
                    print(uttime + ': all nominal. Currently on',currentRun,'which has',nframe,'frames.')

            elif tstamp - lastNew < 60*args.tmax:
                print(uttime + ': all nominal. Currently on',currentRun,'which has',nframe,'frames.')

        else:
            if tstamp - lastNew > 60*args.tmax:
                print(uttime + ': >>>>>>> WARNING: Nothing has changed for ' + \
                    str(int((tstamp-lastNew)/ 60.)) + ' minutes! <<<<<<<<')

    except urllib.error.URLError as err:
        print(uttime + ': ' + str(err) + '; have you started the ATC FileServer?')
    except ultracam.UltracamError as err:
        print(uttime + ': ' + str(err))

    # now wait before polling the server again
    time.sleep(args.wait)
