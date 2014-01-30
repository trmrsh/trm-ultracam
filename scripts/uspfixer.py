#!/usr/bin/env python

usage = \
""" uspfixer is designed to fix a timing problem that was spotted in January
2014 in which an extra rogue timestamp appears in ultraspec runs, shifting all
subsequent ones back by a frame. If not fixed then the time assigned to the
subsequent frames are too early by one frame, which, depending upon the
context, can be highly significant. uspfixer looks for such frames, then
copies the timing bytes backwards by one frame, beginning with the first frame
after the rogue frame. The bytes are copies so the final frame keeps the same
set of byte, and thus will not have a reliable time.  uspfixer runs on all the
runs it finds (defined by the presence of files of the form run###.xml) in the
present working directory. It skips any for which it cannot read the .dat
file. It makes copies of the old data files which have the extension .dat.old.

Note that the final frame of the corrected files will not have the right time.

Simply run the script in the directory of interest, with no arguments.
"""

import argparse, os, re, sys
from trm import ultracam

parser = argparse.ArgumentParser(description=usage)
args = parser.parse_args()

rmat = re.compile('^run\d\d\d\.xml$')

for rpath, rnames, fnames in os.walk('.'):
    runs  = [fname[:-4] for fname in fnames if rmat.match(fname)]
    runs.sort()
    for run in runs:
        try:
            tdat   = ultracam.Rtime(run)
            fsize  = tdat.framesize

            # Find first badframe (if any)
            badFrame = 0
            for nf, time in enumerate(tdat):
                if nf:
                    if (not time[0].good) and \
                       (time[0].reason == 'timestamp too early'):
                        badFrame = nf + 1
                        break
            tdat.close_file()

            if badFrame:
                print run,'has an invalid timestamp in frame',badFrame
                print 'Will now shift times of later frames one forward.'
                with open(run + '.dat.new','wb') as fout:
                    with open(run + '.dat','rb') as fin:

                        # Directly copy the uncorrupted frames leading
                        # up to the failed frame
                        for nf in xrange(badFrame-1):
                            frame = fin.read(fsize)
                            fout.write(frame)

                        # Read the frame with the corrupted timestamp
                        frame1 = bytearray(fin.read(fsize))

                        # Read the succeeding frames and carry out the
                        # timing data transfer
                        while True:
                            frame2 = bytearray(fin.read(fsize))
                            if len(frame2) != fsize: break

                            # Transfer the timing bytes back by a frame and
                            # write the corrected frame out
                            frame1[12:32] = frame2[12:32]
                            fout.write(frame1)
                            frame1 = frame2

                    # At this point the file read has stopped, leaving the
                    # final frame, referenced by frame1, unwritten, so we
                    # write it out and stop. NB we can't fix the time of
                    # this one and just leave the time it has (which has also
                    # been copied to the previous frame) untouched.
                    fout.write(frame1)

#                os.rename(run + '.dat', run + '.dat.old')
                os.rename(run + '.dat.new', run + '.dat')
#                print run,'corrected; old .dat file copied to',run + '.dat.old'
                print run,'corrected'
            else:
                print run,'is OK and has not been changed.'
        except Exception, err:
            print run,'could not be read and has not been changed'
