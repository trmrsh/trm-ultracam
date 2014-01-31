#!/usr/bin/env python

"""
uspfix is designed to fix a timing problem in ULTRASPEC files that was spotted
in January 2014. When the problem occurs (it is intermittent in nature) a null
timestamp appears in one frame of a run (we have no example in which more than
one frame is affected). This is extra in the sense that all the proper
timestamps are shifted back by one frame. We suspect that a null timestamp
enters the FIFO buffer somehow and will be investigating fixing it in the
future. In the meantime, this script fixes any data affected by it.

If the problem occurs, it tends to happens early on, typically, but not
always, on the second frame of a run. This means that all subsequent frames
effectively have timestamps that are too early by one frame. For precision
timing work, this is a serious error. uspfix looks for such frames, then
copies the timing bytes backwards by one frame, beginning with the first frame
after the rogue frame. The bytes are copies so the final frame keeps the same
set of bytes, and thus will not have a reliable timestamp although it may be
possible to recover a proper mid-exposure time because sometimes only
preceding timestamps are needed.

uspfix runs on all the runs it finds (defined by the presence of files of the
form run###.xml) in the present working directory, and works down the
directory hierarchy from there. Thus if you run it in a directory containing
the night-by-night files 2014-02-01, 2014-02-02 etc, it will fix *all* bad
runs in the top-level directory, these directories, and in any sub-directories
they may contain.  It skips any runs for which it cannot read the .dat file
or those with <= 32 bytes per frame (which picks up powerons).  It makes
copies of any files that it changes to run###.dat.old for safety.

Simply run the script in the directory of interest, with no arguments, but
remember that it will run in *all* sub-directories as well. The data must be
write enabled because it is modified by the script. e.g. in unix: "chmod +w
run*.dat" before running the script.

If there are any problems, please contact Tom Marsh at Warwick.
"""

import os, re, shutil

# This to spot rogue timestamps
BLANK = 20*'\x00'

# To match the xml files
rmat = re.compile('^run\d\d\d\.xml$')

# Work through the directory tree starting from the pwd
for rpath, rnames, fnames in os.walk('.'):
    runs  = [fname[:-4] for fname in fnames if rmat.match(fname)]
    runs.sort()

    # Go through all the runs found
    for run in runs:
        root = os.path.join(rpath,run)

        # Check for .dat file
        if not os.path.exists(root + '.dat'):
            print root,'skipped; no ,dat file'
            break

        # Determine the filesize
        fsize = 0
        with open(root + '.xml') as fxml:
            for line in fxml:
                ptr = line.find('framesize')
                if ptr >= 0:
                    fq = ptr + line[ptr:].find('"') + 1
                    sq = fq + line[fq:].find('"')
                    fsize = int(line[fq:sq])

        # PowerOns are 32 bytes
        if fsize <= 32:
            print root,'skipped; probably a poweron'
            continue

        # Search for the timing bug
        badFrame = 0
        with open(root + '.dat','rb') as fin:
            ptr = 12
            nf  = 0
            while True:
                nf += 1
                fin.seek(ptr)
                timing = fin.read(20)
                if len(timing) != 20:
                    break
                if timing == BLANK:
                    badFrame = nf
                    break
                ptr += fsize

        if badFrame:
            # First of all copy the file for safety.
            shutil.copyfile(root + '.dat', root + '.dat.old')

            # Now modify the original. We open it to be able to both read from
            # and write to it in binary mode:
            with open(root + '.dat','r+b') as fio:
                # start from the corrupted frame
                n = badFrame
                while True:
                    # Seek the timing bytes we want (which come immediately
                    # after the frame we are about to modify)
                    fio.seek(12+fsize*n)

                    # Read them, if possible
                    timing = fio.read(20)
                    if len(timing) != 20:
                        break

                    # Write them into preceding location
                    fio.seek(12+fsize*(n-1))
                    fio.write(timing)

                    # Advance to next frame
                    n += 1

            # Report progress
            print root,'corrected; corrupted file copied to',root + '.dat.old'

