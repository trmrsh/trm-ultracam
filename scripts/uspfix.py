from __future__ import absolute_import
from __future__ import print_function
#!/usr/bin/env python

usage = \
"""
If this is your first encounter with this script, please read ALL of the
following. In brief the script corrects a timing bug in ultraspec data and can
be run either as:

 uspfix.py      (corrects .dat files in directory it is run in)

or:

 uspfix.py -r   (runs recursively over directory hierarchy starting from
                 current directory)

You should make the .dat file writeable. e.g. this can be done recursively
as follows:

 find . -name "run*.dat" | xargs chmod +w


In more detail:

uspfix.py is designed to fix a timing problem in ULTRASPEC files that was
spotted in January 2014. When the problem occurs (it is intermittent in
nature) a null timestamp appears in one frame of a run (we have no example in
which more than one frame is affected). This is extra in the sense that all
the proper timestamps are shifted by one frame later in sequence. We suspect
that a null timestamp enters the FIFO buffer somehow and will be investigating
fixing it in the future. In the meantime, this script fixes any data affected
by it. The symptom of the bug is a line like this, as reported e.g. by the
pipeline command "rtplot":

"Ultracam::read_header WARNING: time unreliable: GPS clock not yet synced
since power up"

This is very likely spurious, especially if subsequent times appear to be OK.
If the problem occurs, it tends to happens early on, typically, but not
always, on the second frame of a run. This means that all subsequent frames
effectively have timestamps that are too early by one frame. For precision
timing work, this is a serious error as one can often time events to much
higher precision than the sampling time. uspfix.py works by searching for
frames in which all the timing bytes are set to 0, then copies the timing
bytes backwards by one frame, beginning with the first frame after the rogue
frame. The bytes are copied, so the final frame keeps its timing bytes intact,
and thus will not have a reliable timestamp, although it may be possible to
recover a proper mid-exposure time because sometimes only preceding timestamps
are needed.

uspfix.py runs on all the runs it finds (defined by the presence of files of
the form run###.xml) in the present working directory, and works down the
directory hierarchy from there. Thus if you run it in a directory containing
the night-by-night directories 2014-02-01, 2014-02-02 etc, it will fix *all*
bad runs in the top-level directory, these directories, and in any
sub-directories they may contain.  It skips any runs for which it cannot read
the .dat file, those with <= 32 bytes per frame (which picks up powerons),
those in which the bad frame is the last, since they are uncorrectable, and
finally any in which there is more than one null timestamp. The latter
exception is because we have never encountered such a run, and I don't
therefore want to try correcting them. If you find one of these, please e-mail
me (Tom Marsh, Warwick). uspfix.py makes copies of any files that it changes
to run###.dat.old for safety.

Simply run the script in the directory of interest, with no arguments, but
remember that it will run in *all* sub-directories as well. The data must be
write enabled because they are modified by the script. e.g. in unix:

 find . -name "run*.dat" | xargs chmod +w

to make all the run.dat files in the directory tree writeable before running
the script.

The script reports all the runs that it corrects and those that it cannot read
(e.g. powerons). If you run it twice, it should make no corrections the second
time. The script modifies the data in place bringing with it some possibility
of corrupting data irretrievably. You are advised not the ctrl-C it, although
there is a trap in place in the critical section to prevent this causing
problems. This section only takes a fraction of the total time of the routine,
so it will probably be OK.

Please contact Tom Marsh at Warwick if you encounter problems with this script
and if you encounter any run with more than 1 bad frame as reported by the
script.
"""

import sys, os, re, shutil, signal

if len(sys.argv) > 1:
    if sys.argv[1] == '-r':
        recursive = True
    else:
        print(usage)
        exit(1)
else:
    recursive = False

# This to spot rogue timestamps
BLANK = 20*'\x00'

# To match the xml files
rmat = re.compile('^run\d\d\d\.xml$')

# Compile a list of runs
if recursive:
    runs = []
    for rpath, rnames, fnames in os.walk('.'):
        runs.extend( [os.path.join(rpath,fname[:-4]) for fname in fnames \
                 if rmat.match(fname)] )
else:
    runs = [fname[:-4] for fname in os.listdir('.') \
            if rmat.match(fname)]
runs.sort()

# check that .dat files have write permission
for run in runs:

    # Check for .dat file
    path = run + '.dat'
    if os.path.exists(path) and not os.access(path, os.W_OK):
        print('You do not have write access to',path)
        print('Please fix write permissions on the .dat files and then re-run')
        print('the script. The following command does this recursively:\n')
        print('find . -name "run*.dat" | xargs chmod +w\n')
        exit(1)

# Now do the work
for run in runs:

    # Check for .dat file
    if not os.path.exists(run + '.dat'):
        print(run,'skipped; no .dat file')
        continue

    # Determine the filesize
    fsize = 0
    with open(run + '.xml') as fxml:
        for line in fxml:
            ptr = line.find('framesize')
            if ptr >= 0:
                fq = ptr + line[ptr:].find('"') + 1
                sq = fq + line[fq:].find('"')
                fsize = int(line[fq:sq])

    # PowerOns are 32 bytes
    if fsize <= 32:
        print(run,'skipped; probably a poweron')
        continue

    # Search for bad timestamps, storing the first
    badFrame = 0
    nbad = 0
    with open(run + '.dat','rb') as fin:
        ptr = 12
        nf  = 0
        while True:
            fin.seek(ptr)
            timing = fin.read(20)
            if len(timing) != 20:
                break
            nf += 1
            if timing == BLANK:
                if badFrame == 0: badFrame = nf
                nbad += 1
            ptr += fsize

    # Skip if more than one badframe is found because
    # I don't know how to deal with these properly.
    if nbad > 1:
        print(run,'skipped as it had',nbad,'(>1) bad frames')
        print('The times of this run cannot be trusted!!')
        print('Please e-mail Tom Marsh at Warwick if you encounter this')
        print('in a run that has *some* good times.')
        continue

    if badFrame == nf:
        print(run,'skipped as the bad frame was the final one')
        continue

    if badFrame:
        # First of all copy the file for safety.
        shutil.copyfile(run + '.dat', run + '.dat.old')

        # Block ctrl-C interrupts ...
        signal.signal(signal.SIGINT, signal.SIG_IGN)

        # Now modify the original. We open it to be able to both read from
        # and write to it in binary mode:
        with open(run + '.dat','r+b') as fio:
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
        print(run,'corrected; corrupted file copied to',run + '.dat.old')

        # Unblock ctrl-C interrupts ...
        signal.signal(signal.SIGINT, signal.SIG_DFL)
