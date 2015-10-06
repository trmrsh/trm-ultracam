"""
A few general utility functions
"""
from __future__ import absolute_import
from __future__ import print_function

import struct
import datetime
try:
    import numpy as np
except ImportError:
    print('Failed to import numpy; some routines will fail')

from trm.ultracam.Constants import MAGIC, MJD0, RUN_DATES, \
    RUN_TELS, BIAS_CHANGES, BIAS_LEVELS
from trm.ultracam.UErrors import UltracamError

def write_string(fobj, strng):
    """
    Writes a string in binary format for my C++ code which
    requires first writing the number of characters and then
    the characters

    fobj         -- file object opened for binary output
    strng        -- string to file object opened for binary output
    """
    nchar = len(strng)
    fobj.write(struct.pack('i' + str(nchar) + 's',nchar, strng))

def read_string(fobj, endian=''):
    """
    Reads a string written in binary format by my C++ code

    fobj   -- file object opened for binary input
    endian -- '>' for big-endian, '' for little-endian.
    """
    nchar  = struct.unpack(endian + 'i', fobj.read(4))[0]
    strng  = struct.unpack(endian + str(nchar) + 's', fobj.read(nchar))[0]
    return strng

def check_ucm(fobj):
    """
    Check a file opened for reading in binary mode to see if it is a ucm.

    Returns endian which is a string to be passed
    to later routines indicating endian-ness.
    """

    # read the format code
    fbytes = fobj.read(4)
    fcode  = struct.unpack('i',fbytes)[0]
    if fcode != MAGIC:
        fcode = struct.unpack('>i',fbytes)[0]
        if fcode != MAGIC:
            fobj.close()
            raise UltracamError('check_ucm: could not recognise first 4 bytes of ' +
                                fname + ' as a ucm file')
        endian = '>'
    else:
        endian = '<'
    return endian

def str2mjd(date):
    """
    Returns an MJD given a YYYY-MM-DD date (can also
    be any other separator, but the YYYY, MM and DD
    must come at the exact same positions)
    """
    year  = date[:4]
    month = date[5:7]
    day   = date[8:11]
    return datetime.date(int(year),int(month),int(day)).toordinal()  - MJD0

# helper routine
def mjd2str(mjd, musec=False):
    """
    Converts an MJD to a string.

    mjd   -- a decimal MJD
    musec -- whether to go to fractions of a second or not
    """
    mjd   += MJD0
    imjd   = int(mjd)
    date   = datetime.date.fromordinal(imjd)
    hour   = 24.*(mjd-imjd)
    ihour  = int(hour)
    mins   = 60.*(hour-ihour)
    imins  = int(mins)
    secs   = 60.*(mins-imins)
    isecs  = int(secs)
    if musec:
        musecs = int(1000000*(secs-isecs))
        tim    = datetime.time(ihour, imins, isecs, musecs)
    else:
        tim    = datetime.time(ihour, imins, isecs)
    dtime  = datetime.datetime.combine(date, tim)
    return dtime.isoformat(' ')

def runID(mjd):
    """
    Identifies a run from an MJD. Returns the run ID and telescope.
    Raises an UltracamError if it cannot match the time.
    """
    for dtup in RUN_DATES:
        run_id, start, stop = dtup
        mstart = str2mjd(start)
        mstop  = str2mjd(stop)
        if mstart < mjd and mjd < mstop+1.5:
            return (run_id,RUN_TELS[run_id])

    raise UltracamError('runID: could not identify time = ' + mjd2str(mjd))

def blevs(mjd, mode):
    """
    Returns tuple of tuples containing typical bias levels for each CCD, and
    for each side of each CCD. This has to be done as a function of time and
    readout mode. Rather than enter the time from the file which can be
    subject to error, it is best to generate the MJD from the night of the run.
    This can be done once and reliably for a given run. Some modes have no
    calibration data and no default values. This routine will return None in
    this case.

    If the return is 'def' than def[nc] is a 2-element tuple containing the left
    and right default bias levels for ccd nc. def[1][0] is thus the default value
    for the left-window of the green CCD.
    """
    BMJDS = np.array(BIAS_CHANGES)
    ind   = np.searchsorted(BMJDS, mjd)
    return BIAS_LEVELS[mode][ind]

if __name__ == '__main__':

    mjd2str(56402)
    runID(56402)

    print('test passed')
