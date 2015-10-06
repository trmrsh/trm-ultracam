from __future__ import absolute_import
from __future__ import print_function
#!/usr/bin/env python
#
# Reads 3D FastCam data, dumps to a series of ucm files.
#
# Author: T.R.Marsh

usage = \
"""
fcam3d2ucm converts a 3D FITS file from FastCam into ucm files. These will
be named according to the root of the FastCam file name with a sequence number.
e.g. abc.fits would produce files of the form 'abc_0001.ucm', 'abc_0002.ucm', 
etc. The number of digits is 4 by default but can be varied.
"""

import argparse
import astropy.io.fits as fits
from astropy.time import Time
import trm.ultracam as ucam

parser = argparse.ArgumentParser(description=usage)

# positional
parser.add_argument('fcam',  help='name of FastCam file')

# optional
parser.add_argument('-d', dest='ndigit', type=int,
                    default=4, help='number of digits in file names to store sequence number')

# OK, done with arguments.
args = parser.parse_args()

fname  = args.fcam
ndigit = args.ndigit

if not fname.endswith('.fits'):
    print('ERROR: File must end with .fits')
    exit(1)

# Open file
hdul = fits.open(fname)
head = hdul[0].header

naxis = head['NAXIS']
if naxis != 3:
    print('ERROR: Was expecting a 3D FITS file')
    exit(1)

source = head['SOURCE']
if source != 'FastCam':
    print('ERROR: Source =',source,'but was expecting "FastCam"')
    exit(1)

data = hdul[0].data


# Timing info
exptime = head['EXPTIME'] # exposure time, milliseconds

date = head['DATE'] # date at start [? check]
start_time = Time(date, scale='utc')
mjd = int(start_time.mjd)

ut_start = head['UT']
h,m,s = ut_start.split(':')
hours = float(h) + float(m)/60. + float(s)/3600.

exptime /= 1000.  # convert to seconds

# Get CCD & window parameters
nxmax, nymax = head['CCDSIZE'].split('X')
nxmax, nymax = int(nxmax), int(nymax)
llx, lly = eval(head['IMGORIGIN'])
nx, ny = head['NAXIS1'], head['NAXIS2']
xbin, ybin = head['BINNING'].split('X')
xbin, ybin = int(xbin), int(ybin)

# Few other params
target = head['OBJNAME']
filter = head['FILTER']

# Loop through data

# Offset times to centre of exposures [? check]
hours += exptime/3600./2.

for n, image in enumerate(data):
    win = ucam.Window(image, llx, lly, xbin, ybin)
    utime = ucam.Time(mjd+hours/24., exptime, True, '')
    ccd = ucam.CCD([win,], utime, nxmax, nymax, True, None)

    # Create header
    head = ucam.Uhead()

    head.add_entry('Object',target,ucam.ITYPE_STRING,'Object name')
    head.add_entry('Filter',filter,ucam.ITYPE_STRING,'Filter')
    head.add_entry('Date',date,ucam.ITYPE_STRING,'Date at start of run')
    head.add_entry('UT',ut_start,ucam.ITYPE_STRING,'UT at start of run')
    head.add_entry('UT_date',(mjd,hours),ucam.ITYPE_TIME,'MJD at centre of exposure')
    head.add_entry('Exposure',exptime,ucam.ITYPE_FLOAT,'Exposure time, seconds')
    head.add_entry('Frame','Frame info')
    head.add_entry('Frame.reliable',True,ucam.ITYPE_BOOL,'Time reliable')
    head.add_entry('Instrument','Instrument setup')
    head.add_entry('Instrument.Name','FastCam',ucam.ITYPE_STRING,
                   'Instrument name')
    head.add_entry('Instrument.nblue',1,ucam.ITYPE_INT,'Blue cycle number')

    # Create MCCD 
    mccd = ucam.MCCD([ccd,], head)

    # Write to ucm file
    oname = ('{0:s}_{1:0' + str(ndigit) + 'd}.ucm').format(fname[:-5], n+1)
    mccd.wucm(oname)

    if (n+1) % 10 == 0:
        print('Dumped',oname,'to disk.')

    # Shift time to next exposure
    hours += exptime/3600.
    if hours >= 24.:
        hours -= 24.
        mjd += 1


