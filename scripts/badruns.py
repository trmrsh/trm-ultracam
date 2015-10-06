#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
import six

usage = \
"""
Runs through all stats files, tries to identify possible problem ones. Optionally writes out an html table.
"""

# builtins
import argparse, os, re, sys, traceback

# thirdparty
import pyfits, numpy as np

# mine
from trm import ultracam

def td(strng='', extra=None):
    if extra is None:
        extra = ''
    else:
        extra = ' class='+extra
    if isinstance(strng, six.string_types):
        if strng == '':
            return '<td'+extra+'>&nbsp;</td>'
        else:
            return '<td'+extra+'>' + strng + '</td>'
    elif strng:
        return '<td'+extra+'>X</td>'
    else:
        return '<td>&nbsp;</td>'


def modeOk(data, npix):
    """
    Sees if mode is OK. npix is total number of CCD pixels
    """
    # left side
    if head['NFRAME'] > 1:
        mode  = data.field('model')[:-1]
        nmode = data.field('nmodel')[:-1]
    else:
        mode  = data.field('model')
        nmode = data.field('nmodel')
    sig    = np.sqrt(2.5**2+np.maximum(mode-4000.,0.))

    # next line is a predicted maximum number of pixels at the mode
    # 1.5 is a safety factor
    nmpred = 1.5*(npix/2.)/(np.sqrt(2.*np.pi)*sig)
    ok     = np.any(nmode < nmpred)
    if ok: return (True,0,0,0,0)
    model  = np.median(mode)
    nmodel = np.median(nmode)

    # right side
    if head['NFRAME'] > 1:
        mode  = data.field('moder')[:-1]
        nmode = data.field('nmoder')[:-1]
    else:
        mode  = data.field('moder')
        nmode = data.field('nmoder')
    sig    = np.sqrt(2.5**2+np.maximum(mode-4000.,0.))

    # next line is a predicted maximum number of pixels at the mode
    # 1.5 is a safety factor
    nmpred = 1.5*(npix/2.)/(np.sqrt(2.*np.pi)*sig)
    ok     = np.any(nmode < nmpred)
    if ok: return (True,0,0,0,0)
    moder  = np.median(mode)
    nmoder = np.median(nmode)

    return (False,model,moder,nmodel/(npix/2.),nmoder/(npix/2.))

html_head = \
"""
<!DOCTYPE HTML>
<html lang="en">
<head>
<title>ULTRACAM problem files</title>
<style type = "text/css">
    table, td, th {
      border: 2px solid black;
    } 
    table {
      border-collapse:collapse;
    }
    th {
      background-color:green;
      color:white;
    }
    td.alert {
      background-color:red;
    }
    dt {
      color:blue;
    }
    </style>
<head>
<body>
<h1>ULTRACAM problem runs</h1>
<p>
This page lists ULTRACAM runs which may suffer from problems identified automatically. See at the end for the meaning of the various
possible problems. Bear in mind that these problems are not necessarily fatal; you just need to be on your guard.

<p>
<table border="2">
<tr>
<th>Date</th><th>Run</th><th>r</th><th>g</th><th>b</th><th>Time</th><th>Midnight</th><th>Frame</th><th>Data</th><th>Mode</th><th>Comment</th>
</tr>
"""

html_foot = \
"""
</table>

<p> The problems are: 
<dl>

<dt>r</dt> 
<dd>red CCD out of range.</dd>

<dt>g</dt>
<dd>green CCD out of range.</dd>

<dt>b</dt>
<dd>blue CCD out of range.</dd>  

<dt>Time</dt>
<dd>some sort of timing problem. Could be no
satellites or somehow incorrectly sequenced.</dd>

<dt>Midnight</dt>
<dd>number of "midnight bugs" detected and
corrected. This bug occurs after midnight UT and shows as a 1 day shift in
times within about 10 seconds of midnight. If it occurs in a run that
straddles midnight, probably all is OK, but runs with timing errors can
sometimes trigger a large number of corrections and may well be in error. </dd>

<dt>Frame</dt>
<dd>In some exposure the frame number was different from 
what was expected.</dd>

<dt>Data</dt>
<dd>some problem with the data, most likely a partial frame</dd>

<dt>Mode</dt>
<dd>too many pixels had the same value, even
in the best frame of the entire run (excluding the final one for runs with > 1 frame). 
This is a good indicator of saturation or worse, and is probably the one flag that 
probably indicates that the data are worthless and hence the red background. If this
one is triggered a comment including modal values and fractions is added.
</dd>

</dl>

</body>
</html>
"""

def extract(data, nmax):
    """
    Returns sample of left and right median values. Chops first
    and last frame if possible.
    """
    lm    = data.field('medl')
    rm    = data.field('medr')
    nf    = max(1,len(lm)-2)
    nskip = max(1,nf // nmax)
    
    # need to copy to allow the fits files to be closed.
    if len(lm) > 2:
        l = lm[1::nskip].copy()
        r = rm[1::nskip].copy()
    else:
        l = lm[::nskip].copy()
        r = rm[::nskip].copy()
    return (l,r)

def badRed(data, nmax):
    """"
    Defines "bad" for the red CCD. Tends to err on the side of
    deciding frames are OK more than not OK.
    """

    # get the data
    l,r = extract(data, nmax)

    # 1580 seems a safe lower limit of all measured frames. Something odd
    # if we are below this.
    if (l < 1580).all():
        return True

    if (l > 60000).all():
        return True

    # Worked out from typical scatter in r-l vs l plot which
    # gradually spreads as the mean counts increase but is also
    # offset from zero 
    diff = r-l-70
    bad  = np.abs(diff) > 30+0.05*np.maximum(0,l - 1700)
    
    if bad.all():
        return True
    return False

def badGreen(data, nmax):
    """"
    Defines "bad" for the green CCD. Tends to err on the side of
    deciding frames are OK more than not OK.
    """
    # get the data
    l,r = extract(data, nmax)

    if (l < 1200).all():
        return True

    if (l > 32000).all():
        return True

    # Worked out from typical scatter in r-l vs l plot which
    # gradually spreads as the mean counts increase but is also
    # offset from zero 
    diff = r-l-10
    bad  = np.abs(diff) > 60+0.05*np.maximum(0,l - 1300)
    
    if bad.all():
        return True
    return False

def badBlue(data, nmax):
    """
    Defines "bad" for the blue CCD. Tends to err on the side of
    deciding frames are OK more than not OK.
    """
    # get the data
    l,r = extract(data, nmax)

    if (l < 1000).all():
        return True

    if (l > 32000).all():
        return True

    # Worked out from typical scatter in r-l vs l plot which
    # gradually spreads as the mean counts increase but is also
    # offset from zero 
    diff = r-l-100.
    bad  = np.abs(diff-(5.5/60.)*l) > 70+0.05*np.maximum(0,l - 1500)
    
    if bad.all():
        return True
    return False

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description=usage)

    # optional
    parser.add_argument('--regex', '-r', \
                            help='regular expression for matching specific directories')
    parser.add_argument('--html', action='store_true', help='generate an html summary page')
    parser.add_argument('-m', dest='nmax', type=int, default=20,
                        help='maximum number of frames per run')

    # OK, done with arguments.
    args = parser.parse_args()

    if args.regex is not None:
        rext = re.compile(args.regex)
    else:
        rext = None

    rdir = re.compile('\d\d\d\d-\d\d-\d\d$')
    rmat = re.compile('^run\d\d\d_stats\.fits$')

    if args.html:
        fht = open('problems.html','w')
        fht.write(html_head)

    # accumulate x, y and colour arrays
    for rpath, rnames, fnames in os.walk('.'):

        # search only YYYY-MM-DD directories
        if rdir.search(rpath) and (rext is None or rext.search(rpath)):
        
            fnames.sort()
            stats = [os.path.join(rpath, fname) for fname in fnames if rmat.match(fname)]

            if len(stats):
#                print 'Found',len(stats),'statistics files in directory = ',rpath,''

                for stat in stats:

                    hdul  = pyfits.open(stat)
                    head  = hdul[0].header
                    
                    modeError = True # flag meaning too many pixels have the modal value
                    probs = ''
                    npix  = head['NPIX']

                    # red
                    hdu   = hdul[1]
                    data  = hdu.data
                    badr = badRed(data, args.nmax)
                    if badr:
                        probs += '; red'
                    rok, rmodel, rmoder, rfmodel, rfmoder = modeOk(data, npix)
                    if rok: modeError = False

                    # green
                    hdu   = hdul[2]
                    data  = hdu.data
                    badg  = badGreen(data, args.nmax)
                    if badg:
                        probs += '; green'
                    gok, gmodel, gmoder, gfmodel, gfmoder = modeOk(data, npix)
                    if gok: modeError = False

                    # blue
                    hdu   = hdul[3]
                    data  = hdu.data
                    badb  = badBlue(data, args.nmax)
                    if badb:
                        probs += '; blue'
                    bok, bmodel, bmoder, bfmodel, bfmoder = modeOk(data, npix)
                    if bok: modeError = False

                    # general checks
                    if 'MIDNIGHT' in head and head['MIDNIGHT']:
                        probs += '; mdnght = ' + str(head['MIDNIGHT'])
                        merror = True
                    else:
                        merror = False

                    if 'FERROR' in head and head['FERROR']:
                        probs += '; frame number'
                        ferror = True
                    else:
                        ferror = False

                    if 'DERROR' in head and head['DERROR']:
                        probs += '; data'
                        derror = True
                    else:
                        derror = False

                    if 'TERROR' in head and head['TERROR']:
                        probs += '; timing'
                        terror = True
                    else:
                        terror = False

                    hdul.close()

                    if modeError:
                        probs += '; histogram'

                    if probs != '':
                        print(stat[:-11]+probs)
                        if args.html:
                            fht.write('<tr>' + td(stat[2:-18]) + td(stat[13:-11]))
                            fht.write(td(badr) + td(badg) + td(badb) + td(terror))
                            if merror: 
                                fht.write(td(str(head['MIDNIGHT'])))
                            else:
                                fht.write(td())
                            fht.write(td(ferror) + td(derror) + td(modeError,'alert'))
                            strng = ''
                            if not rok:
                                strng += 'r: ' + str(rmodel) + ', ' + str(rfmodel) + ', ' + \
                                                 str(rmoder) + ', ' + str(rfmoder)
                            if not gok:
                                if len(strng): strng += '; '
                                strng += 'g: ' + str(gmodel) + ', ' + str(gfmodel) + ', ' + \
                                                 str(gmoder) + ', ' + str(gfmoder)
                            if not bok:
                                if len(strng): strng += '; '
                                strng += 'b: ' + str(bmodel) + ', ' + str(bfmodel) + ', ' + \
                                                 str(bmoder) + ', ' + str(bfmoder)
                            fht.write(td(strng) + '</tr>')

if args.html:
    fht.write(html_foot)
