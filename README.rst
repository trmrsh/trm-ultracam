README
======

'trm.ultracam' is a Python module to provide access to raw and ucm ULTRACAM
data files, allow writing to FITS and the like. It is supposed to allow quick
development of scripts, such as may be needed to correct data problems.
Currently it is in "beta testing" and you should not develop anything too 
sophisticated based upon it as I may alter interfaces. This will be the case
until version 1. However, it is quite far along and has some useful stuff,
hence this release.


INSTALLATION
------------

trm.ultracam is "pure python" and tries to be self-contained but 
does require several third-party packages which you may well not have.
Try 'pydoc XXX' for each package XXX to see if you have it. I haven't
been clever about importing these conditionally or anything.
 
In increasing likelihood that you don't have them they are:

  numpy      :  
		numerical Python. If you don't have this already, you
                probably should not be here!

  matplotlib :
                fairly widespread plotting package for Python.


  pyfits     :  
                Python FITS module; I expect to replace this in future
                with astropy

  ppgplot    :  
                (my version.) Get from:
                http://deneb.astro.warwick.ac.uk/phsaap/ppgplot.tar.gz
		This is needed because it is much faster than matplotlib.


A couple of others are used by specific scripts and are not used in the
main API:

  pyds9      :  
                ds9 interface package, used by u2ds9.py, potentially quite 
                a useful script, so you might want this.

  gdata      :  
                google docs Python API. Needed to allow the run checking
                script rchecker.py communicate with the google spreadsheet.
		You won't need this unless you use this script, and if you
                don't know what I am talking about, don't bother.


Once you have these then its the usual::

   python setup.py install

if you are root, or::

   python setup.py install --prefix=<install dir>

if you are not.

If you are primarily interested in the run checker script then try typing
'rchecker.py -h' and then 'rchecker.py -hh' to start with (may need to move 
out of the directory of this README file).

Tom Marsh
