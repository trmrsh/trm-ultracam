#!/usr/bin/env python

"""
trm.ultracam enables Python access to ultracam/spec files. The purpose is to
facilitate development of test and data quality routines. In some aspects it
is reasonably competitive with the pipeline for speed and is easily good
enought for quick scripts.

To access an ultracam raw data file, open it with Rdata which then acts as
a source of the data frames::

  rdat  = Rdata('run045')
  frame = rdat(10)

[trm.ultracam suppressed]. The first line creates an Rdata object with a
connection to the run045 data file. The second reads frame 10 from this
(assuming it exists) returning an MCCD object. This example reads from a local
disk file. You could do the same via the ATC ULTRACAM FileServer as follows::

  rdat  = Rdata('run045',server=True)
  frame = rdat(10)

The 'frame' above is an UCAM object, a sub-class of MCCD, itself a class
designed as a general container of multiple CCD data. For instance::

 red = frame[0]

pulls out the first CCD, the red one, ffrom 'frame'. CCD objects
break down into Window objects, so that::

 for win in red:
   print win 

loops through all Windows of a CCD called 'red'. Window objects
support various operations so that::

 for ccd in frame:
    for win in ccd:
       win -= 100.

would remove 100 from all Windows of all CCDs. MCCD and CCD objects
contain header information in the attribute head::

 print frame.head, frame[0].head

"""

from Constants import *
from Utils import *
from Server import *
from Odict import *
from Window import *
from Time import *
from Uhead import *
from CCD import *
from MCCD import *
from Raw import *
from Log import *
from UErrors import *

