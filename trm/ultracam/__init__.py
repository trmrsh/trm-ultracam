#!/usr/bin/env python

"""
**trm.ultracam** is a module to access ultracam/spec files through Python. 
Its purpose is to facilitate development of test and data quality routines. 
In some aspects it is reasonably competitive with the C++ pipeline for speed 
and is easily good enough for quick scripts.

Good starting points are the :class:`trm.ultracam.Rdata` and :class:`trm.ultracam.MCCD` 
classes.
"""
from __future__ import absolute_import

from .Constants import *
from .Utils import *
from .Server import *
from .Odict import *
from .Window import *
from .Time import *
from .Uhead import *
from .CCD import *
from .MCCD import *
from .Raw import *
from .Log import *
from .UErrors import *

__all__ = ['Rhead', 'Rdata', 'Rtime', 'Rwin', 'CCD', 'MCCD', 'utimer']


