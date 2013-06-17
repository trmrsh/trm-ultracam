#!/usr/bin/env python

"""
enables Python access to ultracam/spec files. The purpose is to facilitate 
development of test and data quality routines. In some aspects it is reasonably
competitive with the pipeline for speed.

To access an ultracam raw data file, open it with Rdata which then acts as
a source of the data frames using get. e.g.

 rdat  = Udata('run045')
 frame = rdat.get(10)

The first line creates an object with a connection to the run045 data. 
The second reads frame 10 from this (assuming it exists) returning an
MCCD object. For more examples of usage, see the script demo.py
"""

# standard imports
import xml.dom.minidom
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm     as cm
import struct
import warnings
import datetime

# less widely known extras
import ppgplot as pg

# kick off with some semi-hidden helper routines
def _write_string(fobj, strng):
    """
    Writes a string in binary format for my C++ code which
    requires first writing the number of characters and then 
    the characters

    fobj         -- file object opened for binary output
    strng        -- string to file object opened for binary output
    """
    nchar = len(strng)
    fobj.write(struct.pack('i' + str(nchar) + 's',nchar, strng)) 

def _read_string(fobj, endian=''):
    """
    Reads a string written in binary format by my C++ code

    fobj   -- file object opened for binary input
    endian -- '>' for big-endian, '' for little-endian.
    """
    nchar  = struct.unpack(endian + 'i', fobj.read(4))[0]
    strng  = struct.unpack(endian + str(nchar) + 's', fobj.read(nchar))[0]
    return strng

def _check_ucm(fobj):
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
            raise UltracamError('ultracam._check_ucm: could not recognise first 4 bytes of ' + fname + ' as a ucm file')
        endian = '>'
    else:
        endian = ''
    return endian

# OK now onto publicly exposed stuff

class Odict(dict):
    """
    A dictionary which stores a key order which it uses for printing. This is a general
    class used as the base for the more specialised Uhead class for storing headers.
    """

    ninset  = 0
    nincr   = 5
    nlength = 20

    def __init__(self, dct = None):
        if dct == None:
            self._keys = []
            dict.__init__(self, {})
        else:
            self._keys = dct.keys()
            dict.__init__(self, dct)

    def __delitem__(self, key):
        dict.__delitem__(self, key)
        self._keys.remove(key)

    def __setitem__(self, key, item):
        dict.__setitem__(self, key, item)
        # 'try' needed to avoid error with pickling with protocol = 2 
        try:
            if key not in self._keys: self._keys.append(key)
        except AttributeError:
            self._keys = [key]

    def __str__(self):
        st    = ''
        inset = ' '*Odict.ninset
        Odict.ninset += Odict.nincr
        for key, val in self.iteritems():
            if isinstance(val, Odict):
                st += ('%s%-' + str(Odict.nlength) + 's\n') % (inset,key)
            else:
                st += ('%s%-' + str(Odict.nlength) + 's ') % (inset,key)
            st += str(val) + '\n'
            
        Odict.ninset -= Odict.nincr
        return st

    def clear(self):
        dict.clear(self)
        self._keys = []

    def copy(self):
        newInstance = Odict()
        newInstance.update(self)
        return newInstance

    def items(self):
        return zip(self._keys, self.values())

    def keys(self):
        return self._keys

    def popitem(self):
        try:
            key = self._keys[-1]
        except IndexError:
            raise KeyError('dictionary is empty')

        val = self[key]
        del self[key]

        return (key, val)

    def setdefault(self, key, failobj = None):
        if key not in self._keys: self._keys.append(key)
        return dict.setdefault(self, key, failobj)

    def update(self, dct):
        for key,val in dct.items():
            self.__setitem__(key,val)

    def values(self):
        return map(self.get, self._keys)

    def __iter__(self):
        for key in self._keys:
            yield key

    def iteritems(self):
        for key in self._keys:
            yield (key, self[key])

    def insert(self, key, item, index):
        """
        Adds a key, value pair just before the element with index = index, unless
        key already exists in which case its value is simply overwritten. Delete
        the element first if you want to re-order. index = 0 inserts
        at the start of the list.
        """
        dict.__setitem__(self, key, item)
        # 'try' needed to avoid error with pickling with protocol = 2 
        try:
            if key not in self._keys: self._keys.insert(index, key)
        except AttributeError:
            self._keys = [key]

class Window(np.ndarray):
    """
    Subclass of numpy.ndarray to represent a window of a CCD. The numpy part
    contains the data. The extra attributes are:
     llx, lly     -- lower-left pixels of window
     xbin, ybin   --  pixel binning factors

    Thus a Window is a numpy array that knows where it is in a detector.
    """

    def __new__(cls, data, llx, lly, xbin, ybin):
        """
        Creates  a Window given some data, a lower-left
        pixel position and binning factors.
        """

        obj = np.asarray(data).view(cls)

        # set the extra attributes
        obj.llx   = llx
        obj.lly   = lly
        obj.xbin  = xbin
        obj.ybin  = ybin

        # Return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # See the scipy web pages on subclassing ndarray objects
        if obj is None: return

        self.llx   = getattr(obj,'llx',0)
        self.lly   = getattr(obj,'lly',0)
        self.xbin  = getattr(obj,'xbin',0)
        self.ybin  = getattr(obj,'ybin',0)
        
    def __array_wrap__(self, obj):
        """
        This is to allow for the case when a single number is returned, e.g. by max()
        In this case rather than returning a Window of a zero-ranked array, we just
        return the number, losing the attributes
        """
        if len(obj.shape) == 0:
            return obj[()]
        else:
            return np.ndarray.__array_wrap__(self, obj)

    def __str__(self):
        ret  = '  llx, lly = ' + str(self.llx) + ', ' + str(self.lly) + \
            '; nx, ny = ' + str(self.shape[1]) + ', ' + str(self.shape[0]) + \
            '; xbin, ybin = ' + str(self.xbin) + ', ' + str(self.ybin) + '\n'
        ret += np.ndarray.__str__(self) + ', dtype = ' + str(self.dtype)
        return ret

    def __eq__(self, other):
        """
        Equality of two windows is defined by their matching binned dimensions
        lower-left pixel and binning-factors.
        """
        if type(other) is type(self):
            return (self.llx == other.llx and self.lly == other.lly and 
                    self.xbin == other.xbin and self.ybin == other.ybin and
                    self.shape == other.shape)
        else:
            return NotImplemented  

    def __neq__(self, other):
        """
        Negation of equality operator.
        """
        if type(other) is type(self):
            return (self.llx != other.llx or self.lly != other.lly or
                    self.xbin != other.xbin or self.ybin != other.ybin or
                    self.shape != other.shape)

        else:
            return NotImplemented  

    def plot(self, vmin, vmax, mpl=False, cmap=cm.binary):
        """
        Elementary intensity plot using either matplotlib's imshow 
        or pgplot's pggray. Typically some setup may be needed 
        before and after this.

        vmin   -- image value for lowest intensity
        vmax   -- image value for highest intensity
        mpl    -- True for matplotlib, otherwise pgplot
        cmap   -- colour map if mpl
        """
        ny, nx = self.shape
        if mpl:
            limits = self.llx-0.5,self.llx+self.xbin*nx-0.5,self.lly-0.5,self.lly+self.ybin*ny-0.5
            plt.imshow(np.asarray(self), cmap=cmap, interpolation='nearest', \
                           vmin=vmin, vmax=vmax, origin='lower', extent=limits)
        else:
            tr = [self.llx-self.xbin,self.xbin,0,self.lly-self.ybin,0,self.ybin]
            pg.pggray(self,0,nx-1,0,ny-1,vmax,vmin,tr)

class Time(object):
    """
    Represents a time for a CCD. Three attributes:

    mjd    -- modified Julian day number
    expose -- exposure time, seconds.
    good   -- is the time thought to be reliable?
    reason -- if good == False, this is the reason.
    """
    def __init__(self, mjd, expose, good, reason):
        self.mjd    = mjd
        self.expose = expose
        self.good   = good
        self.reason = reason

    def __str__(self):
        ret = 'MJD = ' + str(self.mjd) + ', exposure = ' + str(self.expose) + \
            ', status = ' + str(self.good)
        if not self.good:
            ret += ', reason: ' + self.reason
        return ret

    def __repr__(self):
        ret = '[' + str(self.mjd) + ' ' + str(self.expose) + \
            ' ' + str(self.good)
        if not self.good:
            ret += ' ' + self.reason + ']'
        else:
            ret += ']'
        return ret

# Integer type numbers for ucm files. Commented out ones
# are yet to be implemented
ITYPE_DOUBLE    =  0
#ITYPE_CHAR      =  1
ITYPE_INT       =  2
ITYPE_UINT      =  3
#ITYPE_LINT      =  4
#ITYPE_ULINT     =  5
ITYPE_FLOAT     =  6
ITYPE_STRING    =  7
ITYPE_BOOL      =  8
ITYPE_DIR       =  9
#ITYPE_DATE      = 10
ITYPE_TIME      = 11
#ITYPE_POSITION  = 12
ITYPE_DVECTOR   = 13
ITYPE_UCHAR     = 14
#ITYPE_TELESCOPE = 15
ITYPE_USINT     = 16
ITYPE_IVECTOR   = 17
ITYPE_FVECTOR   = 18

TNAME = {ITYPE_DOUBLE : 'double', ITYPE_INT : 'int', ITYPE_UINT : 'uint', 
         ITYPE_FLOAT : 'float', ITYPE_STRING : 'string', ITYPE_BOOL : 'bool',
         ITYPE_DIR : 'directory', ITYPE_TIME : 'time', ITYPE_DVECTOR : 'dvector',
         ITYPE_UCHAR : 'uchar', ITYPE_USINT : 'usint', ITYPE_IVECTOR : 'ivector',
         ITYPE_FVECTOR : 'fvector'}

# ucm magic number
MAGIC           = 47561009

class Uhead(Odict):
    """
    Class for containing headers compatible with ucm files. Each entry is
    keyed on a string of the form 'User.Filter', the dot signifying a
    hierarchy. The class is a sub-class of an ordered dictionary class to
    allow entries to retain an order. Each entry consists of a tuple
    containing a value, a type and a comment, in that order. The type
    corresponds to data types used in ucm files.
    
    The class is subclassed from Odict, a general ordered dictionary class
    supplied as part of the ultracam module. The purpose of subclassing this
    is to control access to the dictionary because of the special structure of
    the keys and values.
    """

    def __init__(self, head=None):
        """
        Constructor, either just a default () or with a header as a dictionary
        (preferably an ordered dictionary). This is tricky to construct
        properly, and it is probably easier to use several lines of
        'add_entry' statements. Whatever you supply will be passed to
        add_entry within a loop iterating over the elements of head.

        head -- a dictionary of key, value pairs, with keys having a
                hierarchical structure with '.'  separators and values each a
                tuple of (value,type,comment) -- see add_entry for more.
        """
        Odict.__init__(self)

    def __setitem__(self, key, value):
        raise UltracamError('Uhead.__setitem__ disabled to prevent invalid items being defined. Use add_entry')
        
    def add_entry(self, *args):
        """
        Adds a new Uhead item, checking the various arguments to reduce the
        chances of problems.  This can have either 2 or 4 arguments. The 4
        argument case is as follows:
        
        key   -- hierarchical string of the form 'User.Filter' where 'User' is a
                 directory or folder of grouped entries. It cannot have blanks
                 and any implied directories must already exists. Thus to set
                 a key 'User.Filter.Wheel', 'User.Filter' would need to exist
                 and be a directory. The existence of the implied 'User' would
                 not be checked in this case, on the assumption that it was
                 checked when 'User.Filter' was created.

        value -- value to associate (will be ignored in the case of
                 directories, but see the 2 argument case below). The nature
                 of the value varies with the itype; see next.

        itype -- one of a range of possible data types. This rather
                 'unpythonic' argument is to address the need to match up with
                 data files and the C++ ULTRACAM pipeline when it comes to
                 writing to disk. Look for integers called 'ITYPE_*' to see
                 the set of potential types. The meaning of most data types is
                 obvious. e.g.  ITYPE_DOUBLE or ITYPE_FLOAT expect floating
                 point numbers. In this case both will be stored as a Python
                 float in memory, but will be saved to disk with different
                 numbers of bytes. Less obvious ones are:

                 ITYPE_TIME -- the corresponding value should be a two-element
                               tuple or list with first an integer for the
                               number of days and then a float for the
                               fraction of a day.


        comment -- comment string with details of the variable.

        If just 2 arguments are given, they will be interpreted as just a key
        and comment for a directory.
        """

        # elementary checks
        if len(args) == 2:
            key, comment = args
            itype = ITYPE_DIR
            value = None
        elif len(args) == 4:
            key, value, itype, comment = args
        else:
            raise UltracamError('Uhead.add_entry: takes either 2 or 4 arguments')

        if not isinstance(key, basestring):
            raise UltracamError('Uhead.add_entry: argument "key" must be a string.')

        if not isinstance(comment, basestring):
            raise UltracamError('Uhead.add_entry: key = ' + key + ': "comment" must be a string.')

        # now look at the key: must have no blanks
        if key.find(' ') > -1:
            raise UltracamError('Uhead.add_entry: key = "' + key + '" contains at least one blank.')

        # if key has a '.' then the part before last dot must already exist
        # and must be a directory. Search in reverse order, as all being well, it
        # should usually be fastest.
        ldot = key.rfind('.')
        if ldot > -1:
            dir = key[:ldot]
            for kold in self.keys()[::-1]:
                if dir == kold and self[kold][1] == ITYPE_DIR:
                    break
            else:
                raise UltracamError('Uhead.add_entry: key = ' + key + 
                                    ': could not locate directory = ' + key[:ldot])

            # determine position of key within Odict. Must add onto 
            # whichever directory it belongs to.
            for index, kold in enumerate(self.keys()):
                if kold.startswith(dir): lind = index

        # the next implicitly check the value: if they can't be converted to the right type,
        # something is wrong.
        if itype == ITYPE_DOUBLE or itype == ITYPE_FLOAT:
            value = float(value)
        elif itype == ITYPE_INT or itype == ITYPE_UINT or \
                itype == ITYPE_UCHAR or itype == ITYPE_USINT:
            value = int(value)
        elif itype == ITYPE_STRING:
            value = str(value)
        elif itype == ITYPE_BOOL:
            value = bool(value)
        elif itype == ITYPE_DIR:
            pass
        elif itype == ITYPE_TIME:
            if len(value) != 2:
                raise UltracamError('Uhead.add_entry: key = ' + key + 
                                ': require a 2-element tuple or list (int,float) for ITYPE_TIME)')
            value[0] = int(value[0])
            value[1] = float(value[1])
        elif itype == ITYPE_DVECTOR:
            if not isinstance(value, np.ndarray) or len(value.shape) != 1:
                raise UltracamError('Uhead.add_entry: key = ' + key + 
                                ': require a 1D numpy.ndarray for ITYPE_DVECTOR)')
            value = value.astype(float64)
        elif itype == ITYPE_IVECTOR:
            if not isinstance(value, np.ndarray) or len(value.shape) != 1:
                raise UltracamError('Uhead.add_entry: key = ' + key + 
                                ': require a 1D numpy.ndarray for ITYPE_IVECTOR)')
            value = value.astype(int)
        elif itype == ITYPE_FVECTOR:
            if not isinstance(value, np.ndarray) or len(value.shape) != 1:
                raise UltracamError('Uhead.add_entry: key = ' + key + 
                                ': require a 1D numpy.ndarray for ITYPE_FVECTOR)')
            value = value.astype(float32)
        else:
            raise UltracamError('Uhead.add_entry: key = ' + key + 
                            ': itype = ' + str(itype) + ' not recognised.')

        # checks passed, finally set item
        if ldot > -1:
            self.insert(key, (value, itype, comment), lind+1)
        else:
            Odict.__setitem__(self, key, (value, itype, comment))

    def value(self, key):
        "Returns the value associated with a given key"
        return self[key][0]

    def itype(self, key):
        "Returns the itype associated with a given key"
        return self[key][1]

    def comment(self, key):
        "Returns the comment associated with a given key"
        return self[key][2]

    def __str__(self):
        ret = ''
        for key, val in self.iteritems():
            ndot  = key.count('.')
            final = key[key.rfind('.')+1:]
            if val[1] == ITYPE_DIR:
                ret  += '\n' + ndot*'  '
                ret  += '%-20s /directory/   %s\n' % (final,val[2])
            else:
                ret  += ndot*'  '
                ret  += '%-20s = %-20s %-12s   %s\n' % \
                    (final,str(val[0]),'/'+TNAME[val[1]]+'/',val[2]) 
        return ret

class CCD(list):
    """
    Class to represent a CCD. Sub-class of 'list'. The list in this
    case is a list of Windows. Added to this are attributes to 
    represent the maximum unbinned dimensions of the CCD, a time
    a flag to say whether the data are good, and a header.
    """
    def __init__(self, wins, time, nxmax, nymax, good, head):
        """
        Creates a new CCD frame.
        
        Arguments:

        wins    -- list of non-overlapping Window objects.
        time    -- a Time representing the central time of the CCD.
        nxmax   -- maximum dimension in X, unbinned pixels.
        nymax   -- maximum dimension in Y, unbinned pixels.
        good    -- True if data are not junk.
        head    -- header. Must be None or a 'Uhead'
        """
        for win in wins:
            if not isinstance(win, Window):
                raise UltracamError('CCD.__init__: all windows must be Windows.')        

        if head is not None and not isinstance(head, Uhead):
            raise UltracamError('CCD.__init__: head should be a Uhead (or None).')        

        list.__init__(self, wins)
        self.time  = time
        self.nxmax = nxmax
        self.nymax = nymax
        self.good  = good
        self.head  = head


    def __eq__(self, other):
        """
        Equality of two CCDs is defined by matching binning factors, 
        maximum dimensions and windows (in order).
        """
        if type(other) is type(self):
            if self.nxmax != other.nymax or len(self) != len(other): return False
        
            # now test for equal windows
            for swin,owin in zip(self,other):
                if swin != owin: return False
            return True
        else:
            return NotImplemented

    def __ne__(self, other):
        """
        Negation of equality operator.
        """
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def anyInt(self):
        """
        Returns True if any of the contributing Windows are based on integers. It can be
        useful for memory and disk space reasons to keep data as 2-byte unsigned integers 
        but cause problems with arithematic operations. This allows you to check. 
        See also 'anyFloat', 'toFloat' and 'toInt'.
        """
        for win in self:
            if issubclass(win.dtype.type,np.integer):
                return True
        return False

    def anyFloat(self):
        """
        Returns True if any of the contributing Windows are based on floats. This is needed
        to evaluate the output type needed when writing to disk.
        """
        for win in self:
            if issubclass(win.dtype.type,np.floating):
                return True
        return False

    def toFloat(self, single=True):
        """
        Converts all Windows to a float type, either single or double
        precision.

        single  -- True to convert to 4-byte floats (else 8-byte)
        """
        for win in self:
            if single:
                win = win.astype(np.float32)
            else:
                win = win.astype(np.float64)

    def toInt(self):
        """
        Converts all Windows to an unsigned 2-byte integer type, rounding
        to the nearest integer. Warnings will be issued if data lies outside
        the 0 to 65535 range, but the conversion will proceed.
        """
        for win in self:
            if win.min() < 0 or win.max() > 65535:
                warnings.warn('CCD.toInt: input data out of range 0 to 65535')
            win = np.rint(win).astype(np.uint16)
        
    def mean(self):
        """
        Returns the mean over all Windows of a CCD
        """
        nelem = 0
        sum   = 0.
        for win in self:
            nelem += win.size
            sum   += win.sum()
        return sum / float(nelem)

    def min(self):
        """
        Returns the minimum over all Windows of a CCD
        """
        minv = None
        for win in self:
            minv = win.min() if minv is None else min(minv, win.min())
        return minv

    def max(self):
        """
        Returns the maximum over all Windows of a CCD
        """
        maxv = None
        for win in self:
            maxv = win.max() if maxv is None else max(maxv, win.max())
        return maxv

    def median(self):
        """
        Returns median over all Windows of a CCD. 
        """

        # generate combined list of all pixels in CCD called 'arr'
        larr = []
        for win in self:
            larr.append(win.flatten())
        arr = np.concatenate(larr)

        return np.median(arr)


    def centile(self, pcent):
        """
        Returns percentile(s) over all Windows of a CCD. Given pcent, this
        routine returns the image level below which pcent percent of the pixel
        values lie.  pcent can be a single number or array-like. In the latter
        case a list of values is returned.

        pcent -- percentile or percentiles (array-like)

        Returns image value or values as a list equivalent to the input 
        percentiles. 
        """

        # check against a string which can look array-like
        if isinstance(pcent, basestring):
            raise UltracamError('CCD.centile: argument "pcent" cannot be a string')

        # generate combined list of all pixels in CCD called 'arr'
        larr = []
        for win in self:
            larr.append(win.flatten())
        arr = np.concatenate(larr)

        return np.percentile(arr,pcent)
    
    def rback(self):
        """
        Removes background from a CCD. Estimates
        background using a median of each window
        separately.
        """
        for win in self:
            win -= np.median(win)

    def plot(self, vmin, vmax, mpl=False, cmap=cm.binary):
        """
        Elementary intensity plot using either matplotlib's imshow 
        or pgplot's pggray. Typically some setup may be needed 
        before and after this. This one simply plots all windows
        using the Window.plot method

        vmin -- image value for lowest intensity
        vmax -- image value for highest intensity
        mpl  -- True for matplotlib, otherwise pgplot
        cmap -- colour map if mpl
        """
        for win in self:
            win.plot(vmin,vmax,mpl,cmap)

    def append(self, item):
        raise NotImplementedError

    # arithematic
    def __iadd__(self, other):
        """
        Adds 'other' to the CCD in place (+=)
        """
        for win,owin in zip(self,other):
            win += owin
        return self

    def __isub__(self, other):
        """
        Subtracts 'other' from the CCD in place (-=)
        """
        for win,owin in zip(self,other):
            win -= owin
        return self

    def __imul__(self, other):
        """
        Multiplies the CCD by 'other' in place (*=)
        """
        for win,owin in zip(self,other):
            win *= owin
        return self

    def __idiv__(self, other):
        """
        Divides the CCD by 'other' in place (/=)
        """
        for win,owin in zip(self,other):
            win /= owin
        return self

    def __add__(self, other):
        """
        Adds 'other' to the CCD (+)
        """
        twins = []
        for win,owin in zip(self,other):
            twins.append(win + owin)
        return CCD(twins, self.time, self.nxmax, self.nymax, self.good and other.good, self.head)

    def __sub__(self, other):
        """
        Subtracts 'other' from the CCD (-)
        """
        twins = []
        for win,owin in zip(self,other):
            twins.append(win - owin)
        return CCD(twins, self.time, self.nxmax, self.nymax, self.good and other.good, self.head)

    def __mul__(self, other):
        """
        Multiplies CCD by 'other' (*)
        """
        twins = []
        for win,owin in zip(self,other):
            twins.append(win * owin)
        return CCD(twins, self.time, self.nxmax, self.nymax, self.good and other.good, self.head)

    def __div__(self, other):
        """
        Divides CCD by 'other' (/)
        """
        twins = []
        for win,owin in zip(self,other):
            twins.append(win / owin)
        return CCD(twins, self.time, self.nxmax, self.nymax, self.good and other.good, self.head)

    def __str__(self):
        """
        Generates readable summary of a CCD
        """

        ret = ''
        if self.head is not None: ret += str(self.head)

        ret += '\nDimensions = ' + str(self.nxmax) + ', ' + str(self.nymax) + \
            ', number of windows = ' + str(len(self)) + ', status = ' + str(self.good) + '\n'

        for nwin,win in enumerate(self):
            ret += '\nWindow number ' + str(nwin+1) + ':\n'
            ret += str(win) + '\n'
        return ret

class MCCD(list):
    """
    Represents multiple CCD frame. The idea is that one has an instrument
    which produces multiple CCDs of data that are intrinsically linked, 
    e.g. in a multi-arm camera where one always gets an image in each of 
    several CCDs. This is to represent ULTRACAM data. There is some common 
    header information, plus the data, represented by a list of CCD objects.

    To allow easy access to the data, like CCD, MCCD is a subclass of list 
    so that a construction like mccd[1] returns the second CCD, mccd[1][0]
    returns the first Window of the second CCD and 

    for ccd in mccd:
      print 'number of windows = ',len(ccd)

    steps through each CCD printing the number of windows.
    """

    def __init__(self, data, head):
        """
        Creates an MCCD.

        Arguments:

          data  -- list of CCD objects

          head -- the header, either None of a Uhead object.
          
        Sets the equivalent attribute head
        """
        if not isinstance(data, list):
            raise UltracamError('MCCC.__init__: data should be a list.')

        for ccd in data:
            if not isinstance(ccd, CCD):
                raise UltracamError('MCCC.__init__: one or more of the elements of data is not a CCD.')

        if head is not None and not isinstance(head, Uhead):
            raise UltracamError('MCCC.__init__: head should be a Uhead (or None).')

        list.__init__(self, data)
        self.head  = head

    @classmethod
    def rucm(cls, fname, flt=True):
        """
        Factory method to produce an MCCD from a ucm file.

        fname -- ucm file name. '.ucm' will be appended if not supplied.

        flt    -- convert to 4-byte floats whatever the input data, or not. ucm
                  files can either contain 4-bytes floats or for reduced disk
                  footprint, unsigned 2-byte integers. If flt=True, either type
                  will end up as float32 internally. If flt=False, the disk type
                  will be retained. The latter is unsafe when arithematic is involved
                  hence the default is to convert to 4-byte floats.

        Exceptions are thrown if the file cannot be found, or an error during the
        read occurs.
        """    

        # Assume it is a file object, if that fails, assume it is
        # the name of a file.
        if not fname.endswith('.ucm'): fname += '.ucm'
        uf = open(fname, 'rb')
        start_format =  _check_ucm(uf)

        # read the header
        lmap = struct.unpack(start_format + 'i', uf.read(4))[0]

        head = Uhead()
        for i in xrange(lmap):
            name    = _read_string(uf, start_format)
            itype   = struct.unpack(start_format + 'i', uf.read(4))[0]
            comment = _read_string(uf, start_format)

            if itype == ITYPE_DOUBLE:
                value = struct.unpack(start_format + 'd', uf.read(8))[0]
            elif itype == ITYPE_INT:
                value = struct.unpack(start_format + 'i', uf.read(4))[0]
            elif itype == ITYPE_UINT:
                value = struct.unpack(start_format + 'I', uf.read(4))[0]
            elif itype == ITYPE_FLOAT:
                value = struct.unpack(start_format + 'f', uf.read(4))[0]
            elif itype == ITYPE_STRING:
                value = _read_string(uf, start_format)
            elif itype == ITYPE_BOOL:
                value = struct.unpack(start_format + 'B', uf.read(1))[0]
            elif itype == ITYPE_DIR:
                value = None
            elif itype == ITYPE_TIME:
                value = struct.unpack(start_format + 'id', uf.read(12))
            elif itype == ITYPE_DVECTOR:
                nvec  = struct.unpack(start_format + 'i', uf.read(4))[0]
                value = struct.unpack(start_format + str(nvec) + 'd', uf.read(8*nvec))
            elif itype == ITYPE_UCHAR:
                value = struct.unpack(start_format + 'c', uf.read(1))[0]
            elif itype == ITYPE_USINT:
                value = struct.unpack(start_format + 'H', uf.read(2))
            elif itype == ITYPE_IVECTOR:
                nvec  = struct.unpack(start_format + 'i', uf.read(4))[0]
                value = struct.unpack(start_format + str(nvec) + 'i', uf.read(4*nvec))
            elif itype == ITYPE_FVECTOR:
                nvec  = struct.unpack(start_format + 'i', uf.read(4))[0]
                value = struct.unpack(start_format + str(nvec) + 'f', uf.read(4*nvec))
            else:
                raise UltracamError('ultracam.MCCD.rucm: do not recognize itype = ' + str(itype))

            # store header information, fast method
            Odict.__setitem__(head, name, (value, itype, comment))
        
        # now for the data
        data  = []
        
        # read number of CCDs
        nccd = struct.unpack(start_format + 'i', uf.read(4))[0]

        for nc in range(nccd):
            # read number of wndows
            nwin = struct.unpack(start_format + 'i', uf.read(4))[0]
            wins  = []
            for nw in range(nwin):
                llx,lly,nx,ny,xbin,ybin,nxmax,nymax,iout = struct.unpack(start_format + '9i', uf.read(36))
                if iout == 0:
                    win = np.fromfile(file=uf, dtype=np.float32, count=nx*ny)
                elif iout == 1:
                    if flt:
                        win = np.fromfile(file=uf, dtype=np.uint16, count=nx*ny).astype(np.float32)
                    else:
                        win = np.fromfile(file=uf, dtype=np.uint16, count=nx*ny)
                else:
                    raise UltracamError('ultracam.MCCD.rucm: iout = ' + str(iout) + ' not recognised')
                win = win.reshape((ny,nx))
                wins.append(Window(win,llx,lly,xbin,ybin))

            data.append(CCD(wins,None,nxmax,nymax,True,None))
        uf.close()

        return cls(data, head)

    def __eq__(self, other):
        """
        Equality operator tests same number of CCDs and that each CCD matches.
        """

        if type(other) is type(self):

            if len(self) != len(other): return False

            for sccd, occd in zip(self,other):
                if len(sccd) != len(occd): return False
            return True
        else:
            return NotImplemented

    def __ne__(self, other):
        """
        Inequality operator based on file formats: same number of CCDs, same
        number of windows per CCD, same binning factors etc.
        """
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def wucm(self, fname):
        """
        Writes to disk in ucm format. The data are saved as 32-bit floats
        or 16-bit unsigned integers according to the internal types.

        fname  -- file to write to. '.ucm' will be appended if necessary.
        """    

        if not fname.strip().endswith('.ucm'):
            fname = fname.strip() + '.ucm'
        uf = open(fname, 'wb')
    
        # write the format code
        uf.write(struct.pack('i',MAGIC))

        # write the header, starting with the number of entries
        lmap = len(self.head)
        uf.write(struct.pack('i',lmap))

        for key,val in self.head.iteritems():

            _write_string(uf, key)

            value, itype, comment = val
            uf.write(struct.pack('i',itype))
            _write_string(uf, comment)

            if itype == ITYPE_DOUBLE:
                uf.write(struct.pack('d', value))
            elif itype == ITYPE_INT:
                uf.write(struct.pack('i', value))
            elif itype == ITYPE_UINT:
                uf.write(struct.pack('I', value))
            elif itype == ITYPE_FLOAT:
                uf.write(struct.pack('f', value))
            elif itype == ITYPE_STRING:
                _write_string(uf, value)
            elif itype == ITYPE_BOOL:
                uf.write(struct.pack('B', value))
            elif itype == ITYPE_DIR:
                pass
            elif itype == ITYPE_TIME:
                uf.write(struct.pack('i', value[0]))
                uf.write(struct.pack('d', value[1]))
            elif itype == ITYPE_DVECTOR:
                uf.write(struct.pack('i', len(value)))
                uf.write(struct.pack(str(len(value))+'d', *value))
            elif itype == ITYPE_UCHAR:
                uf.write(struct.pack('B', value))
            elif itype == ITYPE_USINT:
                uf.write(struct.pack('H', value))
            elif itype == ITYPE_IVECTOR:
                uf.write(struct.pack('i', len(value)))
                uf.write(struct.pack(str(len(value))+'i', *value))
            elif itype == ITYPE_FVECTOR:
                uf.write(struct.pack('i', len(value)))
                uf.write(struct.pack(str(len(value))+'f', *value))
            else:
                raise UltracamError('Hitem: type =' + str(itype) + 'not recognised')

        # number of CCDs
        nccd = len(self)
        uf.write(struct.pack('i', nccd))

        iout = 0 if self.anyFloat() else 1

        for ccd in self:

            # number of windows
            nwin = len(ccd)
            uf.write(struct.pack('i', nwin))

            for win in ccd:
                ny,nx   = win.shape
                uf.write(struct.pack('9i',win.llx,win.lly,nx,ny,win.xbin,win.ybin,
                                     ccd.nxmax,ccd.nymax,iout))
                if iout == 0:
                    win.astype(np.float32).tofile(uf)
                elif iout == 1:
                    if issubclass(win.dtype.type,np.integer):
                        win.astype(np.uint16).tofile(uf)
                    else:
                        np.rint(win).astype(np.uint16).tofile(uf)

        uf.close()

    def rback(self, nc=-1):
        """
        Removes background from one or all CCDs from an MCCD.

        nc  -- the CCD to remove the background from. -1 for all.
        """
        if nc == -1:
            for ccd in self:
                ccd.rback()
        else:
            self.data[nc].rback()

    def anyInt(self):
        """
        Returns True if any of the contributing CCDs are based on integers. It can be
        useful for memory and disk space reasons to keep data as 2-byte unsigned integers 
        but cause problems with arithematic operations. This allows you to check. 
        See also 'anyFloat', 'toFloat' and 'toInt'.
        """
        for ccd in self:
            if ccd.anyInt(): return True
        return False

    def anyFloat(self):
        """
        Returns True if any of the contributing CCDs are based on floats. This is needed
        to evaluate the output type needed when writing to disk.
        """
        for ccd in self:
            if ccd.anyFloat(): return True
        return False

    def toFloat(self, single=True):
        """
        Converts all CCDs to a float type, either single or double
        precision.

        single  -- True to convert to 4-byte floats (else 8-byte)
        """
        for ccd in self:
            ccd.toFloat(single)

    def toInt(self):
        """
        Converts all CCDs to an unsigned 2-byte integer type, rounding
        to the nearest integer. Warnings will be issued if data lies outside
        the 0 to 65535 range, but the conversion will proceed.
        """
        for ccd in self:
            ccd.toInt()

    def max(self):
        """
        Returns a tuple of maximum values, 1 per CCD.
        """
        mx = []
        for ccd in self:
            mx.append(ccd.max())
        return tuple(mx)

    def min(self):
        """
        Returns a tuple of minimum values, 1 per CCD.
        """
        mn = []
        for ccd in self:
            mn.append(ccd.min())
        return tuple(mn)

    def mean(self):
        """
        Returns a tuple of mean values, 1 per CCD.
        """
        mn = []
        for ccd in self:
            mn.append(ccd.mean())
        return tuple(mn)

    def median(self):
        """
        Returns a tuple of median values, 1 per CCD.
        """
        mn = []
        for ccd in self:
            mn.append(ccd.median())
        return tuple(mn)

    # arithematic
    def __iadd__(self, other):
        """
        Adds 'other' to the MCCD in place (+=)
        """
        for ccd,occd in zip(self,other):
            ccd += occd
        return self

    def __isub__(self, other):
        """
        Subtracts 'other' from the MCCD in place (-=)
        """
        for ccd,occd in zip(self,other):
            ccd -= occd
        return self

    def __imul__(self, other):
        """
        Multiplies the MCCD by 'other' in place (*=)
        """
        for ccd,occd in zip(self,other):
            ccd *= occd
        return self

    def __idiv__(self, other):
        """
        Divides the MCCD by 'other' in place (/=)
        """
        for ccd,occd in zip(self,other):
            ccd /= occd
        return self

    def __str__(self):
        ret = ''
        if self.head is not None: ret += str(self.head)
        ret += '\n\nNumber of CCDs = ' + str(len(self)) + '\n'

        for nccd,ccd in enumerate(self):
            ret += '\nCCD number ' + str(nccd+1) + ':\n'
            ret += str(ccd)
        return ret

    def __add__(self, other):
        """
        Adds 'other' to the MCCD (+)
        """
        tccd = []
        for ccd,occd in zip(self,other):
            tccd.append(ccd + occd)
        return MCCD(tccd, self.head)

    def __sub__(self, other):
        """
        Subtract 'other' from the MCCD (-)
        """
        tccd = []
        for ccd,occd in zip(self,other):
            tccd.append(ccd - occd)
        return MCCD(tccd, self.head)

    def __mul__(self, other):
        """
        Multiply 'other' by the MCCD (*)
        """
        tccd = []
        for ccd,occd in zip(self,other):
            tccd.append(ccd * occd)
        return MCCD(tccd, self.head)

    def __div__(self, other):
        """
        Divide MCCD by 'other' from the MCCD (/)
        """
        tccd = []
        for ccd,occd in zip(self,other):
            tccd.append(ccd / occd)
        return MCCD(tccd, self.head)

    def plot(self, vlo=2., vhi=98., nc=-1, method='p', mpl=False, cmap=cm.binary, \
                 close=True, x1=None, x2=None, y1=None, y2=None):
        """
        Plots an MCCD using pgplot or matplotlib if preferred.

         vlo    -- number specifying the lowest level to plot (default as a percentile)
         vhi    -- number specifying the lowest level to plot (default as a percentile)
         nc     -- CCD number (starting from 0, -1 for all)
         method -- how vlo and vhi are to be interpreted. 'p' = percentile, 'a' = automatic (min to max,
                   vlo and vhi are irrelevant), 'd' = direct, i.e. just take the values given.
         mpl    -- True to prefer matplotlib over pgplot (which may not even be an option)
         cmap   -- colour map if using matplotlib
         close  -- close (pgplot) or 'show' (matplotlib) the plot at the end (or not, to allow 
                   you to plot something else, use a cursor etc). In the case of pgplot, this also
                   implies opening the plot at the start, i.e. a self-contained quick plot.
         x1,x2, -- plot limits. Default to 0.5,nxmax+0.5,0.5,nymax+0.5 if not defined. 
         y1,y1 

        Returns the plot range(s) used either as a single 2-element tuple, or
        a list of them, one per CCD plotted.
        """

        if nc == -1:
            nc1 = 0
            nc2 = len(self)
        else:
            nc1 = nc
            nc2 = nc+1
        
        if not mpl:
            if close: pg.pgopen('/xs')
            pg.pgsubp(nc2-nc1,1)

        prange = []
        for nc, ccd in enumerate(self[nc1:nc2]):

            # Determine intensity range to display
            if method == 'p':
                vmin, vmax = ccd.centile((vlo,vhi))
            elif method == 'a':
                vmin, vmax = ccd.min(), ccd.max()
            elif method == 'd':
                vmin, vmax = vlo, vhi
            else:
                raise UltracamError('MCCD.plot: method must be one of p, a or d.')
            
            prange.append((vmin, vmax))

            # start
            nxmax, nymax = ccd.nxmax, ccd.nymax
            x1    = 0.5 if x1 is None else x1
            x2    = nxmax+0.5 if x2 is None else x2
            y1    = 0.5 if y1 is None else y1
            y2    = nymax+0.5 if y2 is None else y2

            if mpl:
                plt.add_subplot(1,nc2-nc1,nc-nc1+1)
                plt.axis('equal')
            else:
                pg.pgpanl(nc-nc1+1,1)
                pg.pgwnad(x1,x2,y1,y2)

            # plot CCD
            ccd.plot(vmin,vmax,mpl,cmap)

            # per-frame finishing-off
            if mpl:
                plt.xlim(x1,x2)
                plt.ylim(y1,y2)
            else:
                pg.pgbox('bcnst',0,0,'bcnst',0,0)
                pg.pglab('X','Y','')

        if close:
            if mpl:
                plt.show()
            else:
                pg.pgclos()

        # return intensity range(s) used
        if len(prange) == 1:
            return prange[0]
        else:
            return tuple(prange)

class _win(object):
    """
    Trivial container class for basic window info
    to help readability of code
    """
    def __init__(self, llx, lly, nx, ny):
        self.llx = llx
        self.lly = lly
        self.nx  = nx
        self.ny  = ny

class Rhead (object):
    """
    Represents essential header info of Ultracam/Ultraspec data read from a
    run###.xml file. 
    """

    def __init__(self, uxml):
        """
        Reads a run###.xml file. UltracamErrors are thrown if some items are not found.
        In some case it will carry on and corresponding attributes are returned as None.

        Arguments:

         uxml     -- xml file name with format run###.xml

         Many attributes are set; here are some of them:

         application  -- data acqusition application template name.

         fname        -- file used to define the format.

         framesize    -- total number of bytes per frame.

         headerwords  -- number of words (2-bytes/word) in timing info at start of a frame.

         instrument   -- instrument name.

         mode         -- a standardised summary of the readout mode derived from the application name.

         speed        -- readout speed.

         user         -- dictionary of user information. Set to None if there was none found.

         win          -- A list of Window objects, one per window. ULTRACAM data is multi-CCD
                         but the windows of each CCD are identical so the information is only stored 
                         once for all CCDs. Each one contains the corrdinates of the lower-left
                         pixel of the window and the binned dimensions

         xbin, ybin   -- pixel binning factors (same for all windows of all CCDs)

         nxmax, nymax -- maximum unbinned dimensions, same for all CCDs.
        """

        self.fname = uxml
        udom = xml.dom.minidom.parse(uxml)

        # Find framesize and headerwords.
        node        = udom.getElementsByTagName('data_status')[0]
        self.framesize   = int(node.getAttribute('framesize'))
        self.headerwords = int(node.getElementsByTagName('header_status')[0].getAttribute('headerwords'))

        # Frame format and other detail.
        node   = udom.getElementsByTagName('instrument_status')[0]
        self.instrument  = node.getElementsByTagName('name')[0].childNodes[0].data
        if self.instrument == 'Ultracam':
            self.instrument = 'ULTRACAM'
            self.nxmax, self.nymax = 1080, 1032
        elif self.instrument == 'Ultraspec':
            self.instrument = 'ULTRASPEC'
            self.nxmax, self.nymax = 1056, 1072
        else:
            raise UltracamError('Rhead.__init__: file = ' + self.fname + ', failed to identify instrument.')

        self.application = [nd for nd in node.getElementsByTagName('application_status') \
                                if nd.getAttribute('id') == 'SDSU Exec'][0].getAttribute('name')

        # get user info, if present
        try:
            nlist = dom.getElementsByTagName('user')
            if len(nlist):
                user = {}
                node = nlist[0]
                for nd in node.childNodes:
                    if nd.nodeType == Node.ELEMENT_NODE and nd.hasChildNodes():
                        user[nd.tagName] = nd.childNodes[0].data
            else:
                user = None
        except Exception, err:
            user = None

        # gather together majority of values
        param = {}
        for nd in node.getElementsByTagName('parameter_status'):
            param[nd.getAttribute('name')] = nd.getAttribute('value')

        try:
            # look for user-supplied information
            nlist = udom.getElementsByTagName('user')
            if len(nlist):
                self.user = {}
                node = nlist[0]
                for nd in node.childNodes:
                    if nd.nodeType == Node.ELEMENT_NODE and nd.hasChildNodes():
                        self.user[nd.tagName] = nd.childNodes[0].data
        except Exception, err:
            self.user = None

        # Translate applications into meaningful mode names
        app = self.application
        if app == 'ap8_250_driftscan' or app == 'ap8_driftscan' or app == 'ap_drift_bin2' or \
                app == 'appl8_driftscan_cfg':
            self.mode    = 'DRIFT'
        elif app == 'ap5_250_window1pair' or app == 'ap5_window1pair' or app == 'ap_win2_bin8' or \
                app == 'ap_win2_bin2' or app == 'appl5_window1pair_cfg':
            self.mode    = '1-PAIR'
        elif app == 'ap5b_250_window1pair' or app == 'appl5b_window1pair_cfg':
            self.mode    = '1-PCLR'
        elif app == 'ap6_250_window2pair' or app == 'ap6_window2pair' or \
                app == 'ap_win4_bin1' or app == 'ap_win4_bin8' or app == 'appl6_window2pair_cfg':
            self.mode    = '2-PAIR'
        elif app == 'ap7_250_window3pair' or app == 'ap7_window3pair' or app == 'appl7_window3pair_cfg':
            self.mode    = '3-PAIR'
        elif app == 'ap3_250_fullframe' or app == 'ap3_fullframe' or app == 'appl3_fullframe_cfg':
            self.mode    = 'FFCLR'
        elif app == 'appl4_frameover_cfg' or app == 'ap4_frameover':
            self.mode    = 'FFOVER'
        elif app == 'ap9_250_fullframe_mindead' or app == 'ap9_fullframe_mindead' or \
                app == 'appl9_fullframe_mindead_cfg':
            self.mode    = 'FFNCLR'
        elif app == 'ccd201_winbin_con' or app == 'ccd201_winbin_cfg':
            if int(param['X2_SIZE']) == 0:
                self.mode    = '1-USPEC'
            else:
                self.mode    = '2-USPEC'
        elif app == 'ccd201_driftscan_cfg':
            self.mode    = 'UDRIFT'
        elif app == 'ap1_poweron' or app == 'ap1_250_poweron' or app == 'ap2_250_poweroff':
            self.mode = 'PONOFF'
            return
        else:
            raise UltracamError('Rhead.__init__: file = ' + self.fname + ' failed to identify application = ' + app)

        # binning factors
        self.xbin = int(param['X_BIN_FAC']) if 'X_BIN_FAC' in param else int(param['X_BIN'])
        self.ybin = int(param['Y_BIN_FAC']) if 'Y_BIN_FAC' in param else int(param['Y_BIN'])

        # Windows. For each one store: x & y coords of lower-left pixel, binned dimensions
        self.win = []
        fsize = 2*self.headerwords
        if self.instrument == 'ULTRACAM':

            self.exposeTime = float(param['EXPOSE_TIME'])
            self.numexp     = int(param['NO_EXPOSURES'])
            self.gainSpeed  = hex(int(param['GAIN_SPEED']))[2:] if 'GAIN_SPEED' in param else None

            if 'V_FT_CLK' in param:
                self.v_ft_clk  = struct.unpack('<B',param['V_FT_CLK'][2])[0]
            elif app == 'appl7_window3pair_cfg':
                self.v_ft_clk  = 140;
            else:
                self.v_ft_clk = 0

            self.nblue    = int(param['NBLUE']) if 'NBLUE' in param else 1

            if self.mode == 'FFCLR' or self.mode == 'FFNCLR':
                self.win.append(_win(  1, 1, 512//self.xbin, 1024//self.ybin))
                self.win.append(_win(513, 1, 512//self.xbin, 1024//self.ybin))
            elif self.mode == 'FFOVER':
                self.win.append(_win(  1, 1, 540//self.xbin, 1032//self.ybin))
                self.win.append(_win(541, 1, 540//self.xbin, 1032//self.ybin))
            else:
                ystart = int(param['Y1_START'])
                xleft  = int(param['X1L_START'])
                xright = int(param['X1R_START'])
                nx     = int(param['X1_SIZE']) // self.xbin
                ny     = int(param['Y1_SIZE']) // self.ybin
                self.win.append(_win(xleft, ystart, nx, ny))
                self.win.append(_win(xright, ystart, nx, ny))
            
            fsize += 12*self.win[-1].nx*self.win[-1].ny

            if self.mode == '2-PAIR' or self.mode == '3-PAIR':
                ystart = int(param['Y2_START'])
                xleft  = int(param['X2L_START'])
                xright = int(param['X2R_START'])
                nx     = int(param['X2_SIZE']) // self.xbin
                ny     = int(param['Y2_SIZE']) // self.ybin
                self.win.append(_win(xleft, ystart, nx, ny))
                self.win.append(_win(xright, ystart, nx, ny))
                fsize += 12*self.win[-1].nx*self.win[-1].ny

            if self.mode == '3-PAIR':
                ystart = int(param['Y3_START'])
                xleft  = int(param['X3L_START'])
                xright = int(param['X3R_START'])
                nx     = int(param['X3_SIZE']) // self.xbin
                ny     = int(param['Y3_SIZE']) // self.ybin
                self.win.append(_win(xleft,ystart,nx,ny))
                self.win.append(_win(xright,ystart,nx,ny))
                fsize += 12*self.win[-1].nx*self.win[-1].ny

        elif self.instrument == 'ULTRASPEC':

            self.exposeTime   = float(param['DWELL'])
            self.numexp   = int(param['NUM_EXPS'])
            self.speed    = ('F' if param['SPEED'] == '0' else \
                                 ('M' if param['SPEED'] == '1' else 'S')) if 'SPEED' in param else None
            self.en_clr   = ('Y' if param['EN_CLR'] == '1' else 'N') if 'EN_CLR' in param else None
            self.hv_gain  = param['HV_GAIN'] if 'HV_GAIN' in param else None
            self.output   = ('N' if param['OUTPUT'] == '0' else 'A') if 'OUTPUT' in param else None

            xstart = int(param['X1_START'])
            ystart = int(param['Y1_START'])
            nx     = int(param['X1_SIZE']) // self.xbin
            ny     = int(param['Y1_SIZE']) // self.ybin
            self.win.append(_win(xstart,ystart,nx,ny))
            fsize += 2*self.win[-1].nx*self.win[-1].ny

            if self.mode == '2-USPEC':
                xstart = int(param['X2_START'])
                ystart = int(param['Y2_START'])
                nx     = int(param['X2_SIZE']) // self.xbin
                ny     = int(param['Y2_SIZE']) // self.ybin
                self.win.append(_win(xstart,ystart,nx,ny))
                fsize += 2*self.win[-1].nx*self.win[-1].ny

        if fsize != self.framesize:
            raise UltracamError('Rhead.__init__: file = ' + self.fname + '. Framesize = ' 
                                + str(self.framesize) + ' clashes with calculated value = ' +  str(fsize))

        # nasty stuff coming up ...
        self.revision  = int(param['REVISION']) if 'REVISION' in param else int(param['VERSION']) if 'VERSION' in param else -1
        self.urevision = int(user['revision']) if user is not None and 'revision' in user else -1
        self.version   = self.urevision if user is not None else self.revision

        if self.headerwords == 16:
            VERSIONS = [100222, 111205, 120716, 120813, 130307]
            if self.version not in VERSIONS:
                raise UltracamError('Rhead.__init__: could not recognise version = ' + str(self.version))

        elif user is not None:
            if self.revision != -1 and self.revision != self.urevision:
                raise UltracamError('Rhead.__init__: clashing revision numbers: ' + str(self.reversion) + ' vs ' + str(self.urevision))

        self.whichRun = ''
        if self.instrument == 'ULTRACAM':
            if user is None:
                self.time_units = 0.001
            else:
                self.time_units = 0.0001
            if 'REVISION' not in param and 'VERSION' not in param and 'V_FT_CLK' not in param:
                self.whichRun = 'MAY2002'
        else:
            if user is not None and self.headerwords == 16 and  self.revision >= 120813:
                self.time_units = 0.0001
            else:
                self.time_units = 0.001


    def isPonoff(self):
        """
        Is the run a power on / off? (no data)
        """
        return self.mode == 'PONOFF'

# Bit masks needed for Meinberg GPS data. 
# See description in read_header.cc in pipeline for more
PCPS_FREER            = 0x01   # DCF77 clock running on xtal, GPS receiver has not verified its position 
PCPS_DL_ENB           = 0x02   # daylight saving enabled 
PCPS_SYNCD            = 0x04   # clock has sync'ed at least once after pwr up 
PCPS_DL_ANN           = 0x08   # a change in daylight saving is announced 
PCPS_UTC              = 0x10   # a special UTC firmware is installed 
PCPS_LS_ANN           = 0x20   # leap second announced, (requires firmware rev. REV_PCPS_LS_ANN_...)
PCPS_IFTM             = 0x40   # the current time was set via PC, (requires firmware rev. REV_PCPS_IFTM_...) 
PCPS_INVT             = 0x80   # invalid time because battery was disconn'd
PCPS_LS_ENB           = 0x0100 # current second is leap second 
PCPS_ANT_FAIL         = 0x0200 # antenna failure 
PCPS_UCAP_OVERRUN     = 0x2000 # events interval too short 
PCPS_UCAP_BUFFER_FULL = 0x4000 # events read too slow 
PCPS_IO_BLOCKED       = 0x8000 # access to microprocessor blocked

class Rdata (Rhead):
    """
    Callable object to represent an Ultracam/spec raw data file. The idea is
    to open the file, and then the object generated can be used to deliver
    frames by specifying a frame number. Frames can be read individually e.g.

      rdat = Rdata('run045')
      fr10 = rdat(10)
      fr11 = rdat()

    reads frame numbers 10 and then 11 in from 'run045'. Frames can also read 
    sequentially. This is in fact a requirement if one wants accurate times when 
    it is often necessary to have read the timestamps of one or more immediately 
    preceding frames. One can do so as follows:

      rdat = Rdata('run045')
      while 1:
        try:
          frm = rdat()
          print 'nccd = ',frm.nccd()
        except:
          break

    or, since Rdata is defined as an iterator, more neatly as:

      for frm in Rdata('run045'):
         print 'nccd = ',dat.nccd()

    Rdata maintains an internal file object that is always at the start of a frame. 
    This enables sequential reads to be swift. If an attempt is made to access a 
    frame that does not exist, this is set to the start of the file (frame 1), and 
    this is the state at the end of the sequential read above for instance.
    """
    def __init__(self, run, nframe=1, flt=True):
        """
        Connects to a raw data file for reading. The file is kept open. 
        The file pointer is set to the start of frame nframe.

        Arguments:

        run     -- as in 'run036'. Will try to access equivalent .xml and .dat
                   files

        nframe  -- frame to position for next read, starting at 1 as the first.

        flt     -- True for reading data in as floats. This is the default for 
                   safety, however the data are stored on disk as unsigned 2-byte 
                   ints. If you are not doing much to the data, and wish to keep
                   them in this form for speed and efficiency, then set flt=False.
                   This parameter is used when iterating through an Rdata. The 
                   __call__ method can override it.
        """
        Rhead.__init__(self, run + '.xml')
        if self.isPonoff():
            raise UltracamError('Rdata.__init__: attempted to read a power on/off')
 
        # Attributes set are:
        #
        # _fobj   -- file object opened on data file
        # _nf     -- next frame to be read
        # _run    -- name of run
        # _flt    -- whether to read as float (else uint16)
        # _tstamp -- list of immediately preceding times
        self._fobj   = open(run + '.dat', 'rb')
        self._nf     = nframe
        self._run    = run
        self._flt    = flt
        self._tstamp = []
        if nframe != 1:
            self._fobj.seek(self.framesize*(nframe-1))
    
    def __iter__(self):
        """
        Generator to allow Rdata to function as an iterator
        """
        try:
            while 1:
                yield self.__call__(flt=self._flt)
        except UendError:
            pass

    def set(self, nframe=1):
        """
        Sets the internal file pointer to point at frame nframe.

        nframe  -- frame number to get, starting at 1. 0 for the last (complete) frame. 'None'
                   will be ignored. A value < 0 will cause an exception. A value greater than
                   the number of frames in the file will work, but will cause an exception to
                   be raised on the next attempted read.
        """

        # position read pointer
        if nframe is not None:
            if nframe < 0:
                raise UltracamError('Data.get: nframe < 0')
            elif nframe == 0:
                self._fobj.seek(0,2)
                fp = self._fobj.tell() 
                nf = fp // self.framesize
                self._fobj.seek(self.framesize*(nf-1)-fp,2)
                self._nf = nf
            elif self._nf != nframe:
                self._fobj.seek(self.framesize*(nframe-1))
                self._nf = nframe

    def __call__(self, nframe=None, flt=None):
        """
        Reads the data of frame nframe (starts from 1) and returns a
        corresponding CCD or MCCD object, depending upon the type of data. If
        nframe is None, just reads whatever frame we are on. Raises an
        exception if it fails to read data.  Resets to start of the file in
        this case. The data are stored internally as either 4-byte floats or
        2-byte unsigned ints.

        nframe -- frame number to get, starting at 1. 0 for the last
                  (complete) frame.

        flt    -- Set True to reading data in as floats. The data are stored on
                  disk as unsigned 2-byte ints. If you are not doing much to
                  the data, and wish to keep them in this form for speed and
                  efficiency, then set flt=False. If None then the value used when
                  consrtucting the MCCD will be used.
        """

        if flt is None: flt = self._flt

        # position read pointer
        self.set(nframe)

        # read timing bytes
        tbytes = self._fobj.read(2*self.headerwords)
        if len(self._tbytes) != 2*self.headerwords:
            self._fobj.seek(0)
            self._nf = 1
            raise UendError('Data.get: failed to read timing bytes')

        # debug
        try:
            print utimer(tbytes, self, self._tstamp)
        except Exception, err:
            print 'timing problem:'
            print err

        # read data
        buff = np.fromfile(self._fobj,'<u2',self.framesize/2-self.headerwords)
        if len(buff) != self.framesize/2-self.headerwords:
            self._fobj.seek(0)
            self._nf = 1
            raise UltracamError('Data.get: failed to read data')

        # move frame counter on by one
        self._nf += 1

        # build header
        head = Uhead()
        head.add_entry('User','Data entered by user at telescope')
        head.add_entry('Instrument','Instrument information')
        head.add_entry('Instrument.instrument',self.instrument,ITYPE_STRING,'Instrument identifier')
        head.add_entry('Instrument.headerwords',self.headerwords,ITYPE_INT,'Number of 2-byte words in timing')
        head.add_entry('Instrument.framesize',self.framesize,ITYPE_INT,'Total number of bytes per frame')
        head.add_entry('Data', 'File information')
        head.add_entry('Data.run',self._run,ITYPE_STRING,'run the frame came from')
        head.add_entry('Data.frame',self._nf,ITYPE_INT,'frame number within run')

        # interpret data
        xbin, ybin = self.xbin, self.ybin
        if self.instrument == 'ULTRACAM':
            # 3 CCDs. Windows come in pairs. Data from equivalent windows come out
            # on a pitch of 6. Some further jiggery-pokery is involved to get the
            # orientation of the frames correct.
            wins1, wins2, wins3 = [],[],[]
            noff = 0
            for wl, wr in zip(self.win[::2],self.win[1::2]):
                npix = 6*wl.nx*wl.ny
                if flt:
                    wins1.append(Window(np.reshape(buff[noff:noff+npix:6],(wl.ny,wl.nx)),
                                        wl.llx,wl.lly,xbin,ybin).astype(np.float32))
                    wins1.append(Window(np.reshape(buff[noff+npix-5:noff:-6],(wr.ny,wr.nx)),
                                        llxr,llyr,xbin,ybin).astype(np.float32))
                    wins2.append(Window(np.reshape(buff[noff+2:noff+npix:6],(wl.ny,wl.nx)),
                                        wl.llx,wl.lly,xbin,ybin).astype(np.float32))
                    wins2.append(Window(np.reshape(buff[noff+npix-3:noff:-6],(wr.ny,wr.nx)),
                                        llxr,llyr,xbin,ybin).astype(np.float32))
                    wins3.append(Window(np.reshape(buff[noff+4:noff+npix:6],(wl.ny,wl.nx)),
                                        wl.llx,wl.lly,xbin,ybin).astype(np.float32))
                    wins3.append(Window(np.reshape(buff[noff+npix-1:noff:-6],(wr.ny,wr.nx)),
                                        llxr,llyr,xbin,ybin).astype(np.float32))
                else:
                    wins1.append(Window(np.reshape(buff[noff:noff+npix:6],(wl.ny,wl.nx)),
                                        wl.llx,wl.lly,xbin,ybin))
                    wins1.append(Window(np.reshape(buff[noff+npix-5:noff:-6],(wr.ny,wr.nx)),
                                        llxr,llyr,xbin,ybin))
                    wins2.append(Window(np.reshape(buff[noff+2:noff+npix:6],(wl.ny,wl.nx)),
                                        wl.llx,wl.lly,xbin,ybin))
                    wins2.append(Window(np.reshape(buff[noff+npix-3:noff:-6],(wr.ny,wr.nx)),
                                        llxr,llyr,xbin,ybin))
                    wins3.append(Window(np.reshape(buff[noff+4:noff+npix:6],(wl.ny,wl.nx)),
                                        wl.llx,wl.lly,xbin,ybin))
                    wins3.append(Window(np.reshape(buff[noff+npix-1:noff:-6],(wr.ny,wr.nx)),
                                        llxr,llyr,xbin,ybin))
                noff += npix

            # Build CCDs
            ccd1 = CCD(wins1, None, self.nxmax, self.nymax, True, None)
            ccd2 = CCD(wins2, None, self.nxmax, self.nymax, True, None)
            ccd3 = CCD(wins3, None, self.nxmax, self.nymax, True, None)

            # Return an MCCD
            return MCCD([ccd1,ccd2,ccd3], head)
        else:
            raise UltracamError('Have yet to implement anything for ' + self.instrument)

    def ntotal(self):
        """
        Returns total number of frames in data file
        """
        self._fobj.seek(0,2)
        ntot = self._fobj.tell() // self.framesize
        self._fobj.seek(self.framesize*(self._nf-1))
        return ntot

    def nframe(self):
        """
        Returns next frame number to be read if reading
        sequentially (starts at 1)
        """
        return self._nf

class Rtime (Rhead):
    """
    Iterator class to enable swift reading of Ultracam times.
    """
    def __init__(self, run, nframe=1):
        """
        Connects to a raw data file for reading. The file is kept open. 
        The file pointer is set to the start of frame nframe.

        Arguments:

        run     -- as in 'run036'. Will try to access equivalent .xml and .dat
                   files

        nframe  -- frame to position for next read, starting at 1 as the first.
        """
        Rhead.__init__(self, run + '.xml')
        if self.isPonoff():
            raise UltracamError('Rtime.__init__: attempted to read a power on/off')
 
        # Attributes set are:
        #
        # _fobj   -- file object opened on data file
        # _nf     -- next frame to be read
        # _run    -- name of run
        self._fobj   = open(run + '.dat', 'rb')
        self._nf     = nframe
        self._run    = run
        if nframe != 1:
            self._fobj.seek(self.framesize*(nframe-1))
    
    def __iter__(self):
        """
        Generator to allow Rdata to function as an iterator
        """
        try:
            while 1:
                yield self.__call__()
        except UendError:
            pass

    def set(self, nframe=1):
        """
        Sets the internal file pointer to point at frame nframe.

        nframe  -- frame number to get, starting at 1. 0 for the last (complete) frame. 'None'
                   will be ignored. A value < 0 will cause an exception. A value greater than
                   the number of frames in the file will work, but will cause an exception to
                   be raised on the next attempted read.
        """

        # position read pointer
        if nframe is not None:
            if nframe < 0:
                raise UltracamError('Data.get: nframe < 0')
            elif nframe == 0:
                self._fobj.seek(0,2)
                fp = self._fobj.tell() 
                nf = fp // self.framesize
                self._fobj.seek(self.framesize*(nf-1)-fp,2)
                self._nf = nf
            elif self._nf != nframe:
                self._fobj.seek(self.framesize*(nframe-1))
                self._nf = nframe

    def __call__(self, nframe=None):
        """
        Returns timing information of frame nframe (starts from 1).

        nframe -- frame number to get, starting at 1. 0 for the last
                  (complete) frame.

        See utimer for gets returned by this.
        """

        # position read pointer
        self.set(nframe)

        # read timing bytes
        tbytes = self._fobj.read(2*self.headerwords)
        if len(tbytes) != 2*self.headerwords:
            self._fobj.seek(0)
            self._nf = 1
            raise UendError('Data.get: failed to read timing bytes')

        tinfo = utimer(tbytes, self, self._nf)
        
        # step to start of next frame
        self._fobj.seek(self.framesize-2*self.headerwords,1)

        # move frame counter on by one
        self._nf += 1

        return tinfo

# Some fixed dates needed by utimer. Put them here so
# that they are only computed once.
DSEC           = 86400
MJD0           = datetime.date(1858,11,17).toordinal()
UNIX           = datetime.date(1970,1,1).toordinal()
DEFDAT         = datetime.date(2000,1,1).toordinal()
MAY2002        = datetime.date(2001,5,12).toordinal()
SEP2002        = datetime.date(2002,9,8).toordinal()
TSTAMP_CHANGE1 = datetime.date(2003,8,1).toordinal()
TSTAMP_CHANGE2 = datetime.date(2005,1,1).toordinal()
TSTAMP_CHANGE3 = datetime.date(2010,3,1).toordinal()
USPEC_CHANGE   = datetime.date(2011,9,21).toordinal()

# ULTRACAM Timing parameters from Vik
INVERSION_DELAY = 110.          # microseconds
HCLOCK          = 0.48          # microseconds
CDS_TIME_FDD    = 2.2           # microseconds
CDS_TIME_FBB    = 4.4           # microseconds
CDS_TIME_CDD    = 10.           # microseconds
SWITCH_TIME     = 1.2           # microseconds

USPEC_FT_ROW    = 14.4e-6       # seconds
USPEC_FT_OFF    = 49.e-6        # seconds
USPEC_CLR_TIME  = 0.0309516     # seconds 

def utimer(tbytes, rhead, fnum):
    """
    Computes the Time corresponding of the most recently read frame, 
    None if no frame has been read. For the Time to be reliable 
    (Time.good == True), several frames might have needed to
    be read and their times passed through tstamp.

     tbytes  -- string of timing bytes

     rhead   -- Rhead header of data file

     fnum    -- frame number we think we are on.



    Returns (time,blueTime,badBlue,gps,nsat,format,vclock_frame,whichRun)

     time         -- the Time as best as can be determined

     blueTime     -- different time for the blue frames of ULTRACAM for nblue > 1

     badBlue      -- blue frame is bad for nblue > 1 (nblue-1 out nblue are bad)

     gps          -- the raw GPS time associated with the frame, no corrections applied
                    other than for various bugs in the timing. i.e. not corrected to 
                    mid-exposure, or for drift mode etc.

     nsat         -- number of satellites (only for early configurations of GPS card, otherwise
                    comes back as -1)

     format       -- the format integer used to translate the timing bytes

     vclock_frame -- vertical clocking time for a whole frame

     whichRun     -- run identifier
    """
    
    # This is an involved routine owing to the various hardware
    # changes and bugs that have cropped up over the years. The 
    # pipeline equivalent routine is read_header.cc
    
    # return immediately if no bytes have been read.
    if tbytes is None: return None

    # start by assuming the time will be good, but many things
    # can wreck this. Once False, it can never return to True
    goodTime = True
    reason   = ''

    if rhead.instrument == 'ULTRASPEC' and rhead.version == -1:
        format = 2
    elif rhead.version == -1 or rhead.version == 70514 or rhead.version == 80127:
        format = 1
    elif rhead.version == 100222 or rhead.version == 110921 or rhead.version == 111205 or \
            rhead.version == 120716 or rhead.version == 120813:
        format = 2
    else:
        raise UltracamError('Rdata._timing: version = ' + str(rhead.version) + ' unrecognised.')
        
    frameNumber = struct.unpack('<I', tbytes[4:8])[0]
    if frameNumber != fnum:
        warnings.warn('ultracam.utimer: run ' + rhead._run + ' expected frame number = ' + 
                      str(fnum) + ' but found ' + str(frameNumber))

    # initialize some attributes of utimer which are used as the equivalent
    # of static variables in C++. They are:
    #
    #  previousFrameNumber  -- frame number of immediately preceding frame read, 0 if none
    #  tstamp               -- list of raw GPS times of preceding frames, [0] most recent
    #  blueTimes            -- list of modified times of preceding frames, [0] most recent

    if hasattr(utimer,'previousFrameNumber'):
        if frameNumber != utimer.previousFrameNumber + 1:
            utimer.tstamp    = []
            utimer.blueTimes = []
    else:
        utimer.tstamp    = []
        utimer.blueTimes = []

    utimer.previousFrameNumber = frameNumber


    if format == 1:
        nsec, nnsec = struct.unpack('<Ii', tbytes[9:17])
        nsat = struct.unpack('<h', tbytes[21:23])[0]
        if nsat <= 2:
            goodTime = False
            reason   = 'too few satellites (' + str(nsat) + ')'

    elif format == 2:

        nexp  = struct.unpack('<I', tbytes[8:12])
        if nexp*rhead.tunits != rhead.expose:
            goodTime = False
            reason = 'XML expose time does not match time in timing bytes.'
        nsec, nnsec = struct.unpack('<Ii', tbytes[12:16])
        nnsec *= 100
        nsat   = None
        tstamp = struct.unpack('<H', tbytes[24:26])[0]
        
        if goodTime and (tstamp & PCPS_ANT_FAIL):
            goodTime = False
            reason   = 'GPS antenna failed'

        if goodTime and (tstamp & PCPS_INVT):
            goodTime = False
            reason   = 'GPS battery disconnected'

        if goodTime and not (tstamp & PCPS_SYNCD):
            goodTime = False
            reason = 'GPS clock not yet synced since power up'

        if goodTime and (tstamp & PCPS_FREER):
            goodTime = False
            reason = 'GPS receiver has not verified its position'

    else:
        raise UltracamError('Rdata.time: format = ' + str(format) + ' not recognised.')

    # One of the bits in the first byte is set if the blue frame is junk. 
    # Unfortunately which bit was set changed hence the check of the format
    fbyte   = struct.unpack('<B',tbytes[0])[0]
    badBlue = rhead.nblue > 1 and \
        ((format == 1 and bool(fbyte & 1<<3)) or (format == 2 and bool(fbyte & 1<<4)))

    def tcon1(offset, nsec, nnsec):
        """
        Convert to MJD
        """
        return offset + float(nsec+nnsec/1.e9)/DSEC - MJD0

    def tcon2(year, month, day, nsec, nnsec):
        """
        Convert to MJD
        """
        return datetime.date(year,month,day).toordinal() + float(nsec+nnsec/1.e9)/DSEC - MJD0

    if format == 1 and nsat == -1:
        goodTime = False
        reason = 'no satellites.'
        mjd = tcon1(DEFDAT, nsec, nnsec)
        if rhead.v_ft_clk > 127:
            vclock_frame = 6.e-9*(40+320*(rhead.v_ft_clk - 128))
	else:
	    vclock_frame = 6.e-9*(40+40*rhead.v_ft_clk)

    else:
        
        if rhead.whichRun == 'MAY2002' and format == 1:
            # Had no date info in the timestamps of this run
            # nsec was offset from 12 May 2002
            mjd = tcon1(MAY2002, nsec, nnsec)

            # Some of the nights ran over the end of the week. Spot
            # from dates prior to the start of the run on 16 May.
            if mjd < MAY2002+4: mjd += 7

            # Correct 10 second error that affected the May 2002 run.
            # Only possible if we have read the previous frames, with
            # time stored in an attribute called tstamp of this function
            if len(utimer.tstamp) and mjd < utimer.tstamp[0]: mjd += 10./DSEC

            # Fix problem with very first night
            if mjd < MAY2002+5.5:
                vclock_frame = 10.0e-6
            else:
                vclock_frame = 24.46e-6

        else:

            # OK now have date info, although we
            # need a stack of special case fixes
            if format == 1:
                day, month, year = struct.unpack('<BBH',tbytes[17:21])
                
                if month == 9 and year == 263: year = 2002
                if year < 2002:
                    mjd = tcon1(SEP2002, nsec, nnsec)

                elif month == 9 and year == 2002:
                    tdiff = datetime.date(year,month,day).toordinal()-SEP2002
                    nweek = tdiff // 7
                    days  = tdiff - 7*nweek
                    if days > 3 and nsec < 2*DSEC:
                        nweek += 1
                    elif days <= 3 and nsec > 5*DSEC:
                        nweek -= 1
                    mjd = tcon1(SEP2002+7*nweek,nsec,nnsec)
                else:
                    mjd = tcon2(year, month, day, nsec % DSEC, nnsec)

            elif format == 2:
                mjd = tcon1(UNIX, nsec, nnsec)

            if mjd > TSTAMP_CHANGE1:
                if rhead.v_ft_clk > 127:
                    vclock_frame = 6.e-9*(40+320*(rhead.v_ft_clk - 128))
                else:
                    vclock_frame = 6.e-9*(40+40*rhead.v_ft_clk)

            else:
                if rhead.v_ft_clk > 127:
                    vclock_frame = 6.e-9*(80+160*(rhead.v_ft_clk - 128))
                else:
                    vclock_frame = 6.e-9*(80+20*rhead.v_ft_clk)

    # 'midnight bug' correction
    if (int(mjd-3) % 7) == ((nsec // DSEC) % 7):
        warnings.warn('ultracam.utimer: run ' + rhead._run + ' midnight bug detected and corrected')
        mjd += 1

    # save this as the raw GPS time.
    gps = mjd

    # next variable determines when the timestamp is assumed to be taken within the 
    # read cycle which has changed a few times owing to various small mishaps.
    defTstamp = mjd < TSTAMP_CHANGE1 or (mjd > TSTAMP_CHANGE2 and mjd < TSTAMP_CHANGE3) 

    # Push time to front of tstamp list
    utimer.tstamp.insert(0,mjd)

    # one extra parameter in addition to those from Vik's
    VCLOCK_STORAGE = vclock_frame
    USPEC_FT_TIME  = 0.0067196 if mjd < USPEC_CHANGE else 0.0149818
    
    if rhead.instrument == 'ULTRACAM':
        if rhead.gainSpeed == 'cdd':
            cds_time = CDS_TIME_CDD
        elif rhead.gainSpeed == 'fbb':
            cds_time = CDS_TIME_FBB
        elif rhead.gainSpeed == 'fdd':
            cds_time = CDS_TIME_FDD
        else:
            raise UltracamError('ultracam.utimer: did not recognize gain speed setting = ' + rhead.gainSpeed)
    elif rhead.instrument == 'ULTRASPEC':
        cdsTime = 0.
        
    VIDEO = SWITCH_TIME + cds_time

    if rhead.instrument == 'ULTRACAM' and \
            (rhead.mode == 'FFCLR' or rhead.mode == 'FFOVER' or rhead.mode == '1-PCLR'):

        # never need more than two times
        if len(utimer.tstamp) > 2: utimer.tstamp.pop()

        if defTstamp:
            mjdCentre  = utimer.tstamp[0]
            mjdCentre += rhead.exposeTime/DSEC/2.
            exposure   = rhead.exposeTime

        else:

            # Time taken to clear CCD
            clear_time   = (1033. + 1027)*vclock_frame

            # Time taken to read CCD (assuming cdd mode)
            if rhead.mode == 'FFCLR':
                readout_time = (1024/rhead.ybin)*(VCLOCK_STORAGE*rhead.ybin + 
                                                  536*HCLOCK + (512/rhead.xbin+2)*VIDEO)/1.e6
            elif rhead.mode == 'FFOVER':
                readout_time = (1032/rhead.ybin)*(VCLOCK_STORAGE*rhead.ybin + 
                                                  540*HCLOCK + (540/rhead.xbin+2)*VIDEO)/1.e6
            else:
                nxb          = rhead.win[1].nx
                nxu          = rhead.xbin*nxb
                nyb          = rhead.win[1].ny
                xleft        = rhead.win[0].llx
                xright       = rhead.win[1].llx + nxu - 1
                diff_shift   = abs(xleft - 1 - (1024 - xright) )
                num_hclocks  =  nxu + diff_shift + (1024 - xright) + 8 \
                    if (xleft - 1 > 1024 - xright) else nxu + diff_shift + (xleft - 1) + 8
                readout_time = nyb*(VCLOCK_STORAGE*serverdata.ybin + 
                                    num_hclocks*HCLOCK + (nxb+2)*VIDEO)/1.e6

            # Frame transfer time
            frameTransfer = 1033.*vclock_frame

            if len(utimer.tstamp) == 1:

                # Case where we have not got a previous timestamp. Hop back over the 
                # readout and frame transfer and half the exposure delay
                mjdCentre  = utimer.tstamp[0]
                mjdCentre -= (frameTransfer+readout_time+serverdata.expose_time/2.)/DSEC
                if goodTime:
                    goodTime = False
                    reason = 'no previous GPS time found in non-default mode'

            else:

                # Case where we have got previous timestamp is somewhat easier and perhaps
                # more reliable since we merely need to step forward over the clear time and
                # half the exposure time.
                mjdCentre  = utimer.tstamp[1]
                mjdCentre += (clear_time + rhead.exposeTime/2.)/DSEC

            exposure = rhead.exposeTime

    elif rhead.instrument == 'ULTRACAM' and \
            (rhead.mode == 'FFNCLR' or rhead.mode == '1-PAIR' 
             or rhead.mode == '2-PAIR' or rhead.mode == '3-PAIR'):

        # never need more than three times
        if len(utimer.tstamp) > 3: utimer.tstamp.pop()

        # Time taken to move 1033 rows.
        frameTransfer = 1033.*vclock_frame

        if rhead.mode == 'FFNCLR':
            readout_time = (1024/rhead.ybin)*(VCLOCK_STORAGE*rhead.ybin + 
                                              536*HCLOCK + (512/rhead.xbin+2)*VIDEO)/1.e6
        else:

            readout_time = 0.
            xbin = rhead.xbin
            ybin = rhead.ybin
            ystart_old = -1
            for wl, wr in zip(rhead.win[::2],rhead.win[1::2]):

                nxu  = xbin*wl.nx
                nyu  = ybin*wl.ny
			
                ystart = wl.lly
                xleft  = wl.llx
                xright = wr.llx + nxu - 1
	  
                if ystart_old > -1:
                    ystart_m = ystart_old
                    nyu_m    = nyu_old
                    y_shift  = (ystart-ystart_m-nyu_m)*VCLOCK_STORAGE
                else:
                    ystart_m = 1
                    nyu_m    = 0
                    y_shift  = (ystart-1)*VCLOCK_STORAGE

                # store for next time
                ystart_old = ystart
                nyu_old    = nyu
			
                # Number of columns to shift whichever window is further 
                # from the edge of the readout to get ready for simultaneous 
                # readout.
                diff_shift = abs(xleft - 1 - (1024 - xright) )

                # Time taken to dump any pixels in a row that come after the ones we want.
                # The '8' is the number of HCLOCKs needed to open the serial register dump gates
                # If the left window is further from the left edge than the right window is from the
                # right edge, then the diffshift will move it to be the same as the right window, and
                # so we use the right window parameters to determine the number of hclocks needed, and
                # vice versa.
                num_hclocks = nxu + diff_shift + (1024 - xright) + 8 \
                    if (xleft - 1 > 1024 - xright) else nxu + diff_shift + (xleft - 1) + 8
			
                # Time taken to read one line. The extra 2 is required to fill the video pipeline buffer
                line_read = VCLOCK_STORAGE*ybin + num_hclocks*HCLOCK + (nxu/xbin+2)*VIDEO
		
                readout_time += y_shift + (nyu/ybin)*line_read

            readout_time /= 1.e6
	          
	if defTstamp:
            if frameNumber == 1:
                mjdCentre  = utimer.tstamp[0]
                exposure   = rhead.exposeTime
                mjdCentre -= (frameTransfer+exposure/2.)/DSEC

            else:
                if len(utimer.tstamp) > 1:
                    texp = DSEC*(utimer.tstamp[0] - utimer.tstamp[1]) - frameTransfer
                    mjdCentre  = utimer.tstamp[1]
                    mjdCentre += texp/DSEC
                    exposure   = texp

                else:
                    texp = readout_time + rhead.exposeTime
                    mjdCentre  = utimer.tstamp[0]
                    mjdCentre -= (frameTransfer+texp/2.)/DSEC
                    exposure   = texp

                    if goodTime:
                        goodTime = False
                        reason = 'could not establish an accurate time without previous GPS timestamp'

        else:
            if frameNumber == 1:
                mjdCentre  = utimer.tstamp[0]
                exposure   = rhead.exposeTime
                mjdCentre -= (frameTransfer+exposure/2.)/DSEC

                if goodTime:
                    goodTime = False
                    reason = 'cannot establish an accurate time for first frame in this mode'

            else:
                
                if len(utimer.tstamp) > 2:
                    texp       = DSEC*(utimer.tstamp[1] - utimer.tstamp[2]) - frameTransfer 
                    mjdCentre  = utimer.tstamp[1]
                    mjdCentre += (rhead.exposeTime - texp/2.)/DSEC
                    exposure   = texp

                elif len(utimer.tstamp) == 2:
                    texp = DSEC*(utimer.tstamp[0] - utimer.tstamp[1]) - frameTransfer 
                    mjdCentre  = utimer.tstamp[1]
                    mjdCentre += (rhead.exposeTime - texp/2.)/DSEC
                    exposure   = texp

                    if goodTime:
                        goodTime = False
                        reason = 'cannot establish an accurate time with only two prior timestamps'

                else:
                    texp       = readout_time + rhead.exposeTime
                    mjdCentre  = utimer.tstamp[0]
                    mjdCentre += (rhead.exposeTime-texp-frameTransfer-texp/2.)/DSEC
                    exposure   = texp

                    if goodTime:
                        goodTime = False
                        reason = 'cannot establish an accurate time with only one prior timestamp'

    elif rhead.instrument == 'ULTRACAM' and rhead.mode == 'DRIFT':

        wl     = rhead.win[0]
        xbin   = rhead.xbin
        ybin   = rhead.xbin
        nxu    = xbin*wl.nx
        nyu    = ybin*wl.ny
        ystart = wl.lly
        xleft  = wl.llx
        wr     = rhead.win[1]
        xright = wl.llx + nxu -1

        # Maximum number of windows in pipeline
        nwins = int((1033./nyu+1.)/2.)
        pipe_shift = int(1033.-(((2.*nwins)-1.)*nyu))

        # Time taken for (reduced) frame transfer, the main advantage of drift mode
        frameTransfer = (nyu + ystart - 1)*vclock_frame
	  
        # Number of columns to shift whichever window is further from the edge of the readout
        # to get ready for simultaneous readout.
        diff_shift = abs(xleft - 1 - (1024 - xright) )

        # Time taken to dump any pixels in a row that come after the ones we want.
        # The '8' is the number of HCLOCKs needed to open the serial register dump gates
        # If the left window is further from the left edge than the right window is from the
        # right edge, then the diffshift will move it to be the same as the right window, and
        # so we use the right window parameters to determine the number of hclocks needed, and
        # vice versa.
        num_hclocks  = nxu + diff_shift + (1024 - xright) + 8 \
            if (xleft - 1 > 1024 - xright) else nxu + diff_shift + (xleft - 1) + 8
			
        # Time taken to read one line. The extra 2 is required to fill the video pipeline buffer
        line_read = VCLOCK_STORAGE*ybin + num_hclocks*HCLOCK + (nxu/xbin+2)*VIDEO
		
        readout_time = ((nyu/ybin)*line_read + pipe_shift*VCLOCK_STORAGE)/1.e6

	# Never need more than nwins+2 times
	if len(utimer.tstamp) > nwins+2: utimer.tstamp.pop()

	if defTstamp:

	    # Pre board change or post-bug fix
	    if len(utimer.tstamp) > nwins:

		texp = utimer.tstamp[nwins-1] - utimer.tstamp[nwins] - frameTransfer
		mjdCentre  = utimer.tstamp[nwins]
		mjdCentre += texp/2./DSEC
		exposure   = texp

	    else:

		# Set to silly value for easy checking
		mjdCentre = DEFDAT
		exposure  = read.exposeTime

		if goodTime:
                    goodTime = False
		    reason = 'too few stored timestamps for drift mode'

	else:

	    if len(utimer.tstamp) > nwins+1:

		texp = utimer.tstamp[nwins] - utimer.tstamp[nwins+1] - frameTransfer
		mjdCentre  = utimer.tstamp[nwins]
		mjdCentre += (rhead.expose_time-texp/2.)/DSEC
		exposure   = texp
                
            elif len(utimer.tstamp) == nwins+1:

		texp       = utimer.tstamp[nwins-1] - utimer.tstamp[nwins] - frameTransfer
		mjdCentre  = utimer.tstamp[nwins]
		mjdCentre += (serverdata.expose_time-texp/2.)/DSEC
		exposure   = texp

		if goodTime:
                    goodTime = False
		    reason = 'too few stored timestamps for drift mode'

            else:
	  
		# Set to silly value for easy checking
		mjdCentre = DEFDAT
		exposure  = rhead.exposeTime

		if goodTime:
                    goodTime = False
		    reason = 'too few stored timestamps for drift mode'

    elif rhead.instrument == 'ULTRASPEC' and \
            (rhead.mode == '1-USPEC' or rhead.mode == '1-USPEC'):
 
	# Avoid excessive accumulation of timestamps.
	if len(utimer.tstamp) > 3: utimer.tstamp.pop()

        if utimer.tstamp[0] < USPEC_CHANGE:  
            
            texp = readout_time + serverdata.expose_time
            mjdCentre = utimer.tstamp[0]
            if rhead.en_clr or frameNumber == 1:
                
                mjdCentre -= rhead.exposeTime/2./DSEC
                exposure   = rhead.exposeTime
                
            elif len(utimer.tstamp) > 1:
                
                texp = utimer.tstamp[0] - utimer.tstamp[1] - USPEC_FT_TIME
                mjdCentre -= text/2./DSEC
                exposure   = texp
                
            else:
                
                # Could be improved with an estimate of the read time
                mjdCentre -= rhead.exposeTime/2./DSEC
                exposure   = rhead.exposeTime

                if goodTime:
                    reason   = 'too few stored timestamps'
                    goodTime = False

        else:

            if rhead.en_clr or frameNumber == 1:
                # Special case for the first frame or if clears are enabled.
                exposure = rhead.exposeTime
                if len(utimer.tstamp) == 1:
                    mjdCentre = utimer.tstamp[0]
                    mjdCentre -= (-USPEC_FT_TIME-serverdata.expose_time/2.)/DSEC
                    if goodTime:
                        reason = 'cannot establish an accurate time without at least 1 prior timestamp'
                        goodTime = False
                else:
                    mjdCentre  = utimer.tstamp[1]
                    mjdCentre += (USPEC_CLR_TIME+serverdata.expose_time/2.)/DSEC

            elif len(utimer.tstamp) > 2:

                # Can backtrack two frames to get a good exposure time.
                texp = utimer.tstamp[1] - utimer.tstamp[2] - USPEC_FT_TIME
                mjdCentre  = utimer.tstamp[1]
                mjdCentre += (serverdata.expose_time-texp/2.)/DSEC
                exposure   = texp

            elif len(utimer.tstamp) == 2:

                # Can only back up one, so estimate of exposure time is
                # actually based on the exposure following the one of
                # interest. Probably not too bad, but technically unreliable
                # as a time.
                texp = utimer.tstamp[0] - utimer.tstamp[1] - USPEC_FT_TIME
                mjdCentre  = utimer.tstamp[1]
                mjdCentre += (serverdata.expose_time-texp/2.)/DSEC
                exposure   = texp

                if goodTime:
                    reason = 'cannot establish an accurate time without at least 2 prior timestamps'
                    goodTime = False
	  
            else:

                # Only one time
                mjdCentre  = utimer.tstamp[0]
                mjdCentre -= (serverdata.expose_time/2.+serverdata.expose_time)/DSEC
                exposure   = serverdata.exposeTime

                if goodTime:
                    reason = 'too few stored timestamps'
                    goodTime = False

    elif rhead.instrument == 'ULTRASPEC' or rhead.mode == 'UDRIFT':

        ybin   = rhead.ybin
        nyu    = ybin*rhead.win[0].ny
        ystart = rhead.win[0].lly;
        nwins  = int(((1037. / nyu) + 1.)/2.)
        frameTransfer = USPEC_FT_ROW*(ystart+ny-1.)+USPEC_FT_OFF

	# Never need more than nwins+2 times
	if len(utimer.tstamp) > nwins+2: utimer.tstamp.pop() 

	# Non-standard mode

	if len(utimer.tstamp) > nwins+1:

	    texp       = utimer.tstamp[nwins] - utimer.tstamp[nwins+1] - frameTransfer
	    mjdCentre  = utimer.tstamp[nwins]
	    mjdCentre += (serverdata.expose_time-texp/2.)/DSEC
	    exposure   = texp
	    
	elif len(utimer.tstamp) == nwins+1:
	    
	    texp          = utimer.tstamp[nwins-1] - utimer.tstamp[nwins] - frameTransfer
	    mjdCentre     = utimer.tstamp[nwins]
	    mjdCentre     = (serverdata.expose_time-texp/2.)/DSEC
	    exposure      = texp
	    if goodTime:
		reason = 'too few stored timestamps'
		goodTime = False
	    
	else:
	  
	    # Set to silly value for easy checking
	    mjdCentre  = DEFDAT
	    exposure   = rhead.exposeTime
	    if goodTime:
		reason = 'too few stored timestamps'
		goodTime = False
  
    # Return some data
    time = Time(mjdCentre, exposure, goodTime, reason)

    if rhead.nblue > 1:

	# The mid-exposure time for the OK blue frames in this case is computed by averaging the 
	# mid-exposure times of all the contributing frames, if they are available.
	utimer.blueTimes.insert(0,time)

	if badBlue:

	    # just pass through the standard time for the junk frames
            blueTime = time

        else:

	    # if any of the contributing times is flagged as unreliable, then
	    # so is the final time. This is also unreliable if any
	    # contributing frame times are missing. Time is calculated as
	    # half-way point between start of first and end of last
	    # contributing exposure.  Corrections are made if there are too
	    # few contributing exposures (even though the final value will
	    # still be flagged as unreliable
	    ncont  = min(rhead.nblue, len(utimer.blueTimes))
	    start  = utimer.blueTimes[ncont-1].mjd - utimer.blueTimes[ncont-1].expose/2./DSEC
	    end    = utimer.blueTimes[0].mjd       + utimer.blueTimes[0].expose/2./DSEC
	    expose = end - start

	    # correct the times
	    ok = ncont == rhead.nblue
            reason = ''
	    if not ok:
		expose *= rhead.nblue/float(ncont)
		start   = end - expose
                reason  = 'not all contributing frames found'
	    else:
		ok = utimer.blueTimes[0].good and utimer.blueTimes[ncont-1].good
                if not ok: reason  = 'time of start or end frame was unreliable'

            blueTime = Time((start+end)/2., time, ok, reason)

        # Avoid wasting memory storing past times
	if len(utimer.blueTimes) > rhead.nblue: utimer.blueTimes.pop()
	    
    else:
        blueTime = time

    return (time,blueTime,badBlue,gps,nsat,format,vclock_frame,rhead.whichRun)

# Exception classes
class UltracamError(Exception):
    """For throwing exceptions from the ultracam module"""
    def __init__(self, value):
        self.value = value
            
    def __str__(self):
        return repr(self.value)

class UendError(UltracamError):
    """
    Exception for the standard way to reach the end of a data 
    file (failure to read the timing bytes). This allows the 
    iterator to die silently in this case while  raising
    exceptions for less anticipated cases.
    """
    def __init__(self, value):
        self.value = value
            
    def __str__(self):
        return repr(self.value)
