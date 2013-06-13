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

# less widely known extras
import ppgplot as pg

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

class Odict(dict):
    """
    A dictionary which stores a key order which it uses for printing
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

    mjd    -- modified Julian day, decimal.
    expose -- exposure time, seconds.
    good   -- is the time thought to be reliable?
    """
    def __init__(self, mjd, expose, good):
        self.mjd    = mjd
        self.expose = expose
        self.good   = good

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
            raise UltracamError('Uhead.__setitem__ takes either 2 or 4 arguments')

        if not isinstance(key, basestring):
            raise UltracamError('Uhead.__setitem__: argument "key" must be a string.')

        if not isinstance(comment, basestring):
            raise UltracamError('Uhead.__setitem__: key = ' + key + ': "comment" must be a string.')

        # now look at the key: must have no blanks
        if key.find(' ') > -1:
            raise UltracamError('Uhead.__setitem__: key = "' + key + '" contains at least one blank.')

        # if key has a '.' then the part before last dot must already exist
        # and must be a directory.
        ldot = key.rfind('.')
        if ldot > -1:
            dir = key[:ldot]
            for kold in self.keys():
                if dir == kold and self[kold][1] == ITYPE_DIR:
                    break
            else:
                raise UltracamError('Uhead.__setitem__: key = ' + key + 
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
                raise UltracamError('Uhead.__setitem__: key = ' + key + 
                                ': require a 2-element tuple or list (int,float) for ITYPE_TIME)')
            value[0] = int(value[0])
            value[1] = float(value[1])
        elif itype == ITYPE_DVECTOR:
            if not isinstance(value, np.ndarray) or len(value.shape) != 1:
                raise UltracamError('Uhead.__setitem__: key = ' + key + 
                                ': require a 1D numpy.ndarray for ITYPE_DVECTOR)')
            value = value.astype(float64)
        elif itype == ITYPE_IVECTOR:
            if not isinstance(value, np.ndarray) or len(value.shape) != 1:
                raise UltracamError('Uhead.__setitem__: key = ' + key + 
                                ': require a 1D numpy.ndarray for ITYPE_IVECTOR)')
            value = value.astype(int)
        elif itype == ITYPE_FVECTOR:
            if not isinstance(value, np.ndarray) or len(value.shape) != 1:
                raise UltracamError('Uhead.__setitem__: key = ' + key + 
                                ': require a 1D numpy.ndarray for ITYPE_FVECTOR)')
            value = value.astype(float32)
        else:
            raise UltracamError('Uhead.__setitem__: key = ' + key + 
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
            if issubclass(win.dtype.type,np.float):
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
                raise Warning('CCD.toInt: input data out of range 0 to 65535')
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

    # arithematic
    def __add__(self, other):
        """
        Adds 'other' to the CCD in place (+=)
        """
        for win,owin in zip(self,other):
            win += owin
        return self

    def __sub__(self, other):
        """
        Subtracts 'other' from the CCD in place (-=)
        """
        for win,owin in zip(self,other):
            win -= owin
        return self

    def __mul__(self, other):
        """
        Multiplies the CCD by 'other' in place (*=)
        """
        for win,owin in zip(self,other):
            win *= owin
        return self

    def __div__(self, other):
        """
        Divides the CCD by 'other' in place (/=)
        """
        for win,owin in zip(self,other):
            win /= owin
        return self

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

    def __eq__(self, other):
        """
        Equality operator tests same number of CCDs and that each CCD matches.
        """

        if type(other) is type(self):

            if self.nccd() != other.nccd(): return False

            for sccd, occd in zip(self,other):
                if len(ccd1) != len(ccd2): return False
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

            write_string(uf, key)

            value, itype, comment = val
            uf.write(struct.pack('i',itype))
            write_string(uf, comment)

            if itype == ITYPE_DOUBLE:
                uf.write(struct.pack('d', value))
            elif itype == ITYPE_INT:
                uf.write(struct.pack('i', value))
            elif itype == ITYPE_UINT:
                uf.write(struct.pack('I', value))
            elif itype == ITYPE_FLOAT:
                uf.write(struct.pack('f', value))
            elif itype == ITYPE_STRING:
                write_string(uf, value)
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
                    win.astype(np.uint16).tofile(uf)

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
    def __add__(self, other):
        """
        Adds 'other' to the MCCD in place (+=)
        """
        for ccd,occd in zip(self,other):
            ccd += occd
        return self

    def __sub__(self, other):
        """
        Subtracts 'other' from the MCCD in place (-=)
        """
        for ccd,occd in zip(self,other):
            ccd -= occd
        return self

    def __mul__(self, other):
        """
        Multiplies the MCCD by 'other' in place (*=)
        """
        for ccd,occd in zip(self,other):
            ccd *= occd
        return self

    def __div__(self, other):
        """
        Divides the MCCD by 'other' in place (/=)
        """
        for ccd,occd in zip(self,other):
            ccd /= occd
        return self

    def __str__(self):
        ret = ''
        if self.head is not None: ret += str(self.head)

        ret = '\n\nNumber of CCDs = ' + str(len(self)) + '\n'

        for nccd,ccd in enumerate(self):
            ret += '\nCCD number ' + str(nccd+1) + ':\n'
            ret += str(ccd)
        return ret

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

class Rhead (object):
    """
    Represents essential header info of Ultracam/Ultraspec data read from a
    run###.xml file. 
    """

    def __init__(self, uxml):
        """
        Reads a run###.xml file. UltracamErrors are thrown if some items are not found.
        In some case it will carry on and corresponding attributes are returned as None.
        See below for the list of attributes set.

        Arguments:

         uxml     -- xml file name with format run###.xml

        Attributes set:

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
        node        = udom.getElementsByTagName('instrument_status')[0]
        self.instrument  = node.getElementsByTagName('name')[0].childNodes[0].data
        if self.instrument == 'Ultracam':
            self.instrument = 'UCM'
            self.nxmax, self.nymax = 1080, 1032
        elif self.instrument == 'Ultraspec':
            self.instrument = 'USP'
            self.nxmax, self.nymax = 1056, 1072
        else:
            raise UltracamError('File = ' + self.fname + ', failed to identify instrument.')

        self.application = [nd for nd in node.getElementsByTagName('application_status') \
                                if nd.getAttribute('id') == 'SDSU Exec'][0].getAttribute('name')
        param = {}
        for nd in node.getElementsByTagName('parameter_status'):
            param[nd.getAttribute('name')] = nd.getAttribute('value')

        try:
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
        elif app == 'appl4_frameover_cfg' or 'ap4_frameover':
            self.mode    = 'FFOVER'
        elif app == 'ap9_250_fullframe_mindead' or app == 'ap9_fullframe_mindead' or \
                app == 'appl9_fullframe_mindead_cfg':
            self.mode    = 'FFNCLR'
        elif app == 'ccd201_winbin_con':
            if int(param['X2_SIZE']) == 0:
                self.mode    = '1-USPEC'
            else:
                self.mode    = '2-USPEC'
        else:
            raise UltracamError('File = ' + self.fname + ' failed to identify application = ' + app)

        # binning factors
        self.xbin = int(param['X_BIN_FAC']) if 'X_BIN_FAC' in param else int(param['X_BIN'])
        self.ybin = int(param['Y_BIN_FAC']) if 'Y_BIN_FAC' in param else int(param['Y_BIN'])

        # Windows. For each one store: x & y coords of lower-left pixel, binned dimensions
        self.win = []
        if self.instrument == 'UCM':

            self.speed   = hex(int(param['GAIN_SPEED']))[2:] if 'GAIN_SPEED' in param else None

            if self.mode == 'FFCLR' or self.mode == 'FFNCLR':
                self.win.append((  1, 1, 512//nx, 1024//ny))
                self.win.append((513, 1, 512//nx, 1024//ny))
            elif self.mode == 'FFOVER':
                self.win.append((  1, 1, 540//nx, 1032//ny))
                self.win.append((541, 1, 540//nx, 1032//ny))
            else:
                ystart = int(param['Y1_START'])
                xleft  = int(param['X1L_START'])
                xright = int(param['X1R_START'])
                nx     = int(param['X1_SIZE'])
                ny     = int(param['Y1_SIZE'])
                self.win.append((xleft, ystart, nx, ny))
                self.win.append((xright, ystart, nx, ny))

            if self.mode == '2-PAIR' or self.mode == '3-PAIR':
                ystart = int(param['Y2_START'])
                xleft  = int(param['X2L_START'])
                xright = int(param['X2R_START'])
                nx     = int(param['X2_SIZE'])
                ny     = int(param['Y2_SIZE'])
                self.win.append((xleft, ystart, nx, ny))
                self.win.append((xright, ystart, nx, ny))

            if self.mode == '3-PAIR':
                ystart = int(param['Y3_START'])
                xleft  = int(param['X3L_START'])
                xright = int(param['X3R_START'])
                nx     = int(param['X3_SIZE'])
                ny     = int(param['Y3_SIZE'])
                self.win.append((xleft,ystart,nx,ny))
                self.win.append((xright,ystart,nx,ny))

        elif self.instrument == 'USP':

            self.speed    = ('F' if param['SPEED'] == '0' else \
                                 ('M' if param['SPEED'] == '1' else 'S')) if 'SPEED' in param else None
            self.en_clr   = ('Y' if param['EN_CLR'] == '1' else 'N') if 'EN_CLR' in param else None
            self.hv_gain  = param['HV_GAIN'] if 'HV_GAIN' in param else None
            self.output   = ('N' if param['OUTPUT'] == '0' else 'A') if 'OUTPUT' in param else None

            xstart = int(param['X1_START'])
            ystart = int(param['Y1_START'])
            nx     = int(param['X1_SIZE'])
            ny     = int(param['Y1_SIZE'])
            self.win.append((xstart,ystart,nx,ny))
            
            if self.mode == '2-USPEC':
                xstart = int(param['X2_START'])
                ystart = int(param['Y2_START'])
                nx     = int(param['X2_SIZE'])
                ny     = int(param['Y2_SIZE'])
                self.win.append((xstart,ystart,nx,ny))

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
        """
        Rhead.__init__(self, run + '.xml')
        self._fobj = open(run + '.dat', 'rb')
        self._nf   = nframe
        self._run  = run
        self._flt  = flt
        if nframe != 1:
            self._fobj.seek(self.framesize*(nframe-1))
    
    def __iter__(self):
        """
        Generator to allow Rdata to function as an iterator
        """
        try:
            while 1:
                yield self.__call__()
        except:
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
        Reads the data of frame nframe (starts from 1) and returns a
        corresponding CCD or MCCD object, depending upon the type of data. If
        nframe is None, just reads whatever frame we are on. Raises an
        exception if it fails to read data.  Resets to start of the file in
        this case. The data are stored internally as either 4-byte floats or
        2-byte unsigned ints.

        nframe -- frame number to get, starting at 1. 0 for the last
                  (complete) frame.
        """

        # position read pointer
        self.set(nframe)

        # read timing bytes
        tbytes = self._fobj.read(2*self.headerwords)
        if len(tbytes) != 2*self.headerwords:
            self._fobj.seek(0)
            self._nf = 1
            raise UltracamError('Data.get: failed to read timing bytes')

        # read data
        buff = np.fromfile(self._fobj,'<i2',self.framesize/2-self.headerwords)
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
        if self.instrument == 'UCM':
            # 3 CCDs. Windows come in pairs. Data from equivalent windows come out
            # on a pitch of 6. Some further jiggery-pokery is involved to get the
            # orientation of the frames correct.
            wins1, wins2, wins3 = [],[],[]
            noff = 0
            for left, right in zip(self.win[::2],self.win[1::2]):
                llxl, llyl, nxl, nyl = left
                llxr, llyr, nxr, nyr = right
                npix = 6*nxl*nyl
                if flt:
                    wins1.append(Window(np.reshape(buff[noff:noff+npix:6],(nyl,nxl)),
                                        llxl,llyl,xbin,ybin).astype(np.float32))
                    wins1.append(Window(np.reshape(buff[noff+npix-5:noff:-6],(nyr,nxr)),
                                        llxr,llyr,xbin,ybin).astype(np.float32))
                    wins2.append(Window(np.reshape(buff[noff+2:noff+npix:6],(nyl,nxl)),
                                        llxl,llyl,xbin,ybin).astype(np.float32))
                    wins2.append(Window(np.reshape(buff[noff+npix-3:noff:-6],(nyr,nxr)),
                                        llxr,llyr,xbin,ybin).astype(np.float32))
                    wins3.append(Window(np.reshape(buff[noff+4:noff+npix:6],(nyl,nxl)),
                                        llxl,llyl,xbin,ybin).astype(np.float32))
                    wins3.append(Window(np.reshape(buff[noff+npix-1:noff:-6],(nyr,nxr)),
                                        llxr,llyr,xbin,ybin).astype(np.float32))
                else:
                    wins1.append(Window(np.reshape(buff[noff:noff+npix:6],(nyl,nxl)),
                                        llxl,llyl,xbin,ybin))
                    wins1.append(Window(np.reshape(buff[noff+npix-5:noff:-6],(nyr,nxr)),
                                        llxr,llyr,xbin,ybin))
                    wins2.append(Window(np.reshape(buff[noff+2:noff+npix:6],(nyl,nxl)),
                                        llxl,llyl,xbin,ybin))
                    wins2.append(Window(np.reshape(buff[noff+npix-3:noff:-6],(nyr,nxr)),
                                        llxr,llyr,xbin,ybin))
                    wins3.append(Window(np.reshape(buff[noff+4:noff+npix:6],(nyl,nxl)),
                                        llxl,llyl,xbin,ybin))
                    wins3.append(Window(np.reshape(buff[noff+npix-1:noff:-6],(nyr,nxr)),
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
        

def _check_ucm(fobj):
    """
    Check a file opened for reading in binary mode to see if it is a ucm.

    Returns endian which is a string to be passed
    to later routines indicating endian-ness. 
    """

    # read the format code
    fbytes = fobj.read(4)
    (fcode,) = struct.unpack('i',fbytes)
    if fcode != MAGIC:
        (fcode,) = struct.unpack('>i',fbytes)
        if fcode != MAGIC:
            fobj.close()
            raise CppError('_open_ucm: could not recognise first 4 bytes of ' + fname + ' as a ucm file')
        endian = '>'
    else:
        endian = ''
    return endian

def _rucm(fnoro):
    """
    Read ucm file from disk

    fnoro  -- either a string containing the name of the file to read from ('.ucm' 
              will be appended if necessary), or a file object opened for reading
              in binary mode. The file is closed on exiting the routine.

    Returns head, data, off, xbin, ybin, nxmax, nymax as needed to construct a Ucm
    """    

    # Assume it is a file object, if that fails, assume it is
    # the name of a file.
    try:
        uf = fnoro
        start_format =  _check_ucm(uf)
    except AttributeError, err:
        uf = open(fnoro, 'rb')
        start_format =  _check_ucm(uf)

    # read the header
    lmap = struct.unpack(start_format + 'i', uf.read(4))[0]

    head = Odict()
    for i in xrange(lmap):
        name = cpp.read_string(uf, start_format)
        itype = struct.unpack(start_format + 'i', uf.read(4))[0]
        comment = cpp.read_string(uf, start_format)

        if itype == ITYPE_DOUBLE:
            value = struct.unpack(start_format + 'd', uf.read(8))[0]
        elif itype == ITYPE_INT:
            value = struct.unpack(start_format + 'i', uf.read(4))[0]
        elif itype == ITYPE_UINT:
            value = struct.unpack(start_format + 'I', uf.read(4))[0]
        elif itype == ITYPE_FLOAT:
            value = struct.unpack(start_format + 'f', uf.read(4))[0]
        elif itype == ITYPE_STRING:
            value = cpp.read_string(uf, start_format)
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
            raise UltracamError('ultracam._rucm: do not recognize itype = ' + str(itype))

        # store header information
        head[name] = {'value' : value, 'comment' : comment, 'type' : itype}
        
    # now for the data
    data  = []
    off   = []
        
    # read number of CCDs
    (nccd,) = struct.unpack(start_format + 'i', uf.read(4))

    for nc in range(nccd):
        # read number of wndows
        (nwin,) = struct.unpack(start_format + 'i', uf.read(4))
        ccd  = []
        coff = []
        for nw in range(nwin):
            llx,lly,nx,ny,xbin,ybin,nxmax,nymax = struct.unpack(start_format + '8i', uf.read(32))
            (iout,) = struct.unpack(start_format + 'i', uf.read(4))
            if iout == 0:
                win = numpy.fromfile(file=uf, dtype=numpy.float32, count=nx*ny)
            elif iout == 1:
                win = numpy.fromfile(file=uf, dtype=numpy.uint16, count=nx*ny)
                win = numpy.cast['float32'](win)
            else:
                raise UltracamError('Ultracam data output type iout = ' + str(iout) + ' not recognised')
            win = win.reshape((ny,nx))
            ccd.append(win)
            coff.append((llx, lly))

        data.append(ccd)
        off.append(coff)
    uf.close()

    return head, data, off, xbin, ybin, nxmax, nymax

# Exception class
class UltracamError(Exception):
    """For throwing exceptions from the ultracam module"""
    def __init__(self, value):
        self.value = value
            
    def __str__(self):
        return repr(self.value)





