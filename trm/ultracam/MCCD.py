"""
Classes to represent multiple CCDs
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import ppgplot as pg

from trm.ultracam.Constants import *
from trm.ultracam.Uhead import Uhead
from trm.ultracam.Window import Window
from trm.ultracam.CCD import CCD
from trm.ultracam.Utils import write_string, read_string, check_ucm
from trm.ultracam.UErrors import UltracamError

class MCCD(object):
    """
    Represents multiple CCD frame. The idea is that one has an instrument
    which produces multiple CCDs of data that are intrinsically linked, 
    e.g. in a multi-arm camera where one always gets an image in each of 
    several CCDs. This is to represent ULTRACAM data. There is some common 
    header information, plus the data, represented by a list of CCD objects.

    Indexed access returns the component CCD objects.
    """

    def __init__(self, data, head):
        """
        Creates an MCCD.

        Arguments:

          data  -- list of CCD objects

          head -- the header, either None or a Uhead object.
          
        Sets the equivalent attribute 'head'
        """
        for ccd in data:
            if not isinstance(ccd, CCD):
                raise UltracamError('MCCC.__init__: one or more of the elements of data is not a CCD.')

        if head is not None and not isinstance(head, Uhead):
            raise UltracamError('MCCC.__init__: head should be a Uhead (or None).')

        self._data = data
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
        start_format = check_ucm(uf)

        # read the header
        lmap = struct.unpack(start_format + 'i', uf.read(4))[0]

        head = Uhead()
        for i in xrange(lmap):
            name    = read_string(uf, start_format)
            itype   = struct.unpack(start_format + 'i', uf.read(4))[0]
            comment = read_string(uf, start_format)

            if itype == ITYPE_DOUBLE:
                value = struct.unpack(start_format + 'd', uf.read(8))[0]
            elif itype == ITYPE_INT:
                value = struct.unpack(start_format + 'i', uf.read(4))[0]
            elif itype == ITYPE_UINT:
                value = struct.unpack(start_format + 'I', uf.read(4))[0]
            elif itype == ITYPE_FLOAT:
                value = struct.unpack(start_format + 'f', uf.read(4))[0]
            elif itype == ITYPE_STRING:
                value = read_string(uf, start_format)
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

    def __len__(self):
        """
        Returns number of CCDs
        """
        return len(self._data)

    def __getitem__(self, i):
        """
        Returns data[i] where data is the internal ndarray data.
        """
        return self._data[i]

    def __setitem__(self, i, ccd):
        """
        Sets the i-th Window
        """
        if not isinstance(ccd, CCD):
            raise UltracamError('MCCD.__setitem__: ccd must be a CCD')
        self._data[i] = ccd

    @property
    def data(self):
        """
        The list of CCDs
        """
        return self._data

    @data.setter
    def data(self, ccds):
        for ccd in ccds:
            if not isinstance(ccd, CCD):
                raise UltracamError('MCCD.data: ccds must be a list of CCDs.')
        self._data = ccds

    @property
    def nccd(self):
        """
        The number of CCDs
        """
        return len(self._data)

    def cropTo(self, mccd):
        """
        Crops the MCCD to match mccd if possible, returns the cropped
        MCCD. Raises an UltracamError if it does not succeed.
        """
        if len(self) != len(mccd):
            raise UltracamError('MCCD.crop: number of CCDs did not match')

        ccds = []
        for ccd, ccdo in zip(self._data, mccd._data):
            ccds.append(ccd.cropTo(ccdo))
        return MCCD(ccds, self.head)

    def __eq__(self, other):
        """
        Equality operator tests same number of CCDs and that each CCD matches.
        """
        
        if len(self) != len(other): return False

        for sccd, occd in zip(self._data,other._data):
            if sccd != occd: 
                return False
        return True

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

        for ccd in self._data:

            # number of windows
            nwin = len(ccd)
            uf.write(struct.pack('i', nwin))

            for win in ccd:
                uf.write(struct.pack('9i',win.llx,win.lly,win.nx,win.ny,win.xbin,win.ybin,
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
            for ccd in self._data:
                ccd.rback()
        else:
            self._data[nc].rback()

    def anyInt(self):
        """
        Returns True if any of the contributing CCDs are based on integers. It can be
        useful for memory and disk space reasons to keep data as 2-byte unsigned integers 
        but cause problems with arithematic operations. This allows you to check. 
        See also 'anyFloat', 'toFloat' and 'toInt'.
        """
        for ccd in self._data:
            if ccd.anyInt(): return True
        return False

    def anyFloat(self):
        """
        Returns True if any of the contributing CCDs are based on floats. This is needed
        to evaluate the output type needed when writing to disk.
        """
        for ccd in self._data:
            if ccd.anyFloat(): return True
        return False

    def toFloat(self, single=True):
        """
        Converts all CCDs to a float type, either single or double
        precision.

        single  -- True to convert to 4-byte floats (else 8-byte)
        """
        for ccd in self._data:
            ccd.toFloat(single)

    def toInt(self):
        """
        Converts all CCDs to an unsigned 2-byte integer type, rounding
        to the nearest integer. Warnings will be issued if data lies outside
        the 0 to 65535 range, but the conversion will proceed.
        """
        for ccd in self._data:
            ccd.toInt()

    def max(self):
        """
        Returns a tuple of maximum values, 1 per CCD.
        """
        mx = []
        for ccd in self._data:
            mx.append(ccd.max())
        return tuple(mx)

    def min(self):
        """
        Returns a tuple of minimum values, 1 per CCD.
        """
        mn = []
        for ccd in self._data:
            mn.append(ccd.min())
        return tuple(mn)

    def mean(self):
        """
        Returns a tuple of mean values, 1 per CCD.
        """
        mn = []
        for ccd in self._data:
            mn.append(ccd.mean())
        return tuple(mn)

    def median(self):
        """
        Returns a tuple of median values, 1 per CCD.
        """
        mn = []
        for ccd in self._data:
            mn.append(ccd.median())
        return tuple(mn)

    # arithematic
    def __iadd__(self, other):
        """
        Adds 'other' to the MCCD in place (+=)
        """
        if isinstance(other, MCCD):
            for ccd,occd in zip(self._data,other._data):
                ccd += occd
        else:
            for ccd in self._data:
                ccd += other
        return self

    def __isub__(self, other):
        """
        Subtracts 'other' from the MCCD in place (-=)
        """
        if isinstance(other, MCCD):
            for ccd,occd in zip(self._data,other._data):
                ccd -= occd
        else:
            for ccd in self._data:
                ccd -= other
        return self

    def __imul__(self, other):
        """
        Multiplies the MCCD by 'other' in place (\*=)
        """
        if isinstance(other, MCCD):
            for ccd,occd in zip(self._data,other._data):
                ccd *= occd
        else:
            for ccd in self._data:
                ccd *= other
        return self

    def __idiv__(self, other):
        """
        Divides the MCCD by 'other' in place (/=)
        """
        if isinstance(other, MCCD):
            for ccd,occd in zip(self._data,other._data):
                ccd /= occd
        else:
            for ccd in self._data:
                ccd /= other
        return self

    def __str__(self):
        ret = ''
        if self.head is not None: ret += str(self.head)
        ret += '\n\nNumber of CCDs = ' + str(len(self)) + '\n'

        for nccd,ccd in enumerate(self._data):
            ret += '\nCCD number ' + str(nccd+1) + ':\n'
            ret += str(ccd)
        return ret

    def __add__(self, other):
        """
        Adds 'other' to the MCCD (+)
        """
        tccd = []
        if isinstance(other, MCCD):
            for ccd,occd in zip(self._data,other._data):
                tccd.append(ccd + occd)
        else:
            for ccd in self._data:
                tccd.append(ccd + other)
        return MCCD(tccd, self.head)

    def __sub__(self, other):
        """
        Subtract 'other' from the MCCD (-)
        """
        tccd = []
        if isinstance(other, MCCD):
            for ccd,occd in zip(self._data,other._data):
                tccd.append(ccd - occd)
        else:
            for ccd in self._data:
                tccd.append(ccd - other)
        return MCCD(tccd, self.head)

    def __mul__(self, other):
        """
        Multiply 'other' by the MCCD (*)
        """
        tccd = []
        if isinstance(other, MCCD):
            for ccd,occd in zip(self._data,other._data):
                tccd.append(ccd * occd)
        else:
            for ccd in self._data:
                tccd.append(ccd * other)
        return MCCD(tccd, self.head)

    def __div__(self, other):
        """
        Divide MCCD by 'other' from the MCCD (/)
        """
        tccd = []
        if isinstance(other, MCCD):
            for ccd,occd in zip(self._data,other._data):
                tccd.append(ccd / occd)
        else:
            for ccd in self._data:
                tccd.append(ccd / other)
        return MCCD(tccd, self.head)
        
    def __truediv__(self,other):
    	"""
    	Divide MCCD by 'other' from the MCCD (/) when future division used
    	"""
    	return self.__div__(other)

    def __radd__(self, other):
        """
        Returns other + MCCD (an MCCD)
        """
        tccds = []
        for ccd in self._data:
            tccds.append(other + ccd)
        return MCCD(tccds, self.head)

    def __rsub__(self, other):
        """
        Returns other - MCCD (an MCCD)
        """
        tccds = []
        for ccd in self._data:
            tccds.append(other - ccd)
        return MCCD(tccds, self.head)

    def __rmul__(self, other):
        """
        Returns other * MCCD (an MCCD)
        """
        tccds = []
        for ccd in self._data:
            tccds.append(other * ccd)
        return MCCD(tccds, self.head)

    def __rdiv__(self, other):
        """
        Returns other / MCCD (an MCCD)
        """
        tccds = []
        for ccd in self._data:
            tccds.append(other / ccd)
        return MCCD(tccds, self.head)

    def plot(self, vlo=2., vhi=98., nc=-1, method='p', mpl=False, cmap=cm.binary, \
                 close=True, x1=None, x2=None, y1=None, y2=None, sepmin=1.):
        """
        Plots an MCCD using pgplot or matplotlib if preferred.

        :Parameters:
         vlo : float
            number specifying the lowest level to plot (default as a percentile)
         vhi : float
            number specifying the lowest level to plot (default as a percentile)
         nc : int
            CCD number (starting from 0, -1 for all)
         method : string
            how vlo and vhi are to be interpreted. 'p' = percentile, 'a' = automatic (min to max,
            vlo and vhi are irrelevant), 'd' = direct, i.e. just take the values given.
         mpl : bool
            True to prefer matplotlib over pgplot (which may not even be an option)
         cmap : matplotlib.cm.binary
            colour map if using matplotlib
         close : bool
            close (pgplot) or 'show' (matplotlib) the plot at the end (or not, to allow 
            you to plot something else, use a cursor etc). In the case of pgplot, this also
            implies opening the plot at the start, i.e. a self-contained quick plot.
         x1 : float
            left-hand plot limit. Defaults to 0.5
         x2 : float
             right-hand plot limit. Defaults to nxmax+0.5
         y1 : float
            lower plot limit. Defaults to 0.5
         y2 : float
             upper plot limit. Defaults to nymax+0.5
         sepmin : float
             minimum separation between intensity limits (> 0 to stop PGPLOT complaining)

        :Returns:
         range(s) : tuple or list
            the plot range(s) used either as a single 2-element tuple, or
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
            if nc2-nc1 > 1: pg.pgsubp(nc2-nc1,1)

        prange = []
        for nc, ccd in enumerate(self._data[nc1:nc2]):

            # Determine intensity range to display
            if method == 'p':
                vmin, vmax = ccd.centile((vlo,vhi))
            elif method == 'a':
                vmin, vmax = ccd.min(), ccd.max()
            elif method == 'd':
                vmin, vmax = vlo, vhi
            else:
                raise UltracamError('MCCD.plot: method must be one of p, a or d.')

            if vmin == vmax:
                vmin -= sepmin/2.
                vmax += sepmin/2.
            prange.append((vmin, vmax))

            # start
            nxmax, nymax = ccd.nxmax, ccd.nymax
            x1    = 0.5 if x1 is None else x1
            x2    = nxmax+0.5 if x2 is None else x2
            y1    = 0.5 if y1 is None else y1
            y2    = nymax+0.5 if y2 is None else y2

            if mpl:
                if nc2-nc1 > 1:
                    plt.subplot(1,nc2-nc1,nc+1)
                plt.axis('equal')
            else:
                if nc2-nc1 > 1: pg.pgpanl(nc-nc1+1,1)
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

    def format(self):
        """
        Returns a string describing the format of the MCCD
        """
        ret = ''
        for nccd, ccd in enumerate(self._data):
            ret += 'CCD ' + str(nccd+1) + ':\n'
            ret += ccd.format()+'\n'
        return ret

class UCAM(MCCD):
    """
    Specialised version of an MCCD which represents ULTRACAM frames.
    Essentially this means that it expects some ULTRACAM-specific 
    features that cannot be assumed for general MCCDs.
    """
    def __init__(self, data, head):
        """
        Imposes some restrictions on the inputs not set by MCCD. There must be
        3 CCDs, and each must have an even number of Windows. 
        """
        if len(data) != 3:
            raise UltracamError('UCAM.__init__: require list of 3 CCDs for data')

        for ccd in data:
            if ccd.nwin % 2 != 0:
                raise UltracamError('UCAM.__init__: all CCDs must have an even number of Windows')

        MCCD.__init__(self, data, head)

    def checkData(self):
        """
        Checks the CCDs for some standard problems with ULTRACAM data. For each
        CCD it returns with a flag that is True if the corresponding CCD is thought
        to have a problem, along with a string describing the issue. This should
        only be applied to raw ULTRACAM images, not ones that have been processed
        in any way and it will only work if they have been read in as integers.

        The problems this routine picks up are:

         -- too low overall level (judged from many previous images)
         -- too high an overall level (loosely based on when peppering occurs)
         -- left, right windows differ too much which sometimes happens
         -- too many pixels with the same value

        It will only report one problem per CCD at most and breaks off checking as 
        soon as it finds an issue, thus it is possible that there are other problems.
        It is chiefly here to enable warnings of possible problems. It tries to
        find the worst problems first.

        Returns a list of tuples, one for each CCD. Each of these consists of a 
        True/False flag and a string. Thus the following code makes sense:

        r,g,b = mccd.checkData()
        if r[0]: print 'Red CCD has a problem: ',r[1]
        if g[0]: print 'Green CCD has a problem: ',g[1]
        if b[0]: print 'Blue CCD has a problem: ',b[1]

        """
        if self[0][0].dtype != np.uint16:
            raise UltracamError('UCAM.checkData: only works with raw unsigned 16-bit int images')

        ret = []
        for nc, ccd in enumerate(self._data):
           for winl, winr in zip(ccd[::2],ccd[1::2]):

               # check the modes
               hist  = np.bincount(winl.flatten())
               nmode = hist[np.argmax(hist)]
               if nmode > winl.size // 4:
                   ret.append((True,'a window has >25% pixels of same value'))
                   break

               hist = np.bincount(winr.flatten())
               nmode = hist[np.argmax(hist)]
               if nmode > winl.size // 4:
                   ret.append((True,'a window has >25% pixels of same value'))
                   break

               # check the medians
               l = winl.median()
               r = winr.median()
               
               if nc == 0:
                   if l < 1580 or r < 1580:
                       ret.append((True,'too low'))
                       break
                   elif l > 55000 or r > 55000:
                       ret.append((True,'too high'))
                       break
                   elif abs(r-l-70) > 30+0.05*max(0,l - 1700):
                       ret.append((True,'left and right too different'))
                       break

               elif nc == 1:
                   if l < 1200 or r < 1200:
                       ret.append((True,'too low'))
                       break
                   elif l > 30000 or r > 30000:
                       ret.append((True,'too high'))
                       break
                   elif abs(r-l-10) > 60+0.05*max(0,l - 1300):
                       ret.append((True,'left and right too different'))
                       break

               elif nc == 2:
                   if l < 1000 or r < 1000:
                       ret.append((True,'too low'))
                       break
                   elif l > 30000 or r > 30000:
                       ret.append((True,'too high'))
                       break
                   elif abs(r-l-100+5.5*l/60.) > 70+0.05*max(0,l - 1500):
                       ret.append((True,'left and right too different'))
                       break

           else:
               # loop traversed without a problem
               ret.append((False,''))

        return ret

if __name__ == '__main__':
    import trm.ultracam.Time as Time
    uhead = Uhead()
    uhead.add_entry('User','User information')    
    uhead.add_entry('User.Filters','ugi', ITYPE_STRING, 'The filters')

    win1 = Window(np.zeros((2,2)),1,2,2,2)    
    win2 = Window(np.zeros((3,3)),100,2,2,2)

    time = Time(55000.2, 20., True, '')

    ccd1 = CCD([win1,win2], time, 1024, 1024, True, uhead)
    ccd2 = CCD([win1+20.,win2+100.], time, 1024, 1024, True, uhead)
    mccd = MCCD([ccd1,ccd2], uhead)
    print 'test passed'

