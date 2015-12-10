"""
Class to represent a sub-window of a CCD
"""
from __future__ import absolute_import
from __future__ import print_function

try:
    import numpy as np
except ImportError:
    print('Failed to import numpy; some routines will fail')

try:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    CMDEF = cm.binary
except ImportError:
    print('Failed to import matplotlib; plotting based on it will fail')
    CMDEF = None

try:
    import ppgplot as pg
except ImportError:
    print('Failed to import ppgplot; plotting based on it will fail')

from trm.ultracam.UErrors import UltracamError

class Window(object):
    """
    Class to represent a window of a CCD. Contains an array
    of data along with the following attributes:

     llx, lly     -- lower-left pixels of window
     xbin, ybin   --  pixel binning factors

    Indexed access to the array along the lines of a numpy.ndarray
    so that for example the following work:

    win = Window(data, 1, 2, 3, 4)
    print win[0][1]
    print win[0,1]
    print win[1:-1,1:-1]

    with the last printing out the data values with the outer edge removed.
    """

    def __init__(self, data, llx, lly, xbin, ybin):
        """
        Creates  a Window given some data, a lower-left
        pixel position and binning factors.

        data  -- a 2D numpy.ndarray
        llx   -- X position of left-most pixels in the window.
        lly   -- Y position of lowest pixels in the window.
        xbin  -- X binning factor
        ybin  -- Y binning factor
        """
        if len(data.shape) != 2:
            raise UltracamError('Window._init__: data must be 2D')
        self._data = data
        self.llx  = llx
        self.lly  = lly
        self.xbin = xbin
        self.ybin = ybin

    def __repr__(self):
        rep = 'Window(data=' + repr(self._data) + ', llx=' + \
            repr(self.llx) + ', lly=' + repr(self.lly) + \
            ', xbin=' + repr(self.xbin) + ', ybin=' + \
            repr(self.ybin) + ')'
        return rep

    @property
    def data(self):
        """
        The data array (2D numpy.ndarray)
        """
        return self._data

    @data.setter
    def data(self, dat):
        if len(dat.shape) != 2:
            raise UltracamError('Window.data: dat must be 2D')
        self._data = dat

    @property
    def nx(self):
        """
        The X dimension in binned pixels
        """
        return self._data.shape[1]

    @property
    def ny(self):
        """
        The Y dimension in binned pixels
        """
        return self._data.shape[0]

    @property
    def dtype(self):
        """
        Type of data stored in the array
        """
        return self._data.dtype

    @property
    def size(self):
        """
        Total number of binned pixels
        """
        return self._data.size

    def astype(self, dtype):
        """
        Returns the data as a numpy.ndarray with data type = dtype, (e.g. np.uint16) 
        rounding if converting from a non-integer to an integer type
        """
        if issubclass(dtype, np.integer) and issubclass(self._data.dtype.type,np.integer):
            return self._data.astype(dtype)
        elif issubclass(dtype, np.integer):
            return np.rint(self._data).astype(dtype)
        else:
            return self._data.astype(dtype)

    def totype(self, dtype):
        """
        Converts the internal data to have data type = dtype, 
        rounding if converting from a non-integer to an integer type
        """
        if issubclass(dtype, np.integer) and issubclass(self._data.dtype.type,np.integer):
            self._data = self._data.astype(dtype)
        elif issubclass(dtype, np.integer):
            self._data = np.rint(self._data).astype(dtype)
        else:
            self._data = self._data.astype(dtype)

    def min(self):
        """
        Returns the minimum value of the Window
        """
        return self._data.min()

    def max(self):
        """
        Returns the maximum value of the Window
        """
        return self._data.max()

    def mean(self):
        """
        Returns the mean value of the Window
        """
        return self._data.mean()

    def median(self):
        """
        Returns the median value of the Window
        """
        return np.median(self._data)

    def flatten(self):
        """
        Returns a 1D version of the data of the Window
        """
        return self._data.flatten()

    def sum(self):
        """
        Returns the sum of the data values of the Window
        """
        return self._data.sum()

    def canCropTo(self, other):
        """
        Determines whether the Window has the correct format to be 
        cropped to match other. For this to be the case, the Window 
        has to equal or exceed 'other' in area, have the same binning
        factors, and its pixels must be in step.
        """

        return (self.xbin == other.xbin and self.ybin == other.ybin and
                self.llx <= other.llx and self.lly <= other.lly and 
                self.llx + self.xbin*self.nx >= other.llx + other.xbin*other.nx and
                self.lly + self.ybin*self.ny >= other.lly + other.ybin*other.ny and
                (self.llx - other.llx) % self.xbin == 0 and
                (self.lly - other.lly) % self.ybin == 0)

    def cropTo(self, other):
        """
        Crops the Window to match the format of 'other',
        returning the cropped Window. The Window itself is
        unchanged.
        
        Raises an UltracamError if it is not possible

        other -- the Window to crop to.
        """
        if self.canCropTo(other):
            x1 = (other.llx-self.llx) // self.xbin
            x2 = x1 + other.nx
            y1 = (other.lly-self.lly) // self.ybin
            y2 = y1 + other.ny
            return Window(self._data[y1:y2,x1:x2], other.llx, other.lly, other.xbin, other.ybin)
        else:
            raise UltracamError('Window.cropTo: Window cannot be cropped to "other"') 

    def trim(self, nleft, nright, nbottom, ntop):
        """
        Clips off rows and columns from a Window, and shifts the
        lower-left pixel correspondingly.

        nleft    -- number of columns on left to remove
        nright   -- number of columns on right to remove
        nbottom  -- number of rows on bottom to remove
        ntop     -- number of rows on top to remove
        """
        if nleft != 0 or nright != 0 or nbottom != 0 or ntop != 0:
            if ntop and nright:
                self._data = self._data[nbottom:-ntop,nleft:-nright]
            elif ntop:
                self._data = self._data[nbottom:-ntop,nleft:]
            elif nright:
                self._data = self._data[nbottom:,nleft:-nright]
            else:
                self._data = self._data[nbottom:,nleft:]
            self.llx += nleft*self.xbin
            self.lly += nbottom*self.ybin

    def plot(self, vmin, vmax, mpl=False, cmap=CMDEF, border=True):
        """
        Elementary intensity plot using either matplotlib's imshow
        or pgplot's pggray. Typically some setup may be needed
        before and after this.

        vmin   -- image value for lowest intensity
        vmax   -- image value for highest intensity
        mpl    -- True for matplotlib, otherwise pgplot
        cmap   -- colour map if mpl
        border -- plot a rectangular border around the outermost pixels or not
        """
        if border:
            x1, x2 = self.llx-0.5,self.llx+self.xbin*self.nx-0.5
            y1, y2 = self.lly-0.5,self.lly+self.ybin*self.ny-0.5

        if mpl:
            limits = self.llx-0.5,self.llx+self.xbin*self.nx-0.5,self.lly-0.5,self.lly+self.ybin*self.ny-0.5
            plt.imshow(self._data, cmap=cmap, interpolation='nearest', \
                           vmin=vmin, vmax=vmax, origin='lower', extent=limits)
            if border:
                plt.plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1])
        else:
            tr = np.array([self.llx-self.xbin,self.xbin,0,self.lly-self.ybin,0,self.ybin])
            pg.pggray(self._data,0,self.nx-1,0,self.ny-1,vmax,vmin,tr)
            if border:
                pg.pgline([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1])

    def __getitem__(self, i):
        """
        Returns data[i] where data is the internal ndarray data.

        Example:

        win  = Window(data,10,10,2,2)
        print win[10:-10,10:-10]

        Will print the specified slice of the stored data.
        """
        return self._data[i]

    def __str__(self):
        ret  = self.format() + '\n'
        ret += str(self._data) + ', dtype = ' + str(self._data.dtype)
        return ret

    def __eq__(self, other):
        """
        Tests quality of two Windows. True if the binned dimensions,
        binning factors and lower-left pixels all match.
        """
        return (self.llx == other.llx and self.lly == other.lly and 
                self.xbin == other.xbin and self.ybin == other.ybin and
                self.nx == other.nx and self.ny == other.ny)

    def __ne__(self, other):
        """
        Negation of the equality operator.
        """
        return not self.__eq__(other)

    def __iadd__(self, other):
        if isinstance(other, Window):
            self._data += other._data
        else:
            self._data += other
        return self

    def __isub__(self, other):
        if isinstance(other, Window):
            self._data -= other._data
        else:
            self._data -= other
        return self

    def __imul__(self, other):
        if isinstance(other, Window):
            self._data *= other._data
        else:
            self._data *= other
        return self

    def __idiv__(self, other):
        if isinstance(other, Window):
            self._data /= other._data
        else:
            self._data /= other
        return self

    def __add__(self, other):
        if isinstance(other, Window):
            return Window(self._data + other._data, self.llx, self.lly, self.xbin, self.ybin)
        else:
            return Window(self._data + other, self.llx, self.lly, self.xbin, self.ybin)

    def __sub__(self, other):
        if isinstance(other, Window):
            return Window(self._data - other._data, self.llx, self.lly, self.xbin, self.ybin)
        else:
            return Window(self._data - other, self.llx, self.lly, self.xbin, self.ybin)

    def __mul__(self, other):
        if isinstance(other, Window):
            return Window(self._data * other._data, self.llx, self.lly, self.xbin, self.ybin)
        else:
            return Window(self._data * other, self.llx, self.lly, self.xbin, self.ybin)

    def __div__(self, other):
        if isinstance(other, Window):
            return Window(self._data / other._data, self.llx, self.lly, self.xbin, self.ybin)
        else:
            return Window(self._data / other, self.llx, self.lly, self.xbin, self.ybin)

    def __truediv__(self,other):
        return self.__div__(other)

    def __radd__(self, other):
        return Window(other + self._data, self.llx, self.lly, self.xbin, self.ybin)

    def __rsub__(self, other):
        return Window(other - self._data, self.llx, self.lly, self.xbin, self.ybin)

    def __rmul__(self, other):
        return Window(other * self._data, self.llx, self.lly, self.xbin, self.ybin)

    def __rdiv__(self, other):
        return Window(other / self._data + other, self.llx, self.lly, self.xbin, self.ybin)

    def __rtruediv__(self,other):
        return self.__rdiv__(other)
        
    def format(self):
        """
        Returns a string defining the format of the Window
        """
        return 'llx, lly = ' + str(self.llx) + ', ' + str(self.lly) + \
            '; nx, ny = ' + str(self.nx) + ', ' + str(self.ny) + \
            '; xbin, ybin = ' + str(self.xbin) + ', ' + str(self.ybin)

if __name__ == '__main__':
    win = Window(np.zeros((10,10)),10,15,2,2)
    win.min()
    win.max()
    print('test passed')
