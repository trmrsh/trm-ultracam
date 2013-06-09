#!/usr/bin/env python

"""
enables Python access to ultracam files. The purpose is to facilitate 
development of test routines and the like.

To access an ultracam raw data file, open it with Rdata which then acts as
a source of the data frames using get. e.g.

 rdat  = Rdata('run045')
 frame = rdat(10)
"""

# standard imports
import xml.dom.minidom
import numpy as np
from scipy.stats import scoreatpercentile
import matplotlib.pyplot as plt
import matplotlib.cm     as cm

# my extras
from trm import subs

class Window(np.ndarray):
    """
    Subclass of numpy.ndarray to represent a window of a CCD. The numpy part
    contains the data. The extra attributes are:

     llx, lly     -- lower-left pixels of window
     xbin, ybin   -- x and y binning factors. 
     nxtot, nytot -- total unbinned dimensions of frame
     good         -- are data good? (i.e. do they mean anything). This is to allow for case
                     where a frame is in effect a placeholder such as happend with blue CCD
                     data from ULTRACAM. [True/False]
    """

    def __new__(cls, data, llx, lly, xbin, ybin, nxtot, nytot, good):
        # create the ndarray instance
        obj = np.asarray(data).view(cls)

        # set the extra attributes
        obj.llx   = llx
        obj.lly   = lly
        obj.xbin  = xbin
        obj.ybin  = ybin
        obj.nxtot = nxtot
        obj.nytot = nytot
        obj.good  = good

        # Return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # See the scipy web pages on subclassing ndarray objects
        if obj is None: return

        self.llx   = getattr(obj,'llx',0)
        self.lly   = getattr(obj,'lly',0)
        self.xbin  = getattr(obj,'xbin',0)
        self.ybin  = getattr(obj,'ybin',0)
        self.nxtot = getattr(obj,'nxtot',0)
        self.nytot = getattr(obj,'nytot',0)
        self.good  = getattr(obj,'good',False)

    def __repr__(self):
        st = np.ndarray.__repr__(self)
        st += '[' + str(self.llx) + ',' + str(self.lly) + ',' + str(self.xbin) + \
            ',' + str(self.ybin) + ',' + str(self.nxtot) + ',' + str(self.nytot) + \
            ',' + str(self.good) + ']'
        return st

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        """
        Two windows are equal if they have the same dimensions, the same lower-left corner
        the same binning factors and the same maximum dimensions. Their good-ness status
        flags do not need to match.
        """
        if type(other) is type(self):
            return (self.llx == other.llx and self.lly == other.lly and 
                    self.xbin == other.xbin and self.ybin == other.ybin and 
                    self.nxtot == other.nxtot and self.nytot == other.nytot)
        else:
            return NotImplemented  

    def __neq__(self, other):
        """
        Two windows are only equal if they have the same dimensions, the same lower-left corner
        the same binning factors and the same maximum dimensions. Their good-ness status
        flags do not need to match.
        """
        if type(other) is type(self):
            return (self.llx != other.llx or self.lly != other.lly or 
                    self.xbin != other.xbin or self.ybin != other.ybin or 
                    self.nxtot != other.nxtot and self.nytot != other.nytot)
        else:
            return NotImplemented  

    def limits(self):
        """
        Returns outermost limits of the Window as tuple (xlo,xhi,ylo,yhi)
        """
        ny,nx = self.shape
        return (self.llx-0.5,self.llx+self.xbin*nx-0.5,self.lly-0.5,self.lly+self.ybin*ny-0.5)

class Time(object):
    """
    Object to represent a time. Two attributes:

    mjd   -- modified Julian day, decimal.
    good  -- is the time thought to be reliable?
    """
    def __init__(self, mjd, good):
        self.mjd  = mjd
        self.good = good

class MCCD(subs.Odict):

    """
    Represents multiple CCD frame. The idea is that the different CCDs are
    intrisically linked, e.g. in a multi-arm camera where one always gets an
    image in each of several CCDs.

    MCCD is a sub-class of an Odict (ordered dictionary). The Odict is used to
    store the header while extra attributes store the data in numpy 2D
    arrays. The Odict is keyed on the header item name with '.' used to
    descend into directory hierarchies. So for instance, if you read a file as
    follows:

    flat = trm.ucm.rucm('flat.ucm')
    
    then

    flat['Site.Observatory']['value']

    will return the observatory name.

    MCCD objects have the following extra attributes:

     data  -- the data. data[nc][nw] is a Window object representing window nw of CCD nc (both starting
              from zero, C-style). data[nc] is a list of all windows on CCD nc.

     times -- the Times of each CCD
    """

    def __init__(self, head, data, times=None):
        """
        Creates an MCCD.

        Arguments:

          head  -- the header, an ordered dictionary keyed on the header name with each entry 
                   being a dictionary with the format:

                   {'value' : value, 'comment' : comment, 'type' : itype}

                   The name can be made hierarchical by using '.'. e.g. entries with
                   names 'Detector', 'Detector.Name' and 'Detector.Type' would create
                   a directory and two sub-items in uinfo.

          data  -- list of list of Windows so that data[nc][nw] represents
                   window nw of CCD nc. 

          times -- list of Time objects, one per CCD.

        """

        if head == None:
            subs.Odict.__init__(self)
        else:
            subs.Odict.__init__(self, head)
        
        self.data  = data
        self.times = times

    def __eq__(self, other):
        """
        Equality operator based on file formats: same number of CCDs,
        same number of windows per CCD, same binning factors etc. 
        """

        if type(other) is type(self):

            if self.nccd() != other.nccd(): return False
            for ccd1,ccd2 in zip(self.data,other.data):
                if len(ccd1) != len(ccd2): return False
                for win1,win2 in zip(ccd1,ccd2):
                    if win1 != win2: return False
            return True
        else:
            return NotImplemented

    def __ne__(self, other):
        """
        Inequality operator based on file formats: same number of CCDs,
        same number of windows per CCD, same binning factors etc. 
        """
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def nccd(self):
        "Returns number of CCDs"
        return len(self.data)

    def nwin(self, nc):
        "Returns number of windows of CCD nc (starting from 1)"
        return len(self.data[nc-1])

    def win(self, nc, nw):
        "Returns window number nw of CCD nc (starting from 1)"
        return self.data[nc-1][nw-1]

    def nxy(self, nc, nw):
        "Returns (ny,nx) tuple of pixels dimensions of window number nw of CCD nc (starting from 1)"
        return self.data[nc-1][nw-1].shape

    def wucm(self, fname):
        """
        Writes out to disk in ucm format

        fname  -- file to write to. '.ucm' will be appended if necessary.
        """    

        if not fname.strip().endswith('.ucm'):
            fname = fname.strip() + '.ucm'
        uf = open(fname, 'wb')
    
        # write the format code
        uf.write(struct.pack('i',MAGIC))

        # write the header, starting with the number of entries
        lmap = len(self)
        uf.write(struct.pack('i',lmap))

        for (key,val) in self.iteritems():

            cpp.write_string(uf, key)
            itype = val['type']
            uf.write(struct.pack('i',itype))
            cpp.write_string(uf, val['comment'])

            if itype == ITYPE_DOUBLE:
                uf.write(struct.pack('d', val['value']))
            elif itype == ITYPE_CHAR:
                raise Exception('Hitem: char not enabled')
            elif itype == ITYPE_INT:
                uf.write(struct.pack('i', val['value']))
            elif itype == ITYPE_UINT:
                uf.write(struct.pack('I', val['value']))
            elif itype == ITYPE_LINT:
                raise Exception('Hitem: linit not enabled')
            elif itype == ITYPE_ULINT:
                raise Exception('Hitem: ulint not enabled')
            elif itype == ITYPE_FLOAT:
                uf.write(struct.pack('f', val['value']))
            elif itype == ITYPE_STRING:
                cpp.write_string(uf, val['value'])
            elif itype == ITYPE_BOOL:
                uf.write(struct.pack('B', val['value']))
            elif itype == ITYPE_DIR:
                pass
            elif itype == ITYPE_DATE:
                raise Exception('Hitem: date not enabled')
            elif itype == ITYPE_TIME:
                uf.write(struct.pack('i', val['value'][0]))
                uf.write(struct.pack('d', val['value'][1]))
            elif itype == ITYPE_POSITION:
                raise Exception('Hitem: position not enabled')
            elif itype == ITYPE_DVECTOR:
                uf.write(struct.pack('i', len(val['value'])))
                uf.write(struct.pack(str(len(val['value']))+'d', *val['value']))
            elif itype == ITYPE_UCHAR:
                uf.write(struct.pack('c', val['value']))
            elif itype == ITYPE_TELESCOPE:
                raise Exception('Hitem: telescope not enabled')
            elif itype == ITYPE_USINT:
                uf.write(struct.pack('H', val['value']))
            elif itype == ITYPE_IVECTOR:
                uf.write(struct.pack('i', len(val['value'])))
                uf.write(struct.pack(str(len(val['value']))+'i', *val['value']))
            elif itype == ITYPE_FVECTOR:
                uf.write(struct.pack('i', len(val['value'])))
                uf.write(struct.pack(str(len(val['value']))+'f', *val['value']))
            else:
                raise Exception('Hitem: type =' + str(itype) + 'not recognised')

        # number of CCDs
        nccd = len(self.data)
        uf.write(struct.pack('i', nccd))

        for nc in range(nccd):

            # number of windows
            nwin = len(self.data[nc])
            uf.write(struct.pack('i', nwin))

            for win in self.data[nc]:
                llx     = win.llx
                lly     = win.lly
                ny,nx   = win.shape
                xbin    = win.xbin
                ybin    = win.ybin
                nxtot   = win.nxtot
                nytot   = win.nytot
                iout    = 0
                uf.write(struct.pack('9i',llx,lly,nx,ny,xbin,ybin,nxtot,nytot,iout))
                numpy.cast['float32'](win).tofile(uf)
        uf.close()

    def min(self, nc):
        """
        Returns the minimum value over all windows of a CCD.

        nccd  -- CCD number, starting at 1.
        """
        data = self.data[nc]-1
        minv = data[0].min()
        for dat in data[1:]:
            minv = min(minv, dat.min())
        return minv

    def max(self, nc):
        """
        Returns the maximum value over all windows of a CCD.

        nccd  -- CCD number, starting at 1.
        """
        data = self.data[nc-1]
        minv = data[0].max()
        for dat in data[1:]:
            maxv = max(maxv, dat.max())
        return maxv

    def centile(self, nc, pcent):
        """
        Returns percentile(s). Given a value pcent, this routine returns
        the image level below which pcent percent of the pixel values lie.
        It can also return multiple values for efficiency.

        nccd  -- the CCD of interest. Starts from 1.
        pcent -- percentile or percentiles (array-like)

        Returns image value or values as a list equivalent to the input percentiles
        """

        # check against a string which can look array-like
        if isinstance(pcent, basestring):
            raise Exception('Rdata.centile: argument "pcent" cannot be a string')

        data = self.data[nc-1]

        # generate combined list of all pixels in CCD called 'arr'
        larr = []
        for dat in data:
            larr.append(dat.flatten())
        arr = np.concatenate(larr)

        # distinguish array like or not
        if hasattr(pcent,'__getitem__'):
            scores = []
            for pc in pcent:
                scores.append(scoreatpercentile(arr,pc))
            return scores
        else:
            return scoreatpercentile(arr,pcent)
    
    def bsub(self, nc=0):
        """
        Removes background from each window of a CCD or CCDs. Estimates
        background using a median.

        nc  -- CCD number starting from 1. 0 for all CCDs.
        """
        if nc == 0:
            for dat in self.data:
                for win in dat:
                    win -= np.median(win)
        else:
            for win in self.data[nc-1]:
                win -= np.median(win)
        


    def plot(self, nc, vlo=2., vhi=98., method='p', cmap=cm.binary, show=True):
        """
        Plots a CCD using matplotlib. 

         nccd   -- CCD number (starting from 1)
         vlo    -- number specifying the lowest level to plot (default as a percentile)
         vhi    -- number specifying the lowest level to plot (default as a percentile)
         method -- how vlo and vhi are to be interpreted. 'p' = percentile, 'a' = automatic (min to max,
                  vlo and vhi are irrelevant), 'd' = direct, i.e. just take the values given.
         show   -- show the plot at the end (or not, to allow you to plot something else, use a cursor
                   etc).

        Returns the plot range used.
        """

        plt.axis('equal')
        if method == 'p':
            vmin, vmax = self.centile(nc,(vlo,vhi))
        elif method == 'a':
            vmin, vmax = self.min(), self.max()
        elif method == 'd':
            vmin, vmax = vlo, vhi
        else:
            raise Exception('Rdata.plot: method must be one of p, a or d.')

        for dat in self.data[nc-1]:
            plt.imshow(np.array(dat), cmap=cmap, interpolation='nearest', vmin=vmin, vmax=vmax, \
                           origin='lower', extent=dat.limits())
        plt.xlim(0.5,self.data[nc-1][0].nxtot)
        plt.ylim(0.5,self.data[nc-1][0].nytot)
        if show: plt.show()
        return vmin, vmax

class Uheader (object):
    """
    Represents essential header info of Ultracam/Ultraspec data read from a
    run###.xml file. 
    """

    def __init__(self, uxml):
        """
        Reads a run###.xml file. Exceptions are thrown if some items are not found.
        In some case it will carry on and corresponding attributes are returned as None.
        See below for the list of attributes set.

        Arguments:

         uxml     -- xml file name with format run###.xml

        Attributes set:

         application -- data acqusition application template name.
         fname       -- file used to define the format.
         framesize   -- total number of bytes per frame.
         headerwords -- number of words (2-bytes/word) in timing info at start of a frame.
         instrument  -- instrument name.
         mode        -- a standardised summary of the readout mode derived from the application name.
         speed       -- readout speed.
         user        -- dictionary of user information. Set to None if there was none found.
         win         -- A list of Window objects, one per window. ULTRACAM data is multi-CCD
                        but the windows of each CCD are identical so the information is only stored 
                        once for all CCDs.
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
            NXTOT, NYTOT = 1080, 1032
        elif self.instrument == 'Ultraspec':
            self.instrument = 'USP'
            NXTOT, NYTOT = 1056, 1072
        else:
            raise Exception('File = ' + self.fname + ', failed to identify instrument.')

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
        elif app == 'ap9_250_fullframe_mindead' or app == 'ap9_fullframe_mindead' or app == 'appl9_fullframe_mindead_cfg':
            self.mode    = 'FFNCLR'
        elif app == 'ccd201_winbin_con':
            if int(param['X2_SIZE']) == 0:
                self.mode    = '1-USPEC'
            else:
                self.mode    = '2-USPEC'
        else:
            raise Exception('File = ' + self.fname + ' failed to identify application = ' + app)

        # binning factors
        xbin = int(param['X_BIN_FAC']) if 'X_BIN_FAC' in param else int(param['X_BIN'])
        ybin = int(param['Y_BIN_FAC']) if 'Y_BIN_FAC' in param else int(param['Y_BIN'])

        # Windows. For each one store:
        #
        # x & y coords of lower-left pixel, binned dimensions, binning factors, dimensions
        # of CCD.
        self.win = []
        if self.instrument == 'UCM':

            self.speed   = hex(int(param['GAIN_SPEED']))[2:] if 'GAIN_SPEED' in param else None

            if self.mode == 'FFCLR' or self.mode == 'FFNCLR':
                self.win.append((  1, 1, 512//nx, 1024//ny, xbin, ybin, NXTOT, NYTOT))
                self.win.append((513, 1, 512//nx, 1024//ny, xbin, ybin, NXTOT, NYTOT))
            elif self.mode == 'FFOVER':
                self.win.append((  1, 1, 540//nx, 1032//ny, xbin, ybin, NXTOT, NYTOT))
                self.win.append((541, 1, 540//nx, 1032//ny, xbin, ybin, NXTOT, NYTOT))
            else:
                ystart = int(param['Y1_START'])
                xleft  = int(param['X1L_START'])
                xright = int(param['X1R_START'])
                nx     = int(param['X1_SIZE'])
                ny     = int(param['Y1_SIZE'])
                self.win.append((xleft, ystart, nx, ny, xbin, ybin, NXTOT, NYTOT))
                self.win.append((xright, ystart, nx, ny, xbin, ybin, NXTOT, NYTOT))

            if self.mode == '2-PAIR' or self.mode == '3-PAIR':
                ystart = int(param['Y2_START'])
                xleft  = int(param['X2L_START'])
                xright = int(param['X2R_START'])
                nx     = int(param['X2_SIZE'])
                ny     = int(param['Y2_SIZE'])
                self.win.append((xleft, ystart, nx, ny, NXTOT, NYTOT))
                self.win.append((xright, ystart, nx, ny, NXTOT, NYTOT))

            if self.mode == '3-PAIR':
                ystart = int(param['Y3_START'])
                xleft  = int(param['X3L_START'])
                xright = int(param['X3R_START'])
                nx     = int(param['X3_SIZE'])
                ny     = int(param['Y3_SIZE'])
                self.win.append((xleft,ystart,nx,ny,NXTOT,NYTOT))
                self.win.append((xright,ystart,nx,ny,NXTOT,NYTOT))

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
            self.win.append((xstart,ystart,nx,ny,NXTOT,NYTOT))
            
            if self.mode == '2-USPEC':
                xstart = int(param['X2_START'])
                ystart = int(param['Y2_START'])
                nx     = int(param['X2_SIZE'])
                ny     = int(param['Y2_SIZE'])
                self.win.append((xstart,ystart,nx,ny,NXTOT,NYTOT))

class Rdata (Uheader):
    """
    Class to represent an Ultracam raw data file. The idea is to open
    the file, and then the object generated can be used to deliver frames 
    by specifying a frame number. An internal frame number is kept to say 
    which frame we are on. This is incremented when a frame is read.

    If accurate times are wanted, it is usually necessary (depending upon
    the readout mode) to read frames sequentially in order to generate a 
    suitable buffer of timestamps from which the time can be deduced. Rdata 
    objects internally store a buffer of timestamps. If this is not done,
    time will be flagged as unreliable.
    """
    def __init__(self, run):
        """
        Connects to a raw data file. The file is kept open. The file pointer
        is set to the start of the file (frame number 1)

        Arguments:

        run -- as in 'run036'. Will try to access equivalent .xml and .dat
               files

        """
        Uheader.__init__(self, run + '.xml')
        self._fobj = open(run + '.dat', 'rb')
        self._nf   = 1

    def get(self, nframe=None):
        """
        Reads the data of frame nframe (starts from 1) and returns a
        corresponding MCCD object. If nframe is None, just reads whatever
        frame we are on. An internal counter is used to keep track of which
        frame we are set to read. Raises Exceptions if it fails to read data.
        Resets to start of the file in this case.

        nframe  -- frame number to get, starting at 1. 0 for the last (complete) frame.
        """

        # position read pointer
        if nframe is not None:
            if nframe < 0:
                raise Exception('Data.get: nframe < 0')
            elif nframe == 0:
                self._fobj.seek(0,2)
                fp = self._fobj.tell() 
                nf = fp // self.framesize
                self._fobj.seek(self.framesize*(nf-1)-fp,2)
                self._nf = nf
            elif self._nf != nframe:
                self._fobj.seek(self.framesize*(nframe-1))
                self._nf = nframe

        # read timing bytes
        tbytes = self._fobj.read(2*self.headerwords)
        if len(tbytes) != 2*self.headerwords:
            self._fobj.seek(0)
            self._nf = 1
            raise Exception('Data.get: failed to read timing bytes')

        # set frame flag
        dok = True

        # read data
        buff = np.fromfile(self._fobj,'<i2',self.framesize/2-self.headerwords)
        if len(buff) != self.framesize/2-self.headerwords:
            self._fobj.seek(0)
            self._nf = 1
            raise Exception('Data.get: failed to read data')

        # interpret data
        if self.instrument == 'UCM':
            # 3 CCDs. Windows come in pairs. Data from equivalent windows come out
            # on a pitch of 6. Some further jiggery-pokery is involved to get the
            # orientation of the frames correct.
            data = [[],[],[]]
            noff = 0
            for left, right in zip(self.win[::2],self.win[1::2]):
                llxl, llyl, nxl, nyl, xbinl, ybinl, nxtotl, nytotl = left
                llxr, llyr, nxr, nyr, xbinr, ybinr, nxtotr, nytotr = right
                npix   = 6*nxl*nyl
                data[0].append(Window(np.reshape(buff[noff:noff+npix:6],(nyl,nxl)),
                                      llxl,llyl,xbinl,ybinl,nxtotl,nytotl,dok))
                data[0].append(Window(np.reshape(buff[noff+npix-5:noff:-6],(nyr,nxr)),
                                      llxr,llyr,xbinr,ybinr,nxtotr,nytotr,dok))
                data[1].append(Window(np.reshape(buff[noff:noff+npix:6],(nyl,nxl)),
                                      llxl,llyl,xbinl,ybinl,nxtotl,nytotl,dok))
                data[1].append(Window(np.reshape(buff[noff+npix-3:noff:-6],(nyr,nxr)),
                                      llxr,llyr,xbinr,ybinr,nxtotr,nytotr,dok))
                data[2].append(Window(np.reshape(buff[noff:noff+npix:6],(nyl,nxl)),
                                      llxl,llyl,xbinl,ybinl,nxtotl,nytotl,dok))
                data[2].append(Window(np.reshape(buff[noff+npix-1:noff:-6],(nyr,nxr)),
                                      llxr,llyr,xbinr,ybinr,nxtotr,nytotr,dok))
                noff += npix
        else:
            raise Exception('Have not implemented anything for ' + self.instrument)
            
        # build header
        head = subs.Odict()
        head['User'] = 'Data entered by user at telescope'
        head['Instrument'] = 'Instrument information'
        head['Instrument.instrument']    = self.instrument
        head['Instrument.headwerwords']  = self.headerwords
        head['Instrument.framesize']     = self.framesize

        # move frame counter on by one
        self._nf += 1
        
        return MCCD(head, data)

    def nframe(self):
        """
        Returns number of complete frames in data file
        """
        self._fobj.seek(0,2)
        nf = self._fobj.tell() // self.framesize
        self._fobj.seek(self.framesize*(self._nf-1))
        return nf

    def next(self):
        """
        Returns next frame number to be read if reading
        sequentially (starts at 1)
        """
        return self._nf
        

# Test section if this is run as a script

if __name__ == '__main__':

    # link to a data file
    rdat = Rdata('../../data/run036')

    # read in a file
    mccd = rdat.get()

    # remove background
    mccd.bsub()

    # plot CCD 2
    mccd.plot(2)
    mccd.plot(1)

    print 'File has',rdat.nframe(),'frames.'

    print 'Next file to read =',rdat.next()
