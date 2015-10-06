"""
Header class for ucm files
"""
from __future__ import absolute_import
from __future__ import print_function

from trm.ultracam.Odict import Odict
from trm.ultracam.Constants import *
from trm.ultracam.UErrors import UltracamError
import six

class Uhead(Odict):
    """
    Class for containing headers compatible with ucm files. Each entry is
    keyed on a string of the form 'User.Filter', the dot signifying a
    hierarchy. The class is a sub-class of an ordered dictionary class to
    allow entries to retain an order.

    Each entry (i.e. what is returned using a key like 'User.Filter') in a
    Uhead consists of a tuple containing a value, a type and a comment, in
    that order. The type corresponds to data types used in ucm files. You can
    simply extract the value from the tuple with index [0], or possibly easier
    to remember, use its 'value' method::

     print uhead['User.Filter'][0]
     print uhead.value('User.Filter')

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

        Args:
          head : a dictionary of key, value pairs, with keys having a
                 hierarchical structure with '.' separators and values
                 each a tuple of (value,type,comment). See add_entry
                 for more.
        """
        Odict.__init__(self)

    def __repr__(self):
        rep = 'Uhead(' + Odict.__repr__(self)[6:]
        return rep

    def __setitem__(self, key, value):
        raise UltracamError('Uhead.__setitem__ disabled to prevent invalid items being defined. Use add_entry')

    def add_entry(self, *args):
        """
        Adds a new Uhead item, checking the various arguments to reduce the
        chances of problems.  This can have either 2 or 4 arguments. The 4
        argument case is as follows:

        Args
          key : hierarchical string of the form 'User.Filter' where 'User'
                is a directory or folder of grouped entries. It cannot have
                blanks and any implied directories must already exists.
                Thus to set a key 'User.Filter.Wheel', 'User.Filter' would
                need to exist and be a directory. The existence of the
                implied 'User' would not be checked in this case, on the
                assumption that it was checked when 'User.Filter' was created.

          value : value to associate (will be ignored in the case of
                  directories, but see the 2 argument case below). The nature
                  of the value varies with the itype; see next.

          itype : one of a range of possible data types. This rather
                  'unpythonic' argument is to address the need to match up
                  with data files and the C++ ULTRACAM pipeline when it
                  comes to writing to disk. Look for integers called
                  'ITYPE_*' to see the set of potential types. The meaning
                  of most data types is obvious. e.g.  ITYPE_DOUBLE or
                  ITYPE_FLOAT expect floating point numbers. In this case
                  both will be stored as a Python float in memory, but
                  will be saved to disk with different numbers of bytes.
                  Less obvious ones are:

                   ITYPE_TIME : the corresponding value should be a two-element
                                tuple or list with first an integer for the
                                number of days and then a float for the
                                number of hours passed.


          comment : comment string with details of the variable.

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

        if not isinstance(key, six.string_types):
            raise UltracamError('Uhead.add_entry: argument "key" must be a string.')

        if not isinstance(comment, six.string_types):
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
            for kold in list(self.keys())[::-1]:
                if dir == kold and self[kold][1] == ITYPE_DIR:
                    break
            else:
                raise UltracamError('Uhead.add_entry: key = ' + key +
                                    ': could not locate directory = ' + key[:ldot])

            # determine position of key within Odict. Must add onto
            # whichever directory it belongs to.
            for index, kold in enumerate(self.keys()):
                if kold.startswith(dir): lind = index

        # the next implicitly check the value: if they can't be converted to
        # the right type, something is wrong.
        print(key, value, itype)
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
            value = (int(value[0]),float(value[1]))
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
        for key, val in six.iteritems(self):
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

if __name__ == '__main__':
    uhead = Uhead()
    uhead.add_entry('User','User information')
    uhead.add_entry('User.Filters','ugi', ITYPE_STRING, 'The filters')
    print('test passed')
