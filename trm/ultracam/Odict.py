#!/usr/bin/env python

"""
Ordered dictionary class for trm.ultracam
"""
from __future__ import absolute_import
from __future__ import print_function
import six
from six.moves import map
from six.moves import zip

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
            self._keys = list(dct.keys())
            dict.__init__(self, dct)

    def __repr__(self):
        rep    = 'Odict({'
        first = True
        for k,v in six.iteritems(self):
            if first:
                rep += "'" + k + "': " + repr(v)
                first = False
            else:
                rep += ", '" + k + "': " + repr(v)
        rep += '})'
        return rep

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
        for key, val in six.iteritems(self):
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
        return list(zip(self._keys, list(self.values())))

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
        for key,val in list(dct.items()):
            self.__setitem__(key,val)

    def values(self):
        return list(map(self.get, self._keys))

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

if __name__ == '__main__':
    od = Odict()
    od['B'] = 1.2
    od['A'] = 2.2
    print('test passed')
