#!/usr/bin/env python

"""
Exception classes
"""
from __future__ import print_function

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

class PowerOnOffError(UltracamError):
    """
    Exception for trying to read a power on/off
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

if __name__ == '__main__':

    UltracamError('a test')
    UendError('a test')
    PowerOnOffError('a test')

    print('test passed')

