"""
Simple class for representing times
"""

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
