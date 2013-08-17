#
# Gets loaded by trm.ultracam
#

import re
from trm.ultracam.UErrors import UltracamError

class Log(object):
    """
    Class to read and store log file data. These come in two formats:

    1) Old style: run, target name, filters, comment
    2) New style: run, comment (target names are in the xml files)

    The class just stores the data in dictionaries 'comment',
    'target', and 'filters'; 'format' is an integer specifying 
    the format as above. 'target' and 'filters' are blank in 
    the case of format == 2. It is assumed in this case that 
    they are present in corresponding xml files. The dictionaries
    are keyed on the run id, i.e. 'run005'
    """

    def __init__(self, fname):
        """
        Constructs a new Log given a file. Makes empty
        dictionaries if none found and reports an error
        """
        self.format  = 2
        self.target  = {}
        self.filters = {}
        self.comment = {}

        rec    = re.compile('file\s+object\s+filter', re.I)
        old    = re.compile('\s*(\S+)\s+(\S+)\s+(.*)$')
        oldii  = re.compile('\s*(\S+)\s*$')
        f  = open(fname)
        for line in f:
            m = rec.search(line)
            if m:
                self.format = 1
                if len(self.comment):
                    raise UltracamError('Error in night log = ' + fname + ', line = ' + line)

            if line.startswith('run'):
                run = line[:6]
                if self.format == 2:
                    self.comment[run] = line[6:].strip()
                else:
                    m = old.search(line[6:])
                    if m:
                        self.target[run]  = m.group(1)
                        self.filters[run] = m.group(2)
                        self.comment[run] = m.group(3)
                    else:
                        m = oldii.search(line[6:])
                        if m:
                            self.target[run]  = m.group(1)
                        else:
                            self.target[run]  = 'UNKNOWN'
                        self.filters[run] = 'UNKNOWN'
                        self.comment[run] = ''

    
