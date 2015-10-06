"""
ATC server access code extra. See Raw.py for more.
"""
from __future__ import absolute_import
from __future__ import print_function

import os
import urllib2

from trm.ultracam.UErrors import UltracamError

# The ATC FileServer recognises various GET requests
# (look for 'action=' in the code) which are accessed
# using urllib2. The following code is to allow urllib2
# to connect to a local version of the FileServer which
# does not work without this code.

proxy_support = urllib2.ProxyHandler({})
opener = urllib2.build_opener(proxy_support)
urllib2.install_opener(opener)

# Get the URL of the FileServer from the environment
# If not set, a later request to the server will raise
# an UltracamError
URL = os.environ['ULTRACAM_DEFAULT_URL'] if 'ULTRACAM_DEFAULT_URL' in os.environ else None

def get_nframe_from_server(run):
    """
    Returns the number of frames in the run via the FileServer
    """

    if URL is None:
        raise UltracamError('get_nframe_from_server: no url for server found.' +
                            ' Have you set the ULTRACAM_DEFAULT_URL environment variable?')
    # get from FileServer
    full_url = URL + run + '?action=get_num_frames'
    resp = urllib2.urlopen(full_url).read()

    # parse the response
    loc = resp.find('nframes="')
    if loc > -1:
        end = resp[loc+9:].find('"')
        return int(resp[loc+9:loc+9+end])
    else:
        raise UltracamError('get_nframe_from_server: failed to parse server response to ' + full_url)

def get_runs_from_server(dir=None):
    """
    Returns with a list of runs from the server

    dir -- name of sub-directory on server
    """
    if URL is None:
        raise UltracamError('get_runs_from_server: no url for server found.' + 
                            ' Have you set the ULTRACAM_DEFAULT_URL environment variable?')
    # get from FileServer
    if dir is None:
        full_url = URL + '?action=dir'
    else:
        full_url = URL + dir + '?action=dir'
    resp = urllib2.urlopen(full_url).read()

    # parse response from server
    ldir = resp.split('<li>')
    runs = [entry[entry.find('>run')+1:entry.find('>run')+7] for entry in ldir
            if entry.find('getdata">run') > -1]
    runs.sort()
    return runs

if __name__ == '__main__':
    print(get_runs_from_server())
    print('test passed')
