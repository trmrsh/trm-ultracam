#!/usr/bin/env python

descr = \
"""
downloads google spreadsheet worksheets. This downloads all of them 
and saves them to a series of files located in the present working
directory. This should be run every once in a while to ensure all
are saved.
"""

# builtin imports
import sys, argparse, os, time, re
import getpass, datetime

# and some extras
import gdata.docs.service
import gdata.spreadsheet.service
from trm import ultracam

class Wsheet(object):
    """
    Connects to a google spreadsheet to allow the worksheets
    to be downloaded
    """
    
    def __init__(self, email, password, ssheet):
        """
        email    -- your gmail address
        password -- password
        ssheet   -- the spreadsheet to access
        """

        # Connect to the document list service
        self._gd_client = gdata.docs.service.DocsService(source="getsheet")
        self._gd_client.ClientLogin(email, password)

        # connect to the google spreadsheet service
        self._gs_client = gdata.spreadsheet.service.SpreadsheetsService()
        self._gs_client.ClientLogin(email, password)

        # track tokens (not sure this is needed)
        self._docs_auth_token   = self._gd_client.GetClientLoginToken()
        self._sheets_auth_token = self._gs_client.GetClientLoginToken()

        self._gd_client.SetClientLoginToken(self._sheets_auth_token)

        # get document feed
        self._gs_feed = self._gs_client.GetSpreadsheetsFeed()
        for entry in self._gs_feed.entry:
            if entry.title.text == ssheet:
                self._skey = entry.id.text.split('/')[-1]
                print 'Found spreadsheet =',ssheet
                break
        else:
            raise Exception('Wsheet: failed to find spreadsheet = ' + ssheet)

        # get its worksheets
        self._wsheets = self._gs_client.GetWorksheetsFeed(self._skey)

    def __iter__(self):
        """
        Iterator to deliver next worksheet.
        Returns (name,table) where name is the
        name of the worksheet amd table are
        its values as a list of lists.
        """

        for entry in self._wsheets.entry:
            wskey = entry.id.text.split('/')[-1]

            # get column names using a cell feed
            cnames = []
            query = gdata.spreadsheet.service.CellQuery()
            query.max_row = '1'
            query.min_row = '1'
            cfeed = self._gs_client.GetCellsFeed(self._skey, wskey, query=query)
            for centry in cfeed.entry:
                cnames.append(centry.content.text)
            keys = [cname.replace(' ','').lower() for cname in cnames]

            # get list feed to access the table data
            lfeed = self._gs_client.GetListFeed(self._skey, wskey)

            table = [cnames,]
            for lentry in lfeed.entry:
                row = [lentry.custom[key].text for key in keys]
                table.append([ent if ent else '' for ent in row])

            yield (entry.title.text, table)

if __name__ == '__main__':

    # ok deal with arguments first
    parser = argparse.ArgumentParser(description=descr, 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # optional
    parser.add_argument('-e', dest='email', help='gmail address (will be prompted if needed)')
    parser.add_argument('-p', dest='password', help='gmail password (will be prompted if needed)')

    # OK, done with arguments.
    args = parser.parse_args()

    # get user name
    user  = getpass.getuser()

    # list of recognised users
    USERS = {'phsaap' : ('trm', 'tom.r.marsh@gmail.com'), 
             'phslax' : ('mcpb', 'mcpbours@gmail.com'), 
             'phsjaf' : ('eb',   'elme.breedt@gmail.com'), 
             }

    # connect to google spreadsheet
    if args.email:
        email = args.email
    elif user in USERS:
        email = USERS[user][1]
    else:
        email = raw_input('gmail address: ')

    if args.password:
        password = args.password
    else:
        password = getpass.getpass('gmail password: ')

    # connect to the spreadsheet
    wsheet = Wsheet(email, password, 'ULTRACAM run IDs')

    # loop through them, reading and writing to disk
    for name, table in wsheet:
        print 'Found worksheet =',name

        # save to disk
        export_filename = name + '.tsv'

        with open(export_filename,'w') as fout:
            for row in table:
                fout.write('\t'.join(row) + '\n')

        print 'Saved to',export_filename

