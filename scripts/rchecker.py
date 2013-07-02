#!/usr/bin/env python

descr = \
"""
For ID-ing ULTRACAM runs. Sequentially displays frames of all runs of a
specified night prompting for user input to define a data type, target 
name, the filters and a comment for each run. Try option -hh when you 
first run the script.
"""

# builtin imports
import sys, argparse, os, time, urllib, urllib2, re
import getpass, datetime, textwrap

# and some extras
from ppgplot import *
from trm import ultracam

hhelp = \
"""
The program is designed to simplify identification of ULTRACAM run types. When
you run it, your job will be to make 4 decisions for every run: (1) what type
of data it represents, (2) what target, (3) what filters, (4) what comment to
make. Default values for (2) and (3) are taken from logs. The script must be
run from a directory containing the night-by-night directories. These can be
written in the form YYYY-MM-DD or YYYY_MM_DD. The directories must contain
.xml and .dat files and also the hand observing logs either YYYY_MM_DD_log.dat
or my possibly-edited copy, YYYY-MM-DD.dat.

Before working through a night, take good look at the logs to try to establish
whether and when any filter changes have occurred. Write down the runs on
which the changes have occurred. It is quite easy for mistakes to occur, so
take care if there have been changes.

The program stores the results in a google spreadsheet called "ULTRACAM run
IDs" or in a local file. The latter can be added later manually using the
import mechanism on the google spreadsheet, appending to what is there already
and using "|" as a custom separation character. The latter means you should
not have any "|" characters in your input or everything will become hopelessly
confused. If you can manage it, the direct upload method is to be preferred.

To use the google spreadsheet upload mechanism will require my sharing it with
you and the installation of the google data API, "gdata". It also requires you
to supply your gmail account info which cannot be hard-coded in the script for
obvious reasons. You can either let it prompt you or supply through -e and -p
options on the command-line. Since the script runs through all runs of a night
this is a negligible overhead.

An approximate bias subtraction is carried out as follows: the date of the
night is used to estimate the bias level based upon past runs of similar date
and the same readout speed "cdd" etc. This avoids arbitrarily removing the
median of the frames or having to identify a suitable bias. If you prefer the
direct data, use the -b switch to turn off the background subtraction. Please
report to me if you ever get warnings about having failed to identify a
default bias value.

The final frame of most ULTRACAM runs is shortened and often not a good
example of the entire run. By default the script will therefore ask if you
want to display it. If you find yourself always hitting <cr> to skip it, then
the -c switch will skip the request altogether and not show the final
frame. The -n switch on the other hand will not ask, but will show it (unless
overridden by -c).

For details of all options, use the -h option.
"""

def pLogs(oneLiners, nr, width, ndelta=None):
    """
    Prints out logs around the current run highlighting the
    current run with '>' and wrapping long comment lines.
    
    oneLiners  -- List of dictionaries in run order. Each
                  dictionary specifies the run, target, filters
                  and a comment.

    nr         -- Index of the current run

    width      -- Number of characters to wrap comments with.

    ndelta     -- Parameter to show how many to display 
                  around the current run. None for the whole 
                  lot, otherwise it starts at nr-ndelta and 
                  ends at nr+ndelta.
    """

    if not hasattr(pLogs,'textWrap'):
        pLogs.textWrap = textwrap.TextWrapper(width=width)

    if ndelta is None:
        nr_start = 0
        nr_end   = len(oneLiners)
    else:
        nr_start = max(0, nr - ndelta)
        nr_end   = min(len(runs), nr_start + 2*ndelta + 1)

    for nrr in xrange(nr_start, nr_end):
        ol = oneLiners[nrr]
        if nrr == nr:
            pLogs.textWrap.initial_indent = \
            '>%6s %-20s %-10s ' % (ol['run'],
                                   ol['target'] if ol['target'] else 'UNKNOWN',
                                   ol['filters'] if ol['filters'] else 'UNKNOWN')
        else:
            pLogs.textWrap.initial_indent = \
            ' %6s %-20s %-10s ' % (ol['run'],
                                   ol['target'] if ol['target'] else 'UNKNOWN',
                                   ol['filters'] if ol['filters'] else 'UNKNOWN')
        pLogs.textWrap.subsequent_indent = 40*' '
        if ol['comment'] == '':
            print pLogs.textWrap.initial_indent
        else:
            print pLogs.textWrap.fill(ol['comment'])
        
if __name__ == '__main__':

    # ok deal with arguments first
    parser = argparse.ArgumentParser(description=descr, 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # positional
    parser.add_argument('night', help='directory with runs to go through, e.g. "2002-05-16"')

    # optional
    parser.add_argument('-f', dest='first', type=int, default=1, help='first run to plot')
    parser.add_argument('-l', dest='local', 
                        help='name of local file to save results rather than google spreadsheet')
    parser.add_argument('-e', dest='email', help='gmail address (will be prompted if needed)')
    parser.add_argument('-p', dest='password', help='gmail password (will be prompted if needed)')
    parser.add_argument('-r', dest='run', help='specific run to plot, e.g. "run003"')
    parser.add_argument('-b', dest='back', action='store_true', help='switch off auto-bias correction')
    parser.add_argument('-c', dest='chop', action='store_true', help='chop final frame of multi-frame runs')
    parser.add_argument('-n', dest='noquery', action='store_true', 
                        help='do not ask whether to show the final frame of multi-frame runs,' + 
                        ' just show it, unless the -c flag is set.')
    parser.add_argument('-w', dest='width', type=int, default=140, 
                        help='character width for log comments wrapping')
    parser.add_argument('-plo', type=float, default=2., help='Lower percentile for intensity display')
    parser.add_argument('-phi', type=float, default=98., help='Upper percentile for intensity display')
    parser.add_argument('-m', dest='max', type=int, default=40, 
                        help='maximum number of frames to plot, spaced evenly over the run')
    parser.add_argument('-t', dest='trim', default='0x0', 
                        help='rows x cols to trim off edges')
    parser.add_argument('-hh', dest='hhelp', action='store_true', 
                        help='print more detailed description of the program and exit')

    # short-circuit if -hh present
    if '-hh' in sys.argv:
        print hhelp
        exit(0)

    # OK, done with arguments.
    args = parser.parse_args()

    # Check arguments
    night   = args.night
    nre     = re.compile(r'(\d\d\d\d)[-_](\d\d)[-_](\d\d)')
    m       = nre.search(night)
    if not m:
        print 'ERROR:',night,'does not have expected YYYY MM DD pattern.'
        exit(1)
    
    year, month, day = m.group(1), m.group(2), m.group(3)
    nstore = year + '.' + month + '.' + day

    if args.trim.find('x') == -1:
        print 'Failed to find an "x" in the -t option argument'
        exit(1)
    else:
        try:
            nrow,ncol = args.trim.split('x')
            nrow = int(nrow)
            ncol = int(ncol)
        except:
            print 'Could not understand -t argument =',args.trim
            print 'Should be of form "1x2", "2x0" etc.'
            exit(1)

    # load log file
    nlog = os.path.join(night, year + '-' + month + '-' + day + '.dat')
    if os.path.isfile(nlog):
        log = ultracam.Log(nlog)
    else:
        olog = os.path.join(night, year + '_' + month + '_' + day + '_log.dat')
        if os.path.isfile(olog):
            log = ultracam.Log(olog)
        else:
            print 'Could not locate log files of the form',nlog,'or',olog
            print 'One of these must be in the run directory for comments'
            print 'and, in old runs, target names and filters as well.'
            exit(1)

    if not os.path.isdir(night):
        print 'ERROR:',night,'is not a directory.'
        exit(1)

    first = args.first
    if first <= 0:
        print 'ERROR: first frame must be > 0'
        exit(1)

    plo = args.plo
    if plo < 0. or plo > 100.:
        print 'ERROR: plo must lie from 0 to 100'
        exit(1)

    phi = args.phi
    if phi < 0. or phi > 100.:
        print 'ERROR: phi must lie from 0 to 100'
        exit(1)

    # generate run list -- a few riders here
    runs = [fname[:-4] for fname in os.listdir(night) if \
                fname.startswith('run') and \
                fname.endswith('.xml') and \
                os.path.isfile(os.path.join(night, fname[:-4] + '.dat'))]

    if len(runs) == 0:
        print 'ERROR: found no runs in directory =',night
        exit(1)
    else:
        print 'Found',len(runs),'runs in directory =',night
    runs.sort()

    # load up targets, filters and comments from the logs
    # these are used as defaults and to provide some context
    # also use this to cut out problem runs
    oneLiners = []
    nruns = []
    for run in runs:
        try:
            rxml  = ultracam.Rhead(os.path.join(night, run + '.xml'))
            if rxml.isPonoff(): continue

            if log.format == 1:
                oneLiners.append({'run' : run, 
                                  'target'  : log.target[run], 
                                  'filters' : log.filters[run],
                                  'comment' : log.comment[run]})
            elif log.format == 2:
                oneLiners.append({'run' : run,
                                  'target'  : rxml.target if rxml.target else None, 
                                  'filters' : rxml.filters if rxml.filters else None, 
                                  'comment' : log.comment[run] if run in log.comment else ''})
            nruns.append(run)
        except ultracam.UltracamError, err:
            print err
            print 'Run',run,'will be skipped.'

    # Transfer over just the OK runs.
    runs = nruns

    # get user name
    user  = getpass.getuser()

    # list of recognised users
    USERS = {'phsaap' : ('trm', 'tom.r.marsh@gmail.com'), 
             'phslax' : ('mcpb', 'mcpbours@gmail.com'), 
             }

    if args.local:
        # store in a local file for later upload
        logFile = args.local if args.local.endswith('.txt') else args.local + '.txt'
        fstore = open(logFile,'a')
    else:

        import gdata.spreadsheet.service

        # write to google spreadsheet
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

        # connect to the google spreadsheet service
        gd_client = gdata.spreadsheet.service.SpreadsheetsService()
        gd_client.email = email
        gd_client.password = password
        gd_client.source = 'rchecker'
        gd_client.ProgrammaticLogin()

        sheets = gd_client.GetSpreadsheetsFeed()
        title  = 'ULTRACAM run IDs'
        for i, entry in enumerate(sheets.entry):
            if entry.title.text == title:
                skey = sheets.entry[i].id.text.split('/')[-1]
                break
        else:
            print 'Failed to find spreadsheet =',title
            exit(1)

        wsheets = gd_client.GetWorksheetsFeed(skey)
        pkey = wsheets.entry[0].id.text.split('/')[-1]
        if not pkey:
            print 'Failed to find page of spreadsheet.'
            exit(1)
    
    # translate user to more meaningful characters if possible
    if user in USERS: user = USERS[user][0]

    # dictionary of recognised data types and short codes for them.
    DTYPES  = {'a' : 'acquisition', 'b' : 'bias', 'd' : 'dark', 
               'f' : 'flat', 's' : 'science', 'j' : 'junk', 
               'p' : 'publicity', 't' : 'technical',
               'q' : 'quit', 'x' : 'standard', 'u' : 'UNSURE',
               'N' : 'next'
               }

    # to identify default bias levels we need an MJD. We do this from the
    # night rather that the data because the latter can be unreliable.
    mjd_bias = ultracam.str2mjd(nstore)

    # stick a blank line in
    print

    # holds defaults
    lastComment   = None
    lastTarget    = None
    lastFilters   = None

    # recognises variable log print requests
    lre = re.compile('l\s+(\d+)')

    # for frame step requests
    sre = re.compile('S\s+([-+]?\d+)')

    # counter to step through the runs and potentially backup
    nr = 0
    while nr < len(runs):
        run = runs[nr]
        runNumber = int(run[3:])
        if (args.run is not None and run != args.run) or runNumber < args.first:
            nr += 1
            continue

        fname  = os.path.join(night, run)
        reply = raw_input('<cr> to display ' + fname + ', q(uit): ')
        if reply == 'q': break

        rdat   = ultracam.Rdata(fname)
        nframe = rdat.ntotal()
        nskip  = nframe // args.max + 1

        # default bias levels
        if args.back:
            def_bias = None
            print 'Will not apply default bias levels.'
        else:
            def_bias = ultracam.blevs(mjd_bias, rdat.gainSpeed)
            if def_bias is None:
                print 'Failed to determine the default bias levels for night =',\
                    night,'readout speed =',rdat.gainSpeed
                print 'Please report to trm'
                reply = raw_input('<cr> to continue without bias correction, q(uit)')
                if reply == 'q': break
            else:
                print 'Will apply default bias levels.'

        # Display the exposures
        pgopen('/xs')
        for nf in xrange(1,nframe+1):

            # next bit of code is to try to minimise the number of bytes read for speed
            # read only the full data when needed, and only the times of the frames
            # leading up to a data frame. Note that in the 'else' statement, it
            # can be assumed that a frame has been read
            if (nf-1) % nskip == 0:
                mccd  = rdat(nf)

                if def_bias is not None:
                    # subtract default levels
                    for nc, ccd in enumerate(mccd):
                        for nw, win in enumerate(ccd):
                            win -= def_bias[nc][nw % 2]

                nnext = nf + nskip

                # The final frame is often an abnormally short exposure and
                # it is often better not to show it. If the number of frames 
                # equals the numexp parameter however, it is OK to show.
                if nf < nframe or rdat.numexp == nframe or nframe == 1:
                    plot = True
                elif not args.chop:
                    if args.noquery:
                        plot = True
                    else:
                        reply = raw_input('plot final frame? [n]: ')
                        plot  = reply == 'y' or reply == 'Y'
                else:
                    plot = False

                if plot:
                    # trim CCDs
                    for ccd in mccd:
                        for nw, win in enumerate(ccd):
                            if nw % 2 == 0:
                                win.trim(ncol,0,nrow,0)
                            else:
                                win.trim(0,ncol,nrow,0)

                    # determine plot ranges
                    prange = mccd.plot(plo,phi,close=False)
                    print fname,', frame ' + str(nf) + ' / ' + str(nframe) + \
                        ', time = ' + ultracam.mjd2str(mccd[0].time.mjd) + ', exp = ' + \
                        str(mccd[0].time.expose) + ', ranges:',prange
                else:
                    print 'Final frame not displayed.'

            else:
                if nf + mccd.head.value('Data.ntmin') >= nnext: rdat.time(nf)
                
        # close plot
        pgclos()

        # blank line
        print

        # print out contextual info with the target name, filters from logs
        # also covering a few preceding and following runs.
        pLogs(oneLiners, nr, args.width, 2)

        # blank line
        print 

        # The data type.
        reply = 'Z'
        while reply not in DTYPES: 
            reply = raw_input('a(qui), b(ias), d(ark), f(lat), s(cie), (flu)x,' + \
                                  ' j(unk), p(ubl), t(ech), u(nsure), h(elp),' + \
                                  ' l(ogs), S(tep), N(ext), q(uit): ')
            ml     = lre.match(reply)
            ms     = sre.match(reply)
            if ml:
                pLogs(oneLiners, nr, args.width, int(ml.group(1)))
                print
            elif reply == 'l':
                pLogs(oneLiners, nr, args.width)
                print
            elif ms:
                nr = max(0, min(nr + int(ms.group(1)), len(runs)-1))
                break
            elif reply == 'h':
                print '\nOptions (case sensitive):\n'
                print '  a --- aquisition, possibly moving, rotating'
                print '  b --- bias frame. No extraneous light or readout funnies.'
                print '  d --- dark frame for calibrating dark count rates.'
                print '  f --- twilight sky flat.'
                print '  j --- junk. Data of absolutely no conceivable use.'
                print '  p --- publicity shots.'
                print '  s --- science data. Largely fixed position.'
                print '  t --- technical. Internals flats, noise tests, etc.'
                print '  u --- unsure, i.e. decide later.'
                print '  x --- flux standard.'
                print '  h --- print this help.'
                print '  l --- print full logs. ("l 3" to get +/- 3 only etc.)'
                print '  S --- step the run: e.g. "S -2" to step back 2 runs, "S 0" to repeat.'
                print '  N --- move to next run, skipping the current one.'
                print '  q --- quit the script\n'
        
        # handle some special cases
        if reply == 'q': 
            break
        elif reply == 'N':
            nr += 1
            continue
        elif ms:
            continue

        dtype = DTYPES[reply]

        # The target name
        target = oneLiners[nr]['target'] if oneLiners[nr]['target'] else lastTarget 
        reply = 'h'
        while reply == 'h' or reply == 'H': 
            if target:
                reply = raw_input('Target: h(elp), q(uit) [' + target + ']: ')
            else:
                reply = raw_input('Target: h(elp), q(uit): ')

            if reply == 'h' or reply == 'H':
                print '\nPlease specify the target. <cr> to get the default taken'
                print 'from the logs. "q" will quit the script.\n'

        if reply == 'q' or reply == 'Q':
            break
        elif reply != '':
            target = reply
        lastTarget = target
        
        # The filters
        filters = oneLiners[nr]['filters'] if oneLiners[nr]['filters'] else lastFilters
        reply = 'h'
        while reply == 'h' or reply == 'H':
            if filters:
                reply = raw_input('Filters: h(elp), q(uit) [' + filters + ']: ')
            else:
                reply = raw_input('Filters: h(elp), q(uit): ')

            if reply == 'h' or reply == 'H':
                print '\nPlease specify the filters. <cr> for the default from'
                print 'the logs\n'

        if reply == 'q' or reply == 'Q':
            break
        elif reply != '':
            filters = reply
        lastFilters = filters

        # The comment
        comment = 'h'
        while comment == 'h' or (dtype == 'j' and comment == ''):
            if lastComment is None:
                comment = raw_input('Comment: h(elp), q(uit): ')
            else:
                comment = raw_input('Comment: h(elp), q(uit) [' + lastComment + ']: ')
            if comment == 'h':
                print '\nPlease type a one-liner comment. It can be blank except in the'
                print 'case of junk files where you must say _something_ about what is'
                print 'wrong with the data. The default is your previous comment, if any.'
                print "If you don't want the default, but want to say nothing then"
                print "either type a space or '' followed by <cr>\n"
            elif comment == '':
                if lastComment is not None: comment = lastComment
            elif comment == "''":
                comment = ''
                lastComment = comment
            
            if comment != 'q' and comment != 'Q' and comment != 'h' and comment != '':
                lastComment = comment

        if comment == 'q' or comment == 'Q':
            break

        # Time stamp
        tstamp = time.strftime('%Y.%m.%dT%H.%M.%S', time.gmtime())

        # Store in a tuple to allow easy output in either of the styles
        row = (('timestamp', tstamp), ('night',nstore),\
                   ('run', run), ('datatype', dtype), ('target', target),\
                   ('filters', filters),('checker',user),\
                   ('comment', comment))

        if args.local:
            # write data to a disk file
            fstore.write(' | '.join([elem[1] for elem in row]) + '\n')
            print 'Written data for',nstore+'/'+run,'to',args.local
        else:
            # upload to google spreadsheet
            entry = gd_client.InsertRow(dict(row), skey, pkey) 
            if isinstance(entry, gdata.spreadsheet.SpreadsheetsList):
                print 'Uploaded data for',nstore+'/'+run,'to',title
            else:
                print 'Died while uploading to spreadsheet.'
                exit(1)
        print
        
        # step one run on
        nr += 1
