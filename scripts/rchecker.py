#!/usr/bin/env python

descr = \
"""
For ID-ing ULTRACAM runs. Sequentially displays frames of all runs of a
specified night prompting for user input to define a data type, target 
name, the filters and a comment for each run. See option -hh for more.
"""

# builtin imports
import sys, argparse, os, time, urllib, urllib2, re, getpass, datetime

# and some extras
from ppgplot import *
from trm import ultracam

hhelp = \
"""
The program either stores data in a google spreadsheet called "ULTRACAM run
IDs" or in a local file. The latter can be added later manually using the
import mechanism on the google spreadsheet, appending to what is there already
and using "|" as a custom separation character. The latter means you should
not have any "|" characters in your input or everything will become hopelessly
confused.

To use the google spreadsheet upload mechanism will require my sharing it with 
you and the installation of the google data API, "gdata". It also requires you 
to supply your gmail account info which cannot be hard-coded in the script for
obvious reasons. You can either let it prompt you or supply through -e and -p
options on the command-line. Since the script runs through all runs of a night
this is a small overhead.

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
"""

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
    parser.add_argument('-plo', type=float, default=2., help='Lower percentile for intensity display')
    parser.add_argument('-phi', type=float, default=98., help='Upper percentile for intensity display')
    parser.add_argument('-m', dest='max', type=int, default=100, 
                        help='maximum number of frames to plot, spaced evenly over the run')
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
    if m:
        year, month, day = m.group(1), m.group(2), m.group(3)
        nstore = year + '.' + month + '.' + day
    else:
        print 'ERROR:',night,'does not have expected YYYY MM DD pattern.'
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
                os.path.isfile(os.path.join(night, fname[:-4] + '.dat')) and \
                (args.run is None or fname[:-4] == args.run) and \
                int(fname[3:-4]) >= args.first]

    if len(runs) == 0:
        print 'ERROR: found no runs in directory =',night
        exit(1)
    else:
        print 'Found',len(runs),'runs in directory =',night
    runs.sort()

    # get user name
    user  = getpass.getuser()

    # list of recognised users
    USERS = {'phsaap' : ('trm', 'tom.r.marsh@gmail.com')}

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
    DTYPES  = {'a' : 'acquisition', 'b' : 'bias', 'f' : 'flat', 's' : 'science',
               'j' : 'junk', 'p' : 'publicity', 't' : 'technical', \
                   'q' : 'quit', 'e' : 'exit'}

    # ordered dictionary of recognised filter combos and short codes
    FILTERS = ultracam.Odict({'r' : 'ugr', 'i' : 'ugi', 'n' : 'ugNaI', 
                              'nrc' : 'uNBF4170RC1', 'z' : 'ugz'})

    # to identify default bias levels we need an MJD. We do this from the
    # night rather that the data because the latter can be unreliable.
    mjd_bias = ultracam.str2mjd(nstore)

    # stick a blank line in
    print

    # Start plot (stays open)
    pgopen('/xs')

    # the next three hold default values
    lastFilterKey = None
    lastTarget    = None
    lastComment   = None

    for run in runs:
        runNumber = int(run[3:])
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
                if nf < nframe or rdat.numexp == nframe:
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
                    prange = mccd.plot(plo,phi,close=False)
                    print fname,', frame ' + str(nf) + ' / ' + str(nframe) + \
                        ', time = ' + ultracam.mjd2str(mccd[0].time.mjd) + ', plot ranges:',prange
                else:
                    print 'Final frame not displayed.'

            else:
                if nf + mccd.head.value('Data.ntmin') >= nnext: rdat.time(nf)


        print
        # OK now we get user input. First the data type.
        reply = 'x'
        while reply not in DTYPES: 
            reply = raw_input('a(qui), b(ias), f(lat), s(ci),' + \
                                  ' j(unk), p(ub), t(ech), h(elp), q(uit), e(xit): ')
            reply = reply.lower()
            if reply == 'h':
                print '\nData type definitions. Note that you really must'
                print 'be certain if you hit "j" for junk:\n'
                print '  a --- aquisition, possibly moving, rotating'
                print '  b --- bias frame. No extraneous light or readout funnies.'
                print '  f --- twilight sky flat.'
                print '  j --- junk. Data of absolutely no conceivable use.'
                print '  p --- publicity shots.'
                print '  s --- science data. Largely fixed position.'
                print '  t --- technical. Internals flats, noise tests, etc.'
                print '  h --- print this help.'
                print '  q --- quit and advance to next run'
                print '  e --- exit altogether\n'

        if reply == 'q':
            continue
        elif reply == 'e':
            break

        dtype = DTYPES[reply]

        # Then the target name
        reply = 'h'
        while reply == 'h' or reply == 'H': 
            if lastTarget is None:
                reply = raw_input('Target: d(efault), h(elp): ')
            else:
                reply = raw_input('Target: d(efault), h(elp), q(unit) [' + lastTarget + ']: ')
                if reply == '': reply = lastTarget
            if reply == 'h' or reply == 'H':
                print '\nPlease specify the target. Just "d" if it is the same as in the on-line logs'
                print 'The last value (other than "d" or "h" will be saved to make it easier to'
                print 're-enter. "q" will quit this run and move to the next.\n'
            elif reply == 'd' or reply == 'D':
                lastTarget = reply
                reply      = 'DEFAULT'
            elif reply != 'q' and reply != 'Q':
                lastTarget = reply

        if reply == 'q' or reply == 'Q':
            continue
        target = reply
        
        # then the filters
        reply = 'x'
        while reply not in FILTERS and len(reply) < 2:
            if lastFilterKey is None:
                reply = raw_input('Filters: r(ugr), i(ugi), h(elp), q(uit): ')
            else:
                reply = raw_input('Filters: r(ugr), i(ugi), h(elp), q(uit): [' + 
                                  FILTERS[lastFilterKey] + ']: ')
                if reply == '':
                    reply = lastFilterKey

            if reply == 'h':
                print '\nPlease specify all three filters. r and i are shorthands for the common'
                print 'combinations "ugr" and "ugi". Otherwise you can type a free format set of'
                print 'filters or preferably one of the following codes:\n'
                for code, filters in FILTERS.iteritems():
                    print '  .... code =',code,'----->',filters
                print

            elif reply in FILTERS or len(reply) > 1:
                lastFilterKey = reply

        if reply == 'q' or reply == 'Q':
            continue
        if reply in FILTERS:
            filters = FILTERS[reply]
        else:
            filters = reply

        # Finally a comment
        comment = 'h'
        while comment == 'h' or (dtype == 'j' and comment == ''):
            if lastComment is None:
                comment = raw_input('Comment: h(elp), q(uit): ')
            else:
                comment = raw_input('Comment: h(elp), q(uit) [' + lastComment + ']: ')
            if comment == 'h':
                print '\nPlease type a one-liner comment. It can be blank except in the'
                print 'case of junk files where you must say _something_ about what is'
                print 'wrong with the data.\n'
            elif comment == '':
                comment = lastComment
            elif comment == "''":
                comment = ''
                lastComment = comment
            elif comment != 'q' and comment != 'Q':
                lastComment = comment

        if comment == 'q' or comment == 'Q':
            continue

        # take a time stamp
        tstamp = time.strftime('%Y.%m.%dT%H.%M.%S', time.gmtime())

        # OK now have all input. Store in a tuple to allow 
        # easy output in either of the styles
        row = (('timestamp', tstamp), ('night',nstore),\
                   ('run', run), ('datatype', dtype), ('target', target),\
                   ('filters', filters),('checker',user),\
                   ('comment', comment))

        if args.local:
            fstore.write(' | '.join([elem[1] for elem in row]) + '\n')
        else:
            entry = gd_client.InsertRow(dict(row), skey, pkey) 
            if isinstance(entry, gdata.spreadsheet.SpreadsheetsList):
                print 'Uploaded data for',nstore+'/'+run,'to',title
            else:
                print 'Died while uploading to spreadsheet.'
                exit(1)
        print


    pgclos()
