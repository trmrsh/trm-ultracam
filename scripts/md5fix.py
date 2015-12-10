#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function


usage = \
"""
Fixes the md5sum for any ru.dat files that have changed following a run of
uspfix.py. It hunts out any files with .dat.old files, if there is a file
MD5SUM_YYYY_MM_DD in the same directory but not an equivalent
MD5SUM_YYYY_MM_DD.old file it then computes the md5sum of the .dat file and
re-creates the MD5SUM_YYYY_MM_DD file. i.e. it only computes the md5sums for
files that have changed. It saves a copy of the old md5sums in a file of the
form MD5SUM_YYYY_MM_DD.old which is used in future runs to avoid
re-computation. If re-computation is thought necessary, you need to move this
file back to its original location without the .old appended.

Run this with no arguments. It operates recursively. but because of the
various files saved, it should be very fast to run again as it won't do
much.
"""

import sys, os, re, hashlib, shutil

if len(sys.argv) > 1:
    print(usage)
    exit(1)

def md5sum(filename, blocksize=65536):
    hash = hashlib.md5()
    with open(filename, 'r+b') as f:
        for block in iter(lambda: f.read(blocksize), ''):
            hash.update(block)
    return hash.hexdigest()

if __name__ == '__main__':

    # To match the xml files
    rmat   = re.compile('^run\d\d\d\.dat.old$')
    md     = re.compile('^MD5SUM_\d\d\d\d_\d\d_\d\d$')
    mdo    = re.compile('^MD5SUM_\d\d\d\d_\d\d_\d\d\.old$')

    for rpath, rnames, fnames in os.walk('.'):
        # search for .dat.old files
        runs = [fname[:-4] for fname in fnames if rmat.match(fname)]
        if len(runs):
            mdf  = [os.path.join(rpath,fname) for fname in fnames if md.match(fname)]
            mdof = [os.path.join(rpath,fname) for fname in fnames if mdo.match(fname)]
            if len(mdf) == 0:
                print('Found .dat.old files in',rpath)
                print('but could not MD%SUM_YYYY_MM_DD file.')
            elif len(mdf) > 1:
                print('Found .dat.old files in',rpath)
                print('but there are multiple MD%SUM_YYYY_MM_DD files =',mdf)
            elif len(mdof):
                print('Found .dat.old files in',rpath)
                print('but there is already an MD%SUM_YYYY_MM_DD.old file =',mdof)
            else:
                mdf = mdf[0]

                # copy the old md5sum file
                shutil.copyfile(mdf, mdf + '.old')
                print('\nCopied',mdf,'to',mdf+'.old')

                # make writeable
                os.chmod(mdf,0o644)
                print('Now correcting md5sums\n')
                with open(mdf, 'w') as fout:
                    with open(mdf + '.old', 'r') as fin:
                        for line in fin:
                            cksum, fname = line.split()
                            if fname in runs:
                                # re-run the md5sum on any files with .dat.old equivalents
                                print('Computing md5sum for',fname)
                                cksum = md5sum(os.path.join(rpath, fname))
                            fout.write(cksum + '  ' + fname + '\n')

