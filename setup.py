from distutils.core import setup, Extension
from distutils.command.sdist import sdist as _sdist

"""Setup script for Python module for ULTRACAM files"""

try:
    from sdist import sdist
    cmdclass = {'sdist': sdist}
except:
    cmdclass = {}

setup(name='trm.ultracam',
      version='0',
      packages = ['trm', 'trm.ultracam'],
      scripts=['scripts/badruns.py', 'scripts/checkrun.py', 'scripts/ctimes.py',
               'scripts/fcam3d2ucm.py', 'scripts/pproc.py', 'scripts/praw.py',
               'scripts/ptimes.py', 'scripts/rchecker.py', 'scripts/tofits.py',
               'scripts/to3dfits.py', 'scripts/utimes.py', 'scripts/ualert.py',
               'scripts/uspchecker.py', 'scripts/uspfix.py', 'scripts/ustats.py',
               'scripts/u2ds9.py', 'scripts/tchecker.py', 'scripts/talert.py',
               'scripts/tnofcorr.py'],

      author='Tom Marsh',
      description="Python module for accessing ULTRACAM files",
      author_email='t.r.marsh@warwick.ac.uk',
      url='http://www.astro.warwick.ac.uk/',
      cmdclass = cmdclass
      )
