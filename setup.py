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
      scripts=['scripts/praw.py','scripts/ustats.py','scripts/pproc.py',
               'scripts/rchecker.py','scripts/ctimes.py','scripts/utimes.py',
               'scripts/badruns.py', 'scripts/checkrun.py'],

      author='Tom Marsh',
      description="Python module for accessing ULTRACAM files",
      author_email='t.r.marsh@warwick.ac.uk',
      url='http://www.astro.warwick.ac.uk/',
      cmdclass = cmdclass
      )
