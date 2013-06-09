from distutils.core import setup, Extension
from distutils.command.sdist import sdist as _sdist

"""Setup script for Python module for ucm files"""

try:
    from sdist import sdist
    cmdclass = {'sdist': sdist}
except:
    cmdclass = {}

setup(name='trm.ucm',
      version='0.9',
      packages = ['trm', 'trm.ucm'],
      scripts=['scripts/pucm', 'scripts/snorm'],

      author='Tom Marsh',
      description="Python module for reading/writing ucm files",
      author_email='t.r.marsh@warwick.ac.uk',
      url='http://www.astro.warwick.ac.uk/',
      cmdclass = cmdclass
      )

