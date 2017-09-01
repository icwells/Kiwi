'''This script will cythonize the Kiwi package.'''

from distutils.core import setup
from Cython.Build import cythonize

# Print blank lines to split output
print()
setup(ext_modules=cythonize("blastResults.pyx"))
print()
setup(ext_modules=cythonize("dbIO.pyx"))
print()
setup(ext_modules=cythonize("flatFileClass.pyx"))
