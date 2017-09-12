'''This script will cythonize the Kiwi package.'''

from distutils.core import setup
from Cython.Build import cythonize

# Print blank lines to split output
print("\n")
setup(ext_modules=cythonize("blastResults.pyx"))
print("\n")
setup(ext_modules=cythonize("dbIO.pyx"))
print("\n")
setup(ext_modules=cythonize("flatFileClass.pyx"))
