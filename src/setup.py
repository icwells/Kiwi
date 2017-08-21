'''This script will cythonize the Kiwi package.'''

from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize("dbIO.pyx"))
setup(ext_modules=cythonize("flatFileClass.pyx"))
