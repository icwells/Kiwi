#!/bin/bash

##############################################################################
# This script will cythonize scripts for the Kiwi package.
# 
# Required programs:	Cython
##############################################################################

# Compile cython scripts and remove build files
cd src/
python setup.py build_ext --inplace
cd ../

mv src/blastResults.*.so bin/blastResults.so
mv src/dbIO.*.so bin/dbIO.so
mv src/flatFileClass.*.so bin/flatFileClass.so

rm -r src/build/
rm src/*.c
