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

# Change names if necessary
mv src/blastResults.*.so src/blastResults.so
mv src/dbIO.*.so src/dbIO.so
mv src/flatFileClass.*.so src/flatFileClass.so

mv src/*.so bin/

rm -r src/build/
rm src/*.c
