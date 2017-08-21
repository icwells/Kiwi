# Kiwi is a series of scripts for managing NCBI Viral RefSeqs in a MySQL database and calling blastx and blastn 
This program is meant specifically for managing an in-lab database, but has been 
uploaded in case the code may prove useful to others.

Copyright 2017 by Shawn Rupp

## Installation
Download the repository:

git clone https://github.com/icwells/Kiwi.git

Most of the scripts are written in python3, but several contain Cython modules which
must be compiled. Cython can be installed from the pypi repository or via Miniconda 
(it is installed by default with the full Anaconda package).

### To install with Miniconda:
conda install cython

### Compiling Kiwi
cd Kiwi/
./install.sh

## Please refer to KiwiReadMe.pdf for more detailed instructions on running the program
