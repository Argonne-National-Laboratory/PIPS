#!/bin/bash

## This script should help to bootstrap PIPS and build thirdparty libraries.

# Number of threads used to build
NUMTHREADS=4


export CC='gcc'
export CXX='g++'
export CFLAGS='-O3'
export CXXFLAGS='-O3'

# exit if a command fails
set -e

cd ./ThirdPartyLibs/ASL
./wgetASL.sh $NUMTHREADS
cd ../..
cd ./ThirdPartyLibs/CBC
./wgetCBC.sh $NUMTHREADS
cd ../..

cd ./ThirdPartyLibs/MA57
if [ ! -f "ma57-3.9.0.tar.gz" ]; then
  echo "ERROR: Please provide ma57 library from http://www.hsl.rl.ac.uk/catalogue/ma57.html and put it in ThirdPartyLibs/MA57 folder."
  exit 1
fi
./installMa57.sh $NUMTHREADS
cd ../..

cd ./ThirdPartyLibs/MA27
if [ ! -f "ma27-1.0.0.tar.gz" ]; then
  echo "ERROR: Please provide ma27 library from http://www.hsl.rl.ac.uk/catalogue/ma27.html and put it in ThirdPartyLibs/MA27 folder."
  exit 1
fi
./installMa27.sh $NUMTHREADS
cd ../..

cd ./ThirdPartyLibs/METIS
./wgetMETIS.sh $NUMTHREADS
cd ../..

exit 0

