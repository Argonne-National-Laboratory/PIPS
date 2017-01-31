#!/bin/bash

# This script should help to bootstrap PIPS.

rm -rf ./ThirdPartyLibs/ASL/solvers
rm -rf ./ThirdPartyLibs/ASL/src
rm -rf ./ThirdPartyLibs/ASL/*.tar.gz
rm -rf ./ThirdPartyLibs/CBC/Cbc-2.9.8
rm -rf ./ThirdPartyLibs/CBC/src
rm -rf ./ThirdPartyLibs/CBC/*.tgz
rm -rf ./ThirdPartyLibs/MA57/ma57*
rm -rf ./ThirdPartyLibs/MA57/src
rm -rf ./ThirdPartyLibs/METIS/metis*
rm -rf ./ThirdPartyLibs/METIS/src


cd ./ThirdPartyLibs/ASL
./wgetASL.sh
cd ../..
cd ./ThirdPartyLibs/CBC
./wgetCBC.sh
cd ../..
cd ./ThirdPartyLibs/MA57
cp ~/ma57-3.9.0.tar.gz .
./installMa57.sh
cd ../..
cd ./ThirdPartyLibs/METIS
./wgetMETIS.sh
cd ../..

rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_ALL=OFF -DBUILD_PIPS_NLP=ON -B. -H..
make -j4


