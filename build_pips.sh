#!/bin/bash

## This script should help to bootstrap PIPS.

# exit if a command fails
set -e

mkdir build
cd build
# CMAKE_BUILD_TYPE is RELEASE or DEBUG
# BUILD_ALL=OFF: Do not build all PIPS variants
# BUILD_PIPS_NLP: Build PIPS_NLP
# BUILD_PIPS_DOC: Build the documentation target (make doc)
# DUMP_ON: Dump 1st stage matrix

cmake -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_ALL=OFF -DBUILD_PIPS_NLP=ON -DBUILD_PIPS_DOC=ON -DDUMP=ON -B. -H..

# Build using 4 processes. 
make -j4

