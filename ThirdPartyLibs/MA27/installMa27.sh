#!/bin/sh
### install MA27

#we assume ma27 is in a tar.gz file
fn=`ls ma27*.tar.gz`
name=`basename ${fn} .tar.gz`
tar -zxvf $fn
ln -s ./${name} ./src

#configure and build ma27
cd src
./configure FFLAGS=-fPIC --prefix=`pwd`
make -j4 install

