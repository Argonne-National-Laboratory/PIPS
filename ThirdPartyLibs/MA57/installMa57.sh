#!/bin/sh
### install MA57
### Note: if using gcc compilers, running this script may give errors
### like "Reference to unimplemented intrinsic" if g77 is installed.
### Run using F77=gfortran ./installMa57.sh instead.

#assume ma57 in tar.gz file
fn=`ls ma57*.tar.gz`
name=`basename ${fn} .tar.gz`
tar -zxvf $fn
ln -s ./${name} ./src

cd src
./configure FFLAGS=-fPIC --with-metis=../../METIS/src/libmetis.a --prefix=`pwd`
make -j$1 install
