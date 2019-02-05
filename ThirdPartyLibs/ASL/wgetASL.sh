#!/bin/sh

echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

echo "### Downloading Ampl Solver Library (ASL):"
coinasl=1.3.0
if wget -O solvers.tar.gz https://github.com/ampl/mp/archive/$coinasl.tar.gz
then
  echo "### ASL: Download Successful.\n"
else
  echo "### ASL: Download Failed.\n"
  exit 1
fi

echo "Unpacking the source code..."

fn=solvers.tar.gz
tar -xzf $fn
mv mp-$coinasl/src/asl/solvers .
rm -rf mp-$coinasl

name=`basename ${fn} .tar.gz`
ln -s ./${name} ./src

chmod +x src/configure
chmod +x src/configurehere

echo "Applying patch for #define filename in asl.h, which is incompatible with mpi.h."
cp ./patch/asl.h ./src
cp ./patch/dtoa.c ./src

cd src
#./configurehere CC='icc' CFLAGS='-O3 -xMIC-AVX512'
#./configurehere CC='gcc' CFLAGS='-O3'
./configurehere
make -j$1
###############################
