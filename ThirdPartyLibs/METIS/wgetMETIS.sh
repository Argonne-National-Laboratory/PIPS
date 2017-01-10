#!/bin/sh
echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

fn=metis-4.0.3.tar.gz
echo "### Downloading Metis:"
if wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/${fn}
then
  echo "### Metis: Download Successful.\n"
else
  echo "### Metis: Download Failed.\n"
  exit 1
fi

name=`basename ${fn} .tar.gz`
tar -zxvf $fn
ln -s ./${name} ./src

#compile metis
cd src
sed -i  "s/\bCOPTIONS =/COPTIONS = -O3 -g -fPIC /g" Makefile.in
make -j64


