#!/bin/sh
fn=Cbc-2.9.8.tgz

echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

echo "### Downloading Cbc:"
if wget http://www.coin-or.org/download/source/Cbc/${fn}
then
  echo "### ${fn}: Download Successful.\n"
else
  echo "### ${fn}: Download Failed.\n"
  exit 1
fi

name=`basename ${fn} .tgz`
tar -zxf $fn
git apply ${name}.patch
ln -s ./${name} ./src

cd src
#./configure --enable-static --prefix=`pwd`
#make -j4 install
#./configure CFLAGS='-O3' CXXFLAGS='-O3' --enable-static --prefix=`pwd`
./configure --enable-static --prefix=`pwd`
make -j$1 install

