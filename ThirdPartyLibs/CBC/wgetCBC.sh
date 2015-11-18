echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

echo "### Downloading Cbc:"
if wget http://www.coin-or.org/download/source/Cbc/Cbc-2.9.5.tgz
then
  echo "### Cbc: Download Successful.\n"
else
  echo "### Cbc: Download Failed.\n"
fi
tar -xzf Cbc-2.9.5.tgz
mv Cbc-2.9.5 src
rm Cbc-2.9.5.tgz

cd src
CWP_TEMP=$(pwd)
./configure --enable-static --disable-dependency-tracking --prefix=${CWP_TEMP}

make install
