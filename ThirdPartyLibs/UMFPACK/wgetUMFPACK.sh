fn=SuiteSparse-4.4.5.tar.gz
echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

echo "### Downloading SuiteSparse:"
if wget http://faculty.cse.tamu.edu/davis/SuiteSparse/${fn}
then
  echo "### SuiteSparse: Download Successful.\n"
else
  echo "### SuiteSparse: Download Failed.\n"
  exit 1
fi

name=`tar -ztf ${fn} |  cut -f1 -d"/" | uniq`
tar -xf ${fn}
ln -s ./${name} ./src

cd src
#need to install openblas-dev pakcage
LD_LIBRARY_PATH=/usr/local/lib/x86_64-linux-gnu
make
