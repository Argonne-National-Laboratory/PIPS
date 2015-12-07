fn=SuiteSparse-4.4.5.tar.gz
echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

echo "### Downloading SuiteSparse:"
#if wget http://faculty.cse.tamu.edu/davis/SuiteSparse/${fn}
#then
#  echo "### SuiteSparse: Download Successful.\n"
#else
#  echo "### SuiteSparse: Download Failed.\n"
#  exit 1
#fi

name=`tar -ztf ${fn} |  cut -f1 -d"/" | uniq`
tar -xf ${fn}
ln -s ./${name} ./src

cd src
#need to install openblas-dev pakcage
#use LD_LIBRARY_PATH is not suggested
#LD_LIBRARY_PATH=/usr/local/lib/x86_64-linux-gnu
cd SuiteSparse_config
conf_file=SuiteSparse_config.mk
rm $conf_file
if [[ "$OSTYPE" == "darwin"* ]]; then
	#echo "ln -s ./SuiteSparse_config/SuiteSparse_config_Mac.mk $conf_file" 
	ln -s SuiteSparse_config_Mac.mk $conf_file
elif [[ "$OSTYPE" == "linux"* ]]; then
	ln -s SuiteSparse_config_linux.mk $conf_file
fi
cd ..
make
