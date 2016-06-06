echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

echo "### Downloading Ampl Solver Library (ASL):"
if wget -O solvers.tar.gz http://www.ampl.com/netlib/ampl/solvers.tgz
then
  echo "### ASL: Download Successful.\n"
else
  echo "### ASL: Download Failed.\n"
  exit 1 
fi

fn=solvers.tar.gz
name=`basename ${fn} .tar.gz`
tar -xzf $fn
ln -s ./${name} ./src

chmod +x src/configure
chmod +x src/configurehere

cd src
./configurehere
make -j4
