echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

echo "### Downloading Ampl Solver Library (ASL):"
if wget -O solvers.tar http://netlib.sandia.gov/cgi-bin/netlib/netlibfiles.tar?filename=netlib/ampl/solvers
then
  echo "### ASL: Download Successful.\n"
else
  echo "### ASL: Download Failed.\n"
  exit 1 
fi

fn=solvers.tar
name=`basename ${fn} .tar`
tar -xf $fn
ln -s ./${name} ./src

chmod +x src/configure
chmod +x src/configurehere

cd src
./configurehere
make
