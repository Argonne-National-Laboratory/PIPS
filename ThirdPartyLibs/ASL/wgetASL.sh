echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

echo "### Downloading Ampl Solver Library (ASL):"
if wget http://netlib.sandia.gov/cgi-bin/netlib/netlibfiles.tar?filename=netlib/ampl/solvers
then
  echo "### ASL: Download Successful.\n"
  mv netlibfiles.tar* amplASL.tar
else
  echo "### ASL: Download Failed.\n"
fi
tar -xf amplASL.tar 
mv solvers src
rm amplASL.tar

chmod +x src/configure
chmod +x src/configurehere

