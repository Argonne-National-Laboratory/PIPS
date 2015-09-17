echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

echo "### Downloading Metis:"
if wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz
then
  echo "### Metis: Download Successful.\n"
else
  echo "### Metis: Download Failed.\n"
fi
tar -xzf metis-4.0.3.tar.gz 
mv metis-4.0.3 src
rm metis-4.0.3.tar.gz

sed -i  "s/\bCOPTIONS =/COPTIONS = -fPIC /g" src/Makefile.in


