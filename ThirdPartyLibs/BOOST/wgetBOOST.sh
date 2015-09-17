echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

echo "### Downloading BOOST:"
if wget http://downloads.sourceforge.net/project/boost/boost/1.58.0/boost_1_58_0.tar.gz
then
  echo "### BOOST: Download Successful.\n"
else
  echo "### BOOST: Download Failed.\n"
fi
tar -xzf boost_1_58_0.tar.gz
mv boost_1_58_0 src
rm boost_1_58_0.tar.gz

