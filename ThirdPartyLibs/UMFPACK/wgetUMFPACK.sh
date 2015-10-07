echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

echo "### Downloading SuiteSparse:"
if wget http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.5.tar.gz
then
  echo "### SuiteSparse: Download Successful.\n"
else
  echo "### SuiteSparse: Download Failed.\n"
fi
tar -xzf SuiteSparse-4.4.5.tar.gz 
mv SuiteSparse src
rm SuiteSparse-4.4.5.tar.gz

