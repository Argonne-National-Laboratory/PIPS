echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

echo "### Downloading ConicBundle:"
if wget https://www-user.tu-chemnitz.de/~helmberg/ConicBundle/CB_v0.3.11.tgz
then
  echo "### ConicBundle: Download Successful.\n"
else
  echo "### ConicBundle: Download Failed.\n"
fi
tar -xzf CB_v0.3.11.tgz 
mv ConicBundle src
rm CB_v0.3.11.tgz 


sed -i  "s/\bECHO.linux   = -e/ECHO.linux   =/g" src/Makefile

