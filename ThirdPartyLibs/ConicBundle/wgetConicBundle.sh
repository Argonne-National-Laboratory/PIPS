fn=CB_v0.3.11.tgz
echo " "
echo "##### Downloading the third party packages for PIPS-NLP:"
echo " "

echo "### Downloading ConicBundle:"
if wget https://www-user.tu-chemnitz.de/~helmberg/ConicBundle/${fn}
then
  echo "### ConicBundle: Download Successful.\n"
else
  echo "### ConicBundle: Download Failed.\n"
  exit 1
fi
name=`tar -tf ${fn} | cut -f1 -d"/" | uniq`
tar -xzf ${fn}
ln -s ${name} ./src 

sed -i  "s/\bECHO.linux   = -e/ECHO.linux   =/g" src/Makefile
cd src
make

