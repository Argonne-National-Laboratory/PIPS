### install MA27

rm *.tar *.tar.gz
mv ma* src
cd src
CWP_TEMP=$(pwd)
./configure FFLAGS=-fPIC --prefix=${CWP_TEMP}
make install

