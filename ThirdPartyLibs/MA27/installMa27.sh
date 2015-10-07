### install MA27

mv ma* src
cd src
CWP_TEMP=$(pwd)
./configure FFLAGS=-fPIC --prefix=${CWP_TEMP}
make install

