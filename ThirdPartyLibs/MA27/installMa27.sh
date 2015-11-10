### install MA27

tar -x --strip-components 1 -v -f *.tar
rm -f *.tar *.tar.gz
CWP_TEMP=$(pwd)/src
./configure FFLAGS=-fPIC --prefix=${CWP_TEMP}
make install

