### install MA57

rm *.tar *.tar.gz
mv ma* src
cd src
CWP_TEMP=$(pwd)
./configure FFLAGS=-fPIC --with-metis=../../METIS/src/libmetis.a --prefix=${CWP_TEMP}
make install


