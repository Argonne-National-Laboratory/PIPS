### install MA27

#we assume ma27 is in a tar.gz file
tar xzvf *.tar.gz

#copy everything from the newly created directory in the current directory
mv ma27-1.0.0/* .
rm -rf ma27-1.0.0
#configure and build ma27
CWP_TEMP=$(pwd)/src
./configure FFLAGS=-fPIC --prefix=${CWP_TEMP}
make install

