### install MA57

#assume ma57 in tar.gz file
fn=`ls ma57*.tar.gz`
name=`basename ${fn} .tar.gz`
tar -zxvf $fn
ln -s ./${name} ./src

cd src
./configure FFLAGS=-fPIC --with-metis=../../METIS/src/libmetis.a --prefix=`pwd`
make install
