# 
# Note that you must install METIS before you install MA57
# Please Get MA57 from HSL (http://www.hsl.rl.ac.uk)
# Note that we need double precision FORTRAN source code.
#
# After you have download MA57, please decompress it in the current folder and rename it as src
# Then use the following commands (for ubuntu) within folder src to install MA57

CWP_TEMP=$(pwd)
./configure FFLAGS=-fPIC --with-metis=../../METIS/src/libmetis.a --prefix=${CWP_TEMP}
make install
cp lib/lib* .
