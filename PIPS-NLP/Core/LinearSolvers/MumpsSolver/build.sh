#/bin/sh

rm *.o a.out

#
#example how to build this simple test driver for testing PIPS' MumpsSolver
#

#linux with gcc, netlib ref blas, mpich
#LIBS="-L/home/petra1/work/projects/pips/PIPS/ThirdPartyLibs/MUMPS/lib/ -L/home/petra1/work/installs/metis-5.1.0/_installation/lib/ -L/home/petra1/work/installs/scalapack_installer/install/lib"
#INCL="-I/home/petra1/work/projects/pips/PIPS/ThirdPartyLibs/MUMPS/include"
#mpicxx -c -DWITHOUT_PIPS  $INCL MumpsSolver.C
#mpicc mumps_test_driver.cpp -DWITHOUT_PIPS  $INCL $LIBS  MumpsSolver.o -ldmumps -lgfortran -lblas -llapack -lmumps_common -lparmetis -lmetis -lmpi  -fopenmp -lmpifort -lpord -lmumps_common  -lscalapack  -lrefblas  -ltmg -lreflapack 

#linux cluster with intel mkl
LIBS="-L/g/g92/petra1/MUMPS_5.1.2/lib -L/lib -L/g/g92/petra1/parmetis-4.0.3/petra_install/lib  -L/g/g92/petra1/parmetis-4.0.3/build/Linux-x86_64/libmetis -L/usr/tce/packages/mkl/mkl-2018.0/mkl/lib/intel64"
INCL="-I/g/g92/petra1/MUMPS_5.1.2/include"
mpicxx -c -DWITHOUT_PIPS  $INCL MumpsSolver.C

MKLROOT="/usr/tce/packages/mkl/mkl-2018.0/mkl/lib/intel64"
LAPACK="-L${MKLROOT} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core"
SCALAP="-L${MKLROOT} -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64"

LIBPAR="${SCALAP} ${LAPACK}"

INCSEQ="-I../libseq"
LIBSEQ="${LAPACK} -L../libseq -lmpiseq"

LIBBLAS="-L${MKLROOT} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core"
LIBOTHERS="-lgomp -lpthread -lm -ldl"

mpicxx mumps_test_driver.cpp -m64 -lmpi -lmpifort -DWITHOUT_PIPS  $INCL $LIBS  MumpsSolver.o -qopenmp -ldmumps $LIBBLAS $LIBPAR -lmumps_common -lparmetis -lmetis -lpord -lmumps_common ${LIBOTHERS}  -lifcore -lmpi -lmpifort

#also export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/tce/packages/mkl/mkl-2018.0/lib/
