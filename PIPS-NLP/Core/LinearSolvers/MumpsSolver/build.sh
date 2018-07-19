#/bin/sh

#example how to build this simple test driver for testing PIPS' MumpsSolver

INCL="-I/home/petra1/work/projects/pips/PIPS/ThirdPartyLibs/MUMPS/include"
LIBS="-L/home/petra1/work/projects/pips/PIPS/ThirdPartyLibs/MUMPS/lib/ -L/home/petra1/work/installs/metis-5.1.0/_installation/lib/ -L/home/petra1/work/installs/scalapack_installer/install/lib"
# -L/export/home/petra1/work/installs/petsc/arch-linux2-c-debug/lib/

LIBS="-L/home/petra1/work/projects/pips/PIPS/ThirdPartyLibs/MUMPS/lib/ -L/home/petra1/work/installs/parmetis-4.0.3/petra_install/lib -L/home/petra1/work/installs/parmetis-4.0.3/metis/lib -L/home/petra1/work/installs/scalapack_installer/install/lib"


rm *.o a.out

mpicxx -c -DWITHOUT_PIPS  $INCL MumpsSolver.C 

mpicxx mumps_test_driver.cpp -O2 -DWITHOUT_PIPS  $INCL $LIBS  MumpsSolver.o -ldmumps -lgfortran -lblas -llapack -lmumps_common -lparmetis -lmetis -lmpi  -fopenmp -lmpifort -lpord -lmumps_common  -lscalapack  -lrefblas  -ltmg -lreflapack 
