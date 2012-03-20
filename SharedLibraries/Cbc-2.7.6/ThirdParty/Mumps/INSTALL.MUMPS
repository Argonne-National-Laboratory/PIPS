Installation of the Mumps linear solver within COIN.

Mumps is a MUltifrontal Massively Parallel sparse direct Solver.

Please visit the home page:

Home page: http://graal.ens-lyon.fr/MUMPS/

and make sure you read the usage conditions at

http://graal.ens-lyon.fr/MUMPS/avail.html

**********************************************************************
Note: It is YOUR RESPONSIBILITY to ensure that you are entitled to
      download and use this third party package.
**********************************************************************

If you are working under Unix/Linux-like system and you have `wget'
installed, you can downloaded the required source code for MUMPS
by executing the included script `get.Mumps' in this directory.

Note, however, that the script might not be looking for the most recent
version of Mumps, and it might fail if it tries to download an outdated
one.  In that case, please check the Mumps home page for the most recent
version number and change the definition of the mumps_ver variable at the
beginning of the get.Mumps script.

If you don't have wget available, follow those instructions:

Download Mumps 4.10:

http://mumps.enseeiht.fr/MUMPS_4.10.0.tar.gz (or try
http://graal.ens-lyon.fr/MUMPS/MUMPS_4.10.0.tar.gz if first URL doesn't work)

and unpack it in this directory:

gunzip MUMPS_4.10.0.tar.gz
tar xf MUMPS_4.10.0.tar

Then, rename the MUMPS_4.10.0 directory just to be MUMPS:

mv MUMPS_4.10.0 MUMPS

**********************************************************************
Additional instructions for building with Msys and MinGW
**********************************************************************

If you are trying to build Mumps in Microsoft Windows using the MSys shell
environment and gnu fortran compilers, you will first need to install the
MinGW gcc 4.2 technology preview in order obtain a Fortran 95 compliant
Fortran compiler. GCC 4.2 is available here:

http://sourceforge.net/project/showfiles.php?group_id=2435&package_id=241304

You need to install the gcc-core, gcc-g++, and gcc-gfortran packages.
Note the executable names in the gcc packages have the suffix "sjlj" and will
not be found by the configuration script of whatever project you are trying to
build, so you must either (1) change their names by deleting the suffix
(beware that if you do this, you may create a conflict with other version of
gcc already installed) or (2) configure with

./configure CC=gcc-sjlj CXX=g++-sjlj F77=gfortran-sjlj

