# What is PIPS-IPM++?

PIPS-IPM++ is a (MPI + OpenMP) parallel interior-point method for doubly bordered block diagonal Linear Programms (QPs are currently not supported). For more information on the current algorithm implemented in PIPS and on how to use it see [here](https://opus4.kobv.de/opus4-zib/files/7432/ip4energy.pdf). Currently, the only general purpose interface to PIPS-IPM++ is via [GAMS](https://www.gams.com/) and thus a GAMS license is required (or you write your own interface).

PIPS-IPM++ is a derivative of the [PIPS](https://github.com/Argonne-National-Laboratory/PIPS) solver originally developed at Argonne National Laboratory.


# INSTALLATION Instructions

Note that at this point we only support installations on Linux systems.

## General installation instructions
1. Install the packages wget, cmake, mpich2, and boost.
You can get them via the following command (xxxx stands for the name of the package):
In Linux(Ubuntu): ```(sudo) apt install xxxx```

Not that PIPS (and its third party libraries) heavily rely on a proper installation of lapack and blas. The correct libraries can be provided to PIPS' cmake in 3 ways (priority as follows):
 * the user can pass the variable MATH_LIBS to cmake passing all necessary linker information for custom defined lapack and blas to used
 * the user can define the environment variable MKLROOT to make PIPS use the Intel MKL libraries (suitable for Intel processing units). Note that this environment variable can and should be set correctly by one of the scrips provided by Intel with each of its installations of MKL (for more information we recommend checking out the Intel MKL installation instructions)
 * the user can not pass/define anything and cmake will try to automatically find blas and lapack on the system (if installed)
Either way, the cmake output shows the finally chosen lapack/blas routines.

2. Check out the PIPS-IPM++ repository from github.

3. Go to the following folders and run the script wgetXXX.sh
ThirdPartyLibs/CBC (deprecated will be removed in the future)
ThirdPartyLibs/METIS
For an example, use command ```sh wgetCBC.sh``` in the folder ThirdPartyLibs/CBC  

4. Download MA27 and/or MA57 from HSL and put the .tar.gz in the correct folder ThirdPartyLibs/MAXX and run the respective install script.
(See ThirdPartyLibs/MA27/README.txt and ThirdPartyLibs/MA57/README.txt for more details.)
(this installation step will become deprecated at some point - in future versions of PIPS it will be sufficient to either have PARDISO or one of the MAXX solvers, currently one always needs an HSL solver even if it is not used)

5. Obtain a PARDISO (best >= 7.2) license and download PARDISO from [here](http://www.pardiso-project.org/).
Got to ThirdPartyLibs and create a folder PARDISO as well as a folder src:
```
cd ThirdPartyLibs
mkdir PARDISO 
mkdir PARDISO/src
```
Copy the correct PARDISO library (the one compiled with GNU) as **pardiso.so** into the folder ThirdPartyLibs/PARDISO/src.
Either copy the lincense (in a file named pardiso.lic) into your home directory or alternatively into the directory pips will be run from or, and this is the recommended way, set the environment variable PARDISO_LIC_PATH to the folder where the pardiso.lic file can be found (see also the [pardiso user guide](https://pardiso-project.org/manual/manual.pdf) ).

We recommend using PARDISO for proper performance of PIPS-IPM++. Currently MA27 and MA57 are not well integrated and tested.


6. You should now be able to compile PIPS-IPM++. Assuming we are trying to install PIPS-IPM++ in the folder PIPSMAINPATH/build, where 
PIPSMAINPATH is the root installation folder, use the following commands in the PIPSMAINPATH
folder to configure and install PIPS:
```{r, engine='bash', withtiming}
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_ALL=FALSE -DBUILD_PIPS_IPM=TRUE
make
```
(alternatively -DCMAKE_BUILD_TYPE=DEBUG)
(alternatively "make -j" to use all your PCs available resources to compile PIPS-IPM++ or "make -jxx" where xx is some max number of threads allowed eg. make -j10)

7. Now the build directory should contain an executable called pipsipmCallbackExample. Run this to run PIPS for the first time and check whether everything is working.

## Running PIPS-IPM++

The actual PIPS-IPM++ executable is given via gmspips. Consult the [best practice guide](https://gitlab.com/beam-me/bpg) Chapter 4.3 in order to learn on how to annotate models in GAMS and how to split them into multiple files so that PIPS-IPM++ can handle them.

1. Assuming you already have a split model in SOMEFOLDER named model0.gdx, model1.gdx, ..., modelN.gdx where N is the number of blocks in your problem.
Then, a typical call to PIPS-IPM++ to solve the model looks like 
```{r, engine='bash', withtiming}
mpirun -np m PIPSMAINPATH/build/gmspips  SOMEFOLDER/model GAMSFOLDER scaleGeo stepLp presolve
```
Note that the actual MPI-Command "mpirun" might differ depending on your system and MPI installation (it could be "srun" etc.). You can always also run PIPS sequentially (leaving out the "mpirun -np n" part) but that would somehow defeat its purpose.
Here GAMSFOLDER is the path your GAMS installations and the ```gams``` executable lie in. The last three arguments are optional but should be provided for best performance.

2. The PIPSIPMpp.opt file is PIPS-IPM++ options file.
There is not yet a detailed description of all adjustable parameters for the options file. Upon start PIPS-IPM++ will look for a file called PIPSIPMpp.opt in the directory where it is run from and read in its options. Each line has the following structure:
```
PARAMETERNAME VALUE TYPE
```
where PARAMETERNAME is the name of the parameter to set, value the value to set it to and type one of "bool/int/double" (depending on whether we are setting a bool, int or double parameter).
 
The interested user can find a load of parameters and some short descriptions in the source files StochOptions.C and QpGenOptions.C but expert knowledge of PIPS-IPM++ is currently required here and in question one should contact the developers (until a concise guide is available).
 
3. The three optional command line arguments turn on certain features. The argument presolve activates PIPS-IPM++ presolving, scaleGeo activates PIPS-IPM++ geometric scaling and stepLp activates the use of primal and dual step lengths in the interior point method. We recommend to activate all of them.

4. The number of threads used by each MPI process in PIPS-IPM++ can be controlled by setting the environment variable OMP_NUM_THREADS. PIPS-IPM++ will complain if this variable is not set. Best performance can usually be achieved by setting
```
export OMP_NUM_THREADS=2
```
but a value of 1 is also fine and generally performance will depend on your instances.

### Profiling and timing for HPC 
PIPS-IPM++ has built-in parallel performance profiling (mostly in the form of detailed timing and extended convergence reporting). To enable this feature, build PIPS with the -DWITH_TIMING option, for example, a typical build command would be
```{r, engine='bash', withtiming}
cmake -DCMAKE_TOOLCHAIN_FILE=../Toolchain.cmake -DWITH_TIMING=ON .. 
```


# LICENSE

PIPS-IPM++ is derivative work of PIPS-IPM.

See LICENSE file.



