[![Build Status](https://travis-ci.org/michel2323/PIPS.svg?branch=develop)](https://travis-ci.org/michel2323/PIPS)

# What is PIPS?

PIPS is a suite of parallel optimization solvers mainly for stochastic optimization problems consisting of the following solvers:
 * PIPS-IPM - parallel MPI+OpenMP interior-point for stochastic LPs and convex QPs
 * PIPS-S   - parallel MPI implementation of the revised simplex method
 * PIPS-NLP - parallel MPI interior-point for structured NLPs

# INSTALLATION Instructions

## General installation instructions
1. Install package wget, cmake, mpich2, and boost.
You can get them via the following command (xxx stands for the name of the package):
In Linux(Ubuntu): apt-get install xxxx

2. Go to the following folders and run the script wgetXXX.sh
ThirdPartyLibs/ASL  
ThirdPartyLibs/CBC 
ThirdPartyLibs/ConicBundle   
ThirdPartyLibs/METIS
For an example, use command "sh wgetASL.sh" in the folder ThirdPartyLibs/ASL  

3. Download MA27 and MA57 from HSL and put the source code in the correct folder. 
(See ThirdPartyLibs/MA27/README.txt and ThirdPartyLibs/MA57/README.txt for more details.)

4. Assuming we are trying to install PIPS in the folder PIPSMAINPATH/build_pips, where 
PIPSMAINPATH is the root installation folder, use the following commands in the PIPSMAINPATH
folder to configure and install PIPS:
mkdir build_pips
cd build_pips
cmake ..
make

5. The build system will install executables from three directories: PIPS-IPM, PIPS-S and PIPS-NLP. 

## Building PIPS-S only can be achieved via 
1. cmake -DBUILD_ALL=OFF -DBUILD_PIPS_S=ON <path_to_CMakeLists.txt>
or 
2. including "set(BUILD_ALL OFF); set(BUILD_PIPS_S ON)" in a toolchain file, then 
invoking
cmake -DCMAKE_TOOLCHAIN_FILE=<path_to_toolchain_file> <path_to_CMakeLists.txt>

Same applies to PIPS-IPM and PIPS-NLP (option names BUILD_PIPS_IPM and BUILD_PIPS_NLP, 
respectively).

### Note regarding the use of Pardiso with PIPS-IPM
Pardiso solver has a special feature for PIPS to compute certain Schur complements efficiently. To use PIPS with Pardiso you will need to obtain the PARDISO library and a license [here](http://www.pardiso-project.org/). The build system of PIPS will build PARDISO-related binaries automatically if the PARDISO library is placed under ThirdPartyLibs/PARDISO/src/libpardiso.so

### Profiling and timing for HPC 
PIPS-IPM + Pardiso offers best-in-class HPC performance. PIPS has built-in parallel performance profiling (mostly in the form of detailed timing and extended convergence reporting). To enable this feature, build PIPS with the -DWITH_TIMING option, for example, a typical build command would be
```{r, engine='bash', withtiming}
cmake -DCMAKE_TOOLCHAIN_FILE=../Toolchain.cmake -DWITH_TIMING=ON .. 
```

# CONTRIBUTIONS

## PIPS-IPM
Developed by:
  * Cosmin G. Petra - Lawrence Livermore National Laboratory

Contributions from:
  * Miles Lubin - while at Argonne National Laboratory

## PIPS-S

Developed by: 
  * Miles Lubin - while at Argonne National Laboratory 
  * Cosmin G. Petra - Lawrence Livermore National Laboratory
  * Julian Hall - U. of Edinburgh
  
Contributions from:
  * Geoffrey Oxberry - Lawrence Livermore National Laboratory


## PIPS-NLP 

Developed by:
 * Naiyuan Chiang - while at Argonne National Laboratory
 * Cosmin G. Petra - Lawrence Livermore 	 National Laboratory
Contributions from:
 * Victor Zavala - Univ. of Wisconsin-Madison


# LICENSE

See LICENSE file.

PIPS-IPM and PIPS-NLP are derivative works of OOQP (http://pages.cs.wisc.edu/~swright/ooqp/) by E. Michael Gertz and Stephen. Wright

# ACKNOWLEDGMENTS

PIPS has been developed under the financial support of: 
- Department of Energy, Office of Advanced Scientific Computing Research
- Department of Energy, Early Career Program 
- Department of Energy, Office of Electricity Delivery and Energy Reliability



