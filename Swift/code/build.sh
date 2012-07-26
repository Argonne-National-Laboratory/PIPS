#!/bin/bash

# BUILD.SH
# Build the leaf package: rounding

# This could be a Makefile but I think it is better
# to use bash as a reference example. -Justin

LEAF_PKG=rounding
LEAF_I=${LEAF_PKG}.i
LEAF_SO=libtcl${LEAF_PKG}.so
LEAF_TCL=${LEAF_PKG}.tcl
LEAF_CXX=${LEAF_PKG}.cpp
LEAF_O=${LEAF_PKG}.o
# The SWIG-generated file:
WRAP_CXX=${LEAF_PKG}_wrap.cxx
WRAP_O=${LEAF_PKG}_wrap.o

PIPS_SRC=${HOME}/collab/PIPS
PIPS_SHARED=${PIPS_SRC}/SharedLibraries
PIPS_BUILD=${HOME}/Public/PIPS.build
APP_LIB_DIRS=( ${PIPS_SHARED}/Cbc-2.7.6/lib
               ${PIPS_SHARED}/PARDISO
               ${PIPS_BUILD}/Input
               ${PIPS_BUILD}/Lagrange )
APP_LIB_NAMES=( stochInput OsiClp Osi Clp ClpRecourseSolver
                CoinUtils pardiso412-GNU443-X86-64 )

# Path to swig-data module
SWIG_DATA=/home/wozniak/exm/apps/swig-data
# Path to MPICH
MPI=/home/wozniak/sfw/mpich2-trunk

check()
{
  CODE=${?}
  if [[ ${CODE} != 0 ]]
  then
    MSG=$1
    echo ${MSG}
    exit ${CODE}
  fi
}

TCLSH=$( which tclsh )
check "Could not find tclsh in PATH!"

TCL_HOME=$( cd $( dirname ${TCLSH} )/.. ; /bin/pwd )
check "Could not find Tcl installation!"

echo "using Tcl in ${TCL_HOME}"

TCL_CONFIG=${TCL_HOME}/lib/tclConfig.sh

[[ -f ${TCL_CONFIG} ]]
check "Could not read tclConfig.sh!"

# This loads many Tcl configuration variables
source ${TCL_CONFIG}
check "tclConfig.sh failed!"

[[ -d ${SWIG_DATA} ]]
check "SWIG_DATA=${SWIG_DATA}: not found!"

[[ -d ${MPI} ]]
check "MPI=${MPI}: not found!"

CFLAGS="-fPIC -g -Wall   "
# Include our stuff:
CFLAGS+="-I . "
CFLAGS+="-I ${SWIG_DATA} "
# Include the PIPS stuff:
CFLAGS+="-I ../../Input "
CFLAGS+="-I ../../SharedLibraries/Cbc-2.7.6/include/coin "
CFLAGS+="-I ../../SolverInterface "
CFLAGS+="-I ../../PIPS-S/Basic "
CFLAGS+="-I ../../Lagrange/RecourseSubproblemSolver "
CFLAGS+="-I ${MPI}/include"

# Compile the functions implementation
echo "g++ ${LEAF_CXX} ..."
g++ ${CFLAGS} -c ${LEAF_CXX} -o ${LEAF_O}
check

# Create the Tcl extension
echo "SWIG ${LEAF_I} ..."
swig -includeall -c++ -tcl ${LEAF_I}
check

# TODO: Figure out why this is necessary:
sed -i 's/Rounding_functions_Init/Tclrounding_functions_Init/' ${WRAP_CXX}

# Compile the Tcl extension
echo "g++ ${WRAP_CXX}"
g++ ${CFLAGS} ${TCL_INCLUDE_SPEC} -c ${WRAP_CXX} -o ${WRAP_O}
check

LINK_ARGS=
for D in ${APP_LIB_DIRS[@]}
do
  LINK_ARGS+="-L ${D} "
done
for N in ${APP_LIB_NAMES[@]}
do
  LINK_ARGS+="-l ${N} "
done
for D in ${APP_LIB_DIRS[@]}
do
  LINK_ARGS+="-Wl,-rpath -Wl,${D} "
done

# Cbc CoinUtils uses libz, libbz2
LINK_ARGS+="-l z -l bz2 "
# PARDISO requires blas, lapack
LINK_ARGS+="-l blas -l lapack"

# Build the Tcl extension as a shared library
echo "LINK ${LEAF_SO} ..."
g++ -shared -o ${LEAF_SO} ${LEAF_PKG}_wrap.o ${LEAF_O} ${LINK_ARGS}
check
echo "created library: ${LEAF_SO}"

# Make the Tcl package index
export LEAF_PKG LEAF_SO LEAF_TCL
${TCLSH} make-package.tcl > pkgIndex.tcl
check
echo "created package."

# Tell the user what they need to do to run this
echo "Set in environment: TURBINE_USER_LIB=${PWD}"

exit 0
