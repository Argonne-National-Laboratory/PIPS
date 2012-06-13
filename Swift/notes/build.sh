#!/bin/bash

# Build the leaf package

# This could be a Makefile but I think it is better
# to use bash as a reference example. -Justin

# TODO: Get name for this leaf package from Lubin. -Justin

LEAF_PKG=functions
LEAF_I=${LEAF_PKG}.i
LEAF_SO=libtcl${LEAF_PKG}.so
LEAF_TCL=${LEAF_PKG}.tcl
# Use the noop implementation
IMPL=noops
LEAF_CXX=${LEAF_PKG}-${IMPL}.C
LEAF_O=${LEAF_PKG}-${IMPL}.o
# The SWIG-generated file:
WRAP_CXX=${LEAF_PKG}_wrap.cxx

# Path to swig-data module
SWIG_DATA=${HOME}/exm/apps/swig-data

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

CFLAGS="-fPIC -g -I . -I ${SWIG_DATA} -Wall"

# Compile the functions implementation
g++ ${CFLAGS} -c ${LEAF_CXX} -o ${LEAF_O}
check

# Create the Tcl extension
swig -includeall -c++ -tcl ${LEAF_I}
check

# TODO: Figure out why this is necessary:
sed -i 's/Functions_Init/Tclfunctions_Init/' ${WRAP_CXX}

# Compile the Tcl extension
g++ ${CFLAGS} ${TCL_INCLUDE_SPEC} -c ${WRAP_CXX}
check

# Build the Tcl extension as a shared library
g++ -shared -o ${LEAF_SO} ${LEAF_PKG}_wrap.o ${LEAF_O}
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
