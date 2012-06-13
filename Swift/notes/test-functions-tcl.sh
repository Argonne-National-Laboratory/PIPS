#!/bin/bash

# User may set VALGRIND=/path/to/valgrind

SWIG_DATA=${HOME}/exm/apps/swig-data
export TCLLIBPATH="${SWIG_DATA} ${PWD}"

${VALGRIND} tclsh test-functions.tcl
[[ ${?} == 0 ]] || exit 1

echo OK
exit 0
