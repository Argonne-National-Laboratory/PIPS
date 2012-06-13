#!/bin/bash

# User may set VALGRIND=/path/to/valgrind

export TCLLIBPATH=${PWD}

${VALGRIND} tclsh test-functions.tcl
[[ ${?} == 0 ]] || exit 1

echo OK
exit 0
