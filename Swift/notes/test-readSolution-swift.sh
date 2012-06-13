#!/bin/bash

# Compile and run Swift script that links to functions

SWIG_DATA=${HOME}/exm/apps/swig-data
export TURBINE_USER_LIB="/home/wozniak/collab/PIPS/Swift/notes ${SWIG_DATA}"

stc -u test-readSolution.{swift,tcl}

turbine -l -n 3 test-readSolution.tcl
