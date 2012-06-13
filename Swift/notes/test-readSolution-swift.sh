#!/bin/bash

# Compile and run Swift script that links to functions

export TURBINE_USER_LIB=/home/wozniak/collab/PIPS/Swift/notes

stc -u test-readSolution.{swift,tcl}

turbine -l -n 3 test-readSolution.tcl
