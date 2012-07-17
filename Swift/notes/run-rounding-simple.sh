#!/bin/bash

stc rounding-simple.{swift,tcl}

SWIG_DATA=${HOME}/exm/apps/swig-data
export TURBINE_USER_LIB="${PWD} ${SWIG_DATA}"

turbine -l -n 3 rounding-simple.tcl
