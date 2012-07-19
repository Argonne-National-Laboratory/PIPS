#!/bin/bash

export TURBINE_HOME=/home/wozniak/Public/turbine-0.0.4-x86_64

SWIG_DATA=${HOME}/exm/apps/swig-data
export TURBINE_USER_LIB="${PWD} ${SWIG_DATA}"

stc -u rounding-simple.{swift,tcl}
[[ ${?} == 0 ]] || exit 1

turbine -l -n 3 rounding-simple.tcl
