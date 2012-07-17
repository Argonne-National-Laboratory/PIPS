#!/bin/bash

export TURBINE_HOME=/home/wozniak/Public/turbine-0.0.4-x86_64
stc rounding-simple.{swift,tcl}

SWIG_DATA=${HOME}/exm/apps/swig-data
export TURBINE_USER_LIB="${PWD} ${SWIG_DATA}"

turbine -l -n 3 rounding-simple.tcl
