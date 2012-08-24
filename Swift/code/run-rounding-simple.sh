#!/bin/bash

SWIG_DATA=${HOME}/exm/apps/swig-data
export TURBINE_USER_LIB="${PWD} ${SWIG_DATA}"

export ADLB_EXHAUST_TIME=1

stc -u rounding-simple.{swift,tcl}
[[ ${?} == 0 ]] || exit 1

START=$( date +%s )
turbine -l -n 6 rounding-simple.tcl |& tee turbine.out
[[ ${?} == 0 ]] || exit 1
STOP=$( date +%s )
DURATION=$(( ${STOP} - ${START} ))
echo "TURBINE TIME: ${DURATION}"
