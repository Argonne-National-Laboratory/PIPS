#!/bin/bash

export TCLLIBPATH=${PWD}

tclsh test-functions.tcl
[[ ${?} == 0 ]] || exit 1

echo OK
exit 0
