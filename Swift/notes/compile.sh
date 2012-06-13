#!/bin/bash

# COMPILE.SH
# Compile functions for syntax sanity checking

g++ -c functions-noops.C -o functions-noops.o
[[ ${?} != 0 ]] && exit 1

echo OK
exit 0
