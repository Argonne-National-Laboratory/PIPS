#!/bin/sh

# Test script called by CMake CTest during 'make test'. Can also be used separately.
# Usage:
# (sh/bash) pipsipmLinkConsTest.sh <pipsipmCallbackExample>
#
# The script checks the output of the PIPS-IPM callback example driver for correctness

check_output()
{
  pipsscmd="$1"

  output=$($pipsscmd 2>&1 | grep "$2")
  
  if [ "$output" = "" ]
  then
    return 0 #no match, return false
  else
    return 1
  fi
}

exe=$1

check_output $exe 'solving finished ... objective value: 14'
if [ "$?" -eq "0" ]
then
  echo 'callback example driver not ok'
  exit 1
fi

echo 'linking constraints test passed'
exit 0


