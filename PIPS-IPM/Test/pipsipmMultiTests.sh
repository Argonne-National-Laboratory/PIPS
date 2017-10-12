#!/bin/sh

# Test script called by CMake CTest during 'make test'. Can also be used separately.
# Usage:
# pipsipmMultiTestsScript.sh <pipsipmRawDriver> <rootDirRawInput>
#
#
# The script checks the output of the PIPS-IPM driver for correctness when applied to a 
# a couple of small test problems (in 'raw input' format)

check_output()
{
  pipsscmd="$1 $2"
  output=$($pipsscmd 2>&1 | grep "$3")
  
  if [ "$output" = "" ]
  then
    return 0 #no match, return false
  else
    return 1
  fi
}

exe=$1

check_output $exe "$2/20data/problemdata 8" 'optimal objective: 250917'
if [ "$?" -eq "0" ]; then
  echo '20data not ok'
  exit 1
fi

check_output $exe "$2/LandSdata/problemdata 8" 'optimal objective: 224.01'
if [ "$?" -eq "0" ]; then
  echo 'LandSdata not ok'
  exit 1
fi

check_output $exe "$2/ssndata/problemdata 8" 'optimal objective: 0.00'
if [ "$?" -eq "0" ]; then
  echo 'ssndata not ok'
  exit 1
fi

#this is valid but slow...
#check_output $exe "$2/stormdata/problemdata 8" 'optimal objective: 1.5500' #1.55008e+07
#if [ "$?" -eq "0" ]; then
#  #in case the output is formated differently
#  check_output $exe "$2/stormdata/problemdata 8" 'optimal objective: 155007'
#  if [ "$?" -eq "0" ]; then
#    echo 'stormdata not ok'
#    exit 1
#  fi
#fi

echo 'all tests passed'
exit 0


