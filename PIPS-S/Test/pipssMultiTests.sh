#/bin/bash

# Test script called by CMake CTest during 'make test'. Can also be used separately.
# Usage:
# pipssMultiTestsScript.sh <pipssRawDriver> <rootDirRawInput>
#
#
# The script checks the output of the PIPS-S driver for correctness when applied to a 
# a couple of small test problems (in 'raw input' format)

check_output()
{
  pipsscmd="$1 $2 8"
  output=$($pipsscmd 2>&1 | grep "$3")
  
  if [ "$output" == "" ]
  then
    return 0
  else
    return 1
  fi
}

echo "$2/20data/problemdata 8"  > log.log

exe=$1

check_output $exe "$2/20data/problemdata 8" 'Optimal!!! Objective value: 250917.31'
if [ "$?" -eq "0" ]; then
  echo '20data not ok'
  exit 1
fi

check_output $exe "$2/LandSdata/problemdata 8" 'Optimal!!! Objective value: 224.01'
if [ "$?" -eq "0" ]; then
  echo 'LandSdata not ok'
  exit 1
fi

check_output $exe "$2/ssndata/problemdata 8" 'Optimal!!! Objective value: 0.0000'
if [ "$?" -eq "0" ]; then
  echo 'ssndata not ok'
  exit 1
fi

check_output $exe "$2/stormdata/problemdata 8" 'Optimal!!! Objective value: 15500797.5'
if [ "$?" -eq "0" ]; then
  echo 'stormdata not ok'
  exit 1
fi

echo 'all tests passed'
exit 0


