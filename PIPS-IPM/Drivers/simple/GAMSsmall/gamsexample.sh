#!/bin/bash

# set default values
debug=false
gams_file="eps"
nblocks="3"
np="1"
scale=""
stepLp=""
presolve=""
memcheck="false"

for i in "$@"
do
case $i in
    -GAMSFILE=*|--GAMSFILE=*)
    gams_file_full="${i#*=}"
    shift # past argument=value
    ;;
    -BLOCKS=*|--BLOCKS=*)
    nblocks="${i#*=}"
    shift # past argument=value
    ;;
    -NP=*|--NP=*)
    np="${i#*=}"
    shift # past argument=value
    ;;
    -SCALE=*|--SCALE=*)
    scale="${i#*=}"
    shift # past argument=value
    ;;
    -STEPLP=*|--STEPLP=*)
    stepLp="${i#*=}"
    shift # past argument=value 
    ;;
    -PRESOLVE=*|--PRESOLVE=*)
    presolve="${i#*=}"
    shift # past argument=value
    ;;
    -DEBUG|--DEBUG)
    debug=true
    shift # past argument=value
    ;;
    -MEMCHECK|--MEMCHECK)
    memcheck=true
    shift # past argument=value
    ;;
    *)
          # unknown option
    ;;
esac
done
    
gams_file=$(basename $gams_file_full)
file_path=$(dirname $gams_file_full)
if [ -f $gams_file_full.gms ]; then
   echo "File $FILE exists."
else
   echo "File $gams_file.gms does not exist. Please give a valid .gms file without the .gms extention. Exiting script."
   exit
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd $file_path
if [ -d "run_$gams_file" ]; then
   rm -r run_$gams_file
fi
mkdir run_$gams_file
if [ -f $file_path/PIPSIPMpp.opt ]; then
   cp $file_path/PIPSIPMpp.opt run_$gams_file/
fi

cd run_$gams_file
gams ../$gams_file --METHOD=PIPS > /dev/null

$DIR/../../../../build_pips/gmschk -g $GAMSSYSDIR -T -X $nblocks $gams_file.gdx > /dev/null

if [ "${stepLp,,}" = "true" ]; then
  stepLp="stepLp"
else
  stepLp=""
fi

if [ "${presolve,,}" = "true" ]; then
  presolve="presolve"
else
  presolve=""
fi

if [ "${scale,,}" = "true" ]; then
  scale="scale"
else
  scale=""
fi

if [ "${debug,,}" = "true" ]; then
  debug="gdb --args"
else
  debug=""
fi

if [ "${memcheck,,}" = "true" ]; then
  memcheck="valgrind"
else
  memcheck=""
fi
mpirun -np $np $memcheck $debug $DIR/../../../../build_pips/gmspips $nblocks $gams_file $GAMSSYSDIR $scale $stepLp $presolve 2>&1 | tee pips.out


