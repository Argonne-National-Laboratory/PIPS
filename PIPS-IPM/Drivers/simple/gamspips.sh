#!/bin/bash

# set default values
built_dir="pipstmp"
regions="12"
to="0.08"
tbsize="8"
np="1"
scale=""
stepLp=""
presolve=""
mins="60"


for i in "$@"
do
case $i in
    -TO=*|--TO=*)
    to="${i#*=}"
    shift # past argument=value
    ;;
    -DIR=*|--DIR=*)
    built_dir="${i#*=}"
    shift # past argument=value
    ;;
    -NBREGIONS=*|--NBREGIONS=*)
    regions="${i#*=}"
    shift # past argument=value
    ;;
    -TBSIZE=*|--TBSIZE=*)
    tbsize="${i#*=}"
    shift # past argument=value
    ;;
    -NP=*|--NP=*)
    np="${i#*=}"
    shift # past argument=value
    ;;
    -MINS=*|--MINS=*)
    mins="${i#*=}"
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
    *)
          # unknown option
    ;;
esac
done

if [ -d "$built_dir" ]; then
  rm -r "$built_dir"
fi

mkdir $built_dir
cd $built_dir


nblocks=$(echo "8760 * $to * (60/$mins) / $tbsize + ((8760*${to}) % ${tbsize} > 0) + 1" | bc)
echo "$nblocks"

gams ../simple4pips.gms --NBREGIONS=$regions --TO=$to --RESOLUTION=\($mins/60\)  --TBSIZE=$tbsize --METHOD=PIPS subsys=../subsysLinux.txt  --SCENBLOCK=-1 > /dev/null
../../../../build_pips/gmschk -g $GAMSSYSDIR -T -X $nblocks allblocksPips.gdx > /dev/null

if [ "$stepLp" = "true" ]; then
  stepLp="stepLp"
else
  stepLp=""
fi

if [ "$presolve" = "true" ]; then
  presolve="presolve"
else
  presolve=""
fi

if [ "$scale" = "true" ]; then
  scale="scale"
elif [ "$scale" = "scaleEqui" ]; then
  scale="scale"
elif [ "$scale" = "scaleGeo" ]; then
  scale="scaleGeo"
elif [ "$scale" = "scaleGeoEqui" ]; then
  scale="scaleGeoEqui"
else
  scale=""
fi

echo "Calling: mpirun -np $np ../../../../build_pips/gmspips $nblocks allblocksPips $GAMSSYSDIR $scale $stepLp $presolve"
mpirun -np $np ../../../../build_pips/gmspips $nblocks allblocksPips $GAMSSYSDIR $scale $stepLp $presolve 2>&1 | tee pips.out


