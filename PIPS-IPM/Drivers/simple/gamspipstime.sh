#!/bin/bash

# set default values
built_dir="pipstmp"
regions="12"
to="0.08"
tbsize="8"
np="1"


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

nblocks=$(echo "8760 * $to / $tbsize + ((8760*${to}) % ${tbsize} > 0) + 1" | bc)
echo "$nblocks"

gams ../simple4pips.gms --NBREGIONS=$regions --TO=$to --RESOLUTION=\(60/60\)  --TBSIZE=$tbsize --METHOD=PIPS subsys=../subsysLinux.txt  --SCENBLOCK=-1 > /dev/null
../../../../build_pips/gmschk -g $GAMSSYSDIR -T -X $nblocks allblocksPips.gdx > /dev/null
time mpirun -np $np ../../../../build_pips/gmspips $nblocks allblocksPips $GAMSSYSDIR  2>&1 | tee pips.out

