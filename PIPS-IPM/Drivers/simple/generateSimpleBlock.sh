#!/bin/bash
mkdir blk$PMI_RANK
SIMPLEPATH=$HOME/PIPS/PIPS-IPM/Drivers/simple
ARGS="--TO=0.0003 --NBREGIONS=5 --nbscen=5 --CPLEXBENDERS=0 " 
ARGS+="subsys=$SIMPLEPATH/subsysLinux.txt fileStem=simple$PMI_RANK "
ARGS+="lo=2 sd=blk$PMI_RANK optdir=blk$PMI_RANK"
echo Generating block $PMI_RANK
gams $SIMPLEPATH/simple $ARGS --METHOD=spexplicitde --SCENBLOCK=$PMI_RANK
if [ $? != 0 ] ; then
    echo Error generating block $PMI_RANK
else
    echo Done generating block $PMI_RANK
fi
rm -rf blk$PMI_RANK
