#!/bin/bash
if [ x$1 = x ]; then
   echo "Usage gensub.sh nScenarios"
else
   mkdir rank$PMI_RANK
   pushd rank$PMI_RANK
   SIMPLEPATH=$HOME/PIPS/PIPS-IPM/Drivers/simple
   ARGS="--TO=0.0003 --NBREGIONS=5 --nbscen=$1 --CPLEXBENDERS=0 " 
   ARGS+="subsys=$SIMPLEPATH/subsysLinux.txt "
   ARGS+="lo=2 --SCENBLOCK=0 --METHOD=spexplicitde"
   rankp1=$(($PMI_RANK + 1))
   nSize=$(($1 / $PMI_SIZE))
   nStart=$((($nSize * PMI_RANK) + 1))
   if [ $rankp1 -eq $PMI_SIZE ] ; then # last process
     nEnd=$1
   else
     nEnd=$(($nSize * (PMI_RANK + 1)))
   fi
   echo Generating blocks $nStart to $nEnd for $PMI_RANK
   gams $SIMPLEPATH/simple $ARGS --SCENLISTSTART=$nStart --SCENLISTEND=$nEnd
   if [ $? != 0 ] ; then
       echo Error generating blocks for $PMI_RANK
   else
       echo Done generating blocks for $PMI_RANK
   fi
   popd
fi
