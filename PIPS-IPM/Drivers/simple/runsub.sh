#!/bin/bash
if [ x$1 = x ]; then
   echo "Usage runsub.sh nScenarios"
else
   if [ x$GAMSSYSDIR = x ]; then
      GAMSSYSDIR=$HOME/PIPS/PIPS-IPM/Drivers/gams
   fi
   pushd rank$PMI_RANK
   nScenp1=$(( $1 + 1 ))
   $HOME/PIPS/build/gmspips $nScenp1 block $GAMSSYSDIR
   popd
fi
