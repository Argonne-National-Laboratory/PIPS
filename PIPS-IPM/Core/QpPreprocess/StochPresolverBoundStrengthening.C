/*
 * StochPresolverBoundStrengthening.C
 *
 *  Created on: 28.05.2018
 *      Author: bzfuslus
 */

#include "StochPresolverBoundStrengthening.h"

StochPresolverBoundStrengthening::StochPresolverBoundStrengthening(PresolveData& presData)
: StochPresolverBase(presData)
{
 // todo
}

StochPresolverBoundStrengthening::~StochPresolverBoundStrengthening()
{
 // todo
}


bool StochPresolverBoundStrengthening::applyPresolving(int& nelims)
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   return true;
}


