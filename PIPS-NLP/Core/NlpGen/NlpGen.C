/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "NlpGen.h"
#include "NlpGenData.h"
#include "NlpGenVars.h"
#include "NlpGenResiduals.h"

#include <cassert>

NlpGen::NlpGen( long long nx_, long long my_, long long mz_ )
{
  nx = nx_;
  my = my_;
  mz = mz_;
}


Residuals * NlpGen::makeResiduals( Data * prob_in )
{
  NlpGenData * prob = (NlpGenData *) prob_in;

  return new NlpGenResiduals( la,
			     nx, my, mz,
			     prob->ixlow, prob->ixupp,
			     prob->iclow, prob->icupp );
}


Variables * NlpGen::makeVariables( Data * prob_in )
{
  NlpGenData * prob = (NlpGenData *) prob_in;

  return new NlpGenVars( la,
			nx, my, mz,
			prob->ixlow, prob->ixupp,
			prob->iclow, prob->icupp );
}



/*
Variables * NlpGen::makePrimalVector( Data * prob_in, )
{
  NlpGenData * prob = (NlpGenData *) prob_in;
  v 	= OoqpVectorHandle( la->newVector( nx ) );

  return new NlpGenVars( la,
			nx, my, mz,
			prob->ixlow, prob->ixupp,
			prob->iclow, prob->icupp, nxL, nxU, nsL, nsU );
}
*/


void NlpGen::copyXSYZ_fromArray( OoqpVector& vec_xsyz, double* array_in, const int nb_col)
{
  assert("not implement yet." &&0);
}

void NlpGen::copyXSYZ_toArray( OoqpVector& vec_xsyz, double* array_in, const int nb_col)
{
  assert("not implement yet." &&0);
}

