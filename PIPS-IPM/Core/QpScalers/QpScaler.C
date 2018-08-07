/*
 * QpScaler.C
 *
 *  Created on: 19.12.2017
 *      Author: Daniel Rehfeldt
 */

#include "QpScaler.h"
#include <algorithm>
#include "pipsdef.h"

#ifndef NBEBUG
#include <algorithm>
#endif

QpScaler::QpScaler(Data * prob, bool bitshifting)
: Scaler(prob, bitshifting)
{
   QpGenData* qpprob = dynamic_cast<QpGenData*>(prob);

   vec_rowscaleQ = NULL;
   vec_rowscaleA = NULL;
   vec_rowscaleC = NULL;
   vec_colscale = NULL;
   factor_objscale = 1.0;

   Q = qpprob->Q;
   A = qpprob->A;
   C = qpprob->C;
   obj = qpprob->g;
   bA = qpprob->bA;
   bux = qpprob->bux; // upper bound of x
   blx = qpprob->blx; // lower bound of x
   rhsC = qpprob->bu; // RHS of C
   lhsC = qpprob->bl; // LHS of C

   factor_objscale = 1.0;
}

double QpScaler::getOrigObj(double objval)
{
   assert(vec_colscale != NULL);
   assert(factor_objscale > 0.0);

   return (objval / factor_objscale);
}

void QpScaler::applyScaling()
{
   // todo scale Q

   // scale A and rhs
   A->ColumnScale(*vec_colscale);
   A->RowScale(*vec_rowscaleA);
   bA->componentMult(*vec_rowscaleA);

   // scale C and lhs, rhs
   C->ColumnScale(*vec_colscale);
   C->RowScale(*vec_rowscaleC);
   rhsC->componentMult(*vec_rowscaleC);
   lhsC->componentMult(*vec_rowscaleC);

   // scale ub and lb of x
   bux->componentDiv(*vec_colscale);
   blx->componentDiv(*vec_colscale);
}

double QpScaler::maxRowRatio(OoqpVector& maxvecA, OoqpVector& maxvecC, OoqpVector& minvecA, OoqpVector& minvecC)
{
   A->getRowMinMaxVec(true, true, NULL, minvecA);
   A->getRowMinMaxVec(false, true, NULL, maxvecA);
   C->getRowMinMaxVec(true, true, NULL, minvecC);
   C->getRowMinMaxVec(false, true, NULL, maxvecC);

#ifndef NDEBUG
   int j;
   double max;

   maxvecA.max(max, j);
   assert(max < 0 || max == A->abmaxnorm());

   maxvecC.max(max, j);
   assert(max < 0 || max == C->abmaxnorm());
#endif

   OoqpVector* const ratiovecA = maxvecA.clone();
   OoqpVector* const ratiovecC = maxvecC.clone();

   ratiovecA->copyFrom(maxvecA);
   ratiovecC->copyFrom(maxvecC);

   ratiovecA->divideSome(minvecA, minvecA);
   ratiovecC->divideSome(minvecC, minvecC);

   int i;
   double maxratio;
   ratiovecA->max(maxratio, i);

   PIPSdebugMessage("max column ratio A: %f", maxratio);

   double maxvalC;
   ratiovecC->max(maxvalC, i);

   PIPSdebugMessage("max column ratio C: %f", maxvalC);

   if( maxvalC > maxratio )
      maxratio = maxvalC;

   delete ratiovecA;
   delete ratiovecC;

   return maxratio;
}

double QpScaler::maxColRatio(OoqpVector& maxvec, OoqpVector& minvec)
{
   A->getColMinMaxVec(true, true, NULL, minvec);
   C->getColMinMaxVec(true, false, NULL, minvec);

   A->getColMinMaxVec(false, true, NULL, maxvec);
   C->getColMinMaxVec(false, false, NULL, maxvec);

#ifndef NDEBUG
   int j;
   double max;

   maxvec.max(max, j);

   assert(max < 0 || max == std::max(A->abmaxnorm(), C->abmaxnorm()));
#endif

   OoqpVector* const ratiovec = maxvec.clone();

   ratiovec->copyFrom(maxvec);

   ratiovec->divideSome(minvec, minvec);

   int i;
   double maxratio;
   ratiovec->max(maxratio, i);
   assert(maxratio >= 0.0);

   PIPSdebugMessage("max column ratio: %f", maxratio);

   delete ratiovec;

   return maxratio;
}

QpScaler::~QpScaler()
{
   delete vec_rowscaleQ;
   delete vec_rowscaleA;
   delete vec_rowscaleC;
   delete vec_colscale;
}
