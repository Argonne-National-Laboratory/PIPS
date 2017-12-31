/*
 * QpScaler.C
 *
 *  Created on: 19.12.2017
 *      Author: Daniel Rehfeldt
 */

#include "QpScaler.h"
#include <algorithm>

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
}

void QpScaler::applyScaling(QpGenData * prob)
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

   // scale obj
   assert(factor_objscale > 0.0);
   obj->componentMult(*vec_colscale);
   obj->scalarMult(factor_objscale);
}

double QpScaler::maxRowRatio(OoqpVector& maxvecA, OoqpVector& maxvecC, OoqpVector& minvecA, OoqpVector& minvecC)
{
   A->getRowMaxVec(maxvecA, NULL, true);
   A->getRowMinVec(minvecA, NULL, true);
   C->getRowMaxVec(maxvecC, NULL, true);
   C->getRowMinVec(minvecC, NULL, true);

   //TODO
   //OoqpVector tmpA = maxvecA;
   //OoqpVector* tmpA ()

   maxvecA.divideSome(minvecA, minvecA);
   maxvecC.divideSome(minvecC, minvecC);

   int index;
   double maxratio;
   maxvecA.max(maxratio, index);
   assert(maxratio >= 0.0);

   double maxvalC;
   maxvecC.max(maxvalC, index);
   assert(maxvalC >= 0.0);

   if( maxvalC > maxratio )
      maxratio = maxvalC;

   return maxratio ;
}

double QpScaler::maxColRatio() //( OoqpVector& colvecmax, OoqpVector& colvecmin )
{
   double maxratio = 0.0;

   return maxratio;
}

QpScaler::~QpScaler()
{
   delete vec_rowscaleQ;
   delete vec_rowscaleA;
   delete vec_rowscaleC;
   delete vec_colscale;
}
