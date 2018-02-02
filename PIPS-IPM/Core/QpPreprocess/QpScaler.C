/*
 * QpScaler.C
 *
 *  Created on: 19.12.2017
 *      Author: Daniel Rehfeldt
 */

#include "../QpPreprocess/QpScaler.h"

#include <algorithm>

QpScaler::QpScaler(Data * prob, bool bitshifting)
: Scaler(prob, bitshifting)
{
   const QpGenData* qpprob = dynamic_cast<QpGenData*>(prob);

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

   std::cout << "Anorm before " << A->abmaxnorm() << std::endl;
   std::cout << "Cnorm before " << C->abmaxnorm() << std::endl;

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

   std::cout << "Anorm after " << A->abmaxnorm() << std::endl;
   std::cout << "Cnorm after " << C->abmaxnorm() << std::endl;
}

double QpScaler::maxRowRatio(OoqpVector& maxvecA, OoqpVector& maxvecC, OoqpVector& minvecA, OoqpVector& minvecC)
{
   A->getRowMinMaxVec(true, true, NULL, minvecA);
   A->getRowMinMaxVec(false, true, NULL, maxvecA);
   C->getRowMinMaxVec(true, true, NULL, minvecC);
   C->getRowMinMaxVec(false, true, NULL, maxvecC);

   OoqpVectorHandle ratiovecA(maxvecA.clone());
   OoqpVectorHandle ratiovecC(maxvecC.clone());

   ratiovecA->copyFrom(maxvecA);
   ratiovecC->copyFrom(maxvecC);

   ratiovecA->divideSome(minvecA, minvecA);
   ratiovecC->divideSome(minvecC, minvecC);

   int i;
   double maxratio;
   ratiovecA->max(maxratio, i);
   assert(maxratio >= 0.0);

   double maxvalC;
   ratiovecC->max(maxvalC, i);
   assert(maxvalC >= 0.0);

   if( maxvalC > maxratio )
      maxratio = maxvalC;

   return maxratio;
}

double QpScaler::maxColRatio(OoqpVector& maxvec, OoqpVector& minvec)
{
   A->getColMinMaxVec(true, true, NULL, minvec);
   C->getColMinMaxVec(true, false, NULL, minvec);

   A->getColMinMaxVec(false, true, NULL, maxvec);
   C->getColMinMaxVec(false, false, NULL, maxvec);

   OoqpVectorHandle ratiovec(maxvec.clone());

   ratiovec->copyFrom(maxvec);

   int i;
   double m;
   ratiovec->max(m, i);
   std::cout << "COLmaxvec " << m << std::endl;

   ratiovec->divideSome(minvec, minvec);

   int index;
   double maxratio;
   ratiovec->max(maxratio, index);
   assert(maxratio >= 0.0);

   std::cout << "colmaxratio " << maxratio << std::endl;

   return maxratio;
}

QpScaler::~QpScaler()
{
}
