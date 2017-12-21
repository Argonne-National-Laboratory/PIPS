/*
 * QpScaler.C
 *
 *  Created on: 19.12.2017
 *      Author: Daniel Rehfeldt
 */

#include "QpScaler.h"

QpScaler::QpScaler(QpGenData * prob)
: Scaler(prob)
{
   vec_rowscaleA = NULL;
   vec_rowscaleC = NULL;
   vec_colscale = NULL;
   factor_objscale = 1.0;

   GenMatrixHandle A = prob->A;
   GenMatrixHandle C = prob->C;
   OoqpVectorHandle obj = prob->g;
   OoqpVectorHandle bA = prob->bA;
   OoqpVectorHandle bux = prob->bux; // upper bound of x
   OoqpVectorHandle blx = prob->blx; // lower bound of x
   OoqpVectorHandle rhsC = prob->bu; // RHS of C
   OoqpVectorHandle lhsC = prob->bl; // LHS of C
}

QpScaler::~QpScaler()
{
   delete vec_rowscaleA;
   delete vec_rowscaleC;
   delete vec_colscale;
}

void QpScaler::applyScaling(QpGenData * prob)
{
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

