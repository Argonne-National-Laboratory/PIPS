/*
 * QpScaler.C
 *
 *  Created on: 19.12.2017
 *      Author: Daniel Rehfeldt
 */

//#define PIPS_DEBUG
#include <algorithm>

#include "QpScaler.h"
#include "StochVector.h"
#include "QpGenData.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "pipsdef.h"
#include "pipsport.h"

QpScaler::QpScaler(Data * prob, bool bitshifting)
: Scaler(prob, bitshifting), scaling_applied(false)
{
   QpGenData* qpprob = dynamic_cast<QpGenData*>(prob);

   vec_rowscaleQ = nullptr;
   vec_rowscaleA = nullptr;
   vec_rowscaleC = nullptr;
   vec_colscale = nullptr;

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

double QpScaler::getObjUnscaled(double objval) const
{
   assert(vec_colscale != nullptr);
   assert(factor_objscale > 0.0);

   if( scaling_applied )
      return (objval / factor_objscale);
   else
      return objval;
}

Variables* QpScaler::getVariablesUnscaled(const Variables& vars) const
{
   QpGenVars* qp_vars = new QpGenVars(dynamic_cast<const QpGenVars&>(vars)); 
   unscaleVars(*qp_vars);

   return qp_vars;
};

Residuals* QpScaler::getResidualsUnscaled(const Residuals& resids) const
{
   QpGenResiduals* qp_resids = new QpGenResiduals(dynamic_cast<const QpGenResiduals&>(resids));
   unscaleResids(*qp_resids);

   return qp_resids;
};

void QpScaler::unscaleVars( Variables& vars ) const
{
   if( !scaling_applied )
      return;

   // todo : Q
   assert(problem);
   assert(vec_colscale);
   assert(vec_rowscaleA);
   assert(vec_rowscaleC);

   QpGenVars& qp_vars = dynamic_cast<QpGenVars&>(vars); 

   qp_vars.x->componentMult(*vec_colscale);
   qp_vars.s->componentDiv(*vec_rowscaleC);
   qp_vars.y->componentMult(*vec_rowscaleA);
   qp_vars.z->componentMult(*vec_rowscaleC);

   qp_vars.v->componentMult(*vec_colscale);
   qp_vars.gamma->componentDiv(*vec_colscale);
   qp_vars.w->componentMult(*vec_colscale);
   qp_vars.phi->componentDiv(*vec_colscale);
   qp_vars.t->componentDiv(*vec_rowscaleC);
   qp_vars.lambda->componentMult(*vec_rowscaleC);
   qp_vars.u->componentDiv(*vec_rowscaleC);
   qp_vars.pi->componentMult(*vec_rowscaleC);
}

void QpScaler::unscaleResids( Residuals& resids ) const
{
   if( !scaling_applied )
      return;

   assert(problem);
   assert(vec_colscale);
   assert(vec_rowscaleA);
   assert(vec_rowscaleC);
   
   QpGenResiduals& qp_resids = dynamic_cast<QpGenResiduals&>(resids);

   qp_resids.rQ->componentDiv(*vec_colscale);
   qp_resids.rA->componentDiv(*vec_rowscaleA);
   qp_resids.rC->componentDiv(*vec_rowscaleC);
   qp_resids.rz->componentMult(*vec_rowscaleC);

   if( qp_resids.getNxlow() > 0)
      qp_resids.rv->componentMult(*vec_colscale);      

   if( qp_resids.getNxupp() > 0)
      qp_resids.rw->componentMult(*vec_colscale);

   if( qp_resids.getMclow() > 0 )
      qp_resids.rt->componentDiv(*vec_rowscaleC);

   if( qp_resids.getMcupp() > 0)
      qp_resids.ru->componentDiv(*vec_rowscaleC);
   // nothing to to for rgamma, rphi, rlambda, rpi;

   // gap is scaling resistant

   qp_resids.recomputeResidualNorm();
}

OoqpVector* QpScaler::getPrimalUnscaled(const OoqpVector& solprimal) const
{
   assert(problem && vec_colscale);
   OoqpVector* unscaledprimal = solprimal.cloneFull();

   // unscale primal
   if( scaling_applied )
      unscaledprimal->componentMult(*vec_colscale);

   return unscaledprimal;
}

OoqpVector* QpScaler::getDualEqUnscaled(const OoqpVector& soldual) const
{
   assert(problem && vec_rowscaleA);
   OoqpVector* unscaleddual = soldual.cloneFull();

   // unscale dual
   if( scaling_applied )
      unscaleddual->componentMult(*vec_rowscaleA);

   return unscaleddual;
}

OoqpVector* QpScaler::getDualIneqUnscaled(const OoqpVector& soldual) const
{
   assert(problem && vec_rowscaleC);
   OoqpVector* unscaleddual = soldual.cloneFull();

   // unscale dual
   if( scaling_applied )
      unscaleddual->componentMult(*vec_rowscaleC);

   return unscaleddual;
}

OoqpVector* QpScaler::getDualVarBoundsUppUnscaled(const OoqpVector& soldual) const
{
   assert(problem && vec_colscale);
   OoqpVector* unscaleddual = soldual.cloneFull();

   // unscale primal
   if( scaling_applied )
      unscaleddual->componentDiv(*vec_colscale);

   return unscaleddual;
}

OoqpVector* QpScaler::getDualVarBoundsLowUnscaled(const OoqpVector& soldual) const
{
   assert(problem && vec_colscale);
   OoqpVector* unscaleddual = soldual.cloneFull();

   // unscale primal
   if( scaling_applied )
      unscaleddual->componentDiv(*vec_colscale);

   return unscaleddual;
}


void QpScaler::applyScaling()
{
   // todo scale Q

   doObjScaling();

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

   scaling_applied = true;
}

double QpScaler::maxRowRatio(OoqpVector& maxvecA, OoqpVector& maxvecC, OoqpVector& minvecA, OoqpVector& minvecC, const OoqpVector* colScalevec)
{
   A->getRowMinMaxVec(true, true, colScalevec, minvecA);
   A->getRowMinMaxVec(false, true, colScalevec, maxvecA);
   C->getRowMinMaxVec(true, true, colScalevec, minvecC);
   C->getRowMinMaxVec(false, true, colScalevec, maxvecC);

#ifndef NDEBUG

   if( !colScalevec )
   {
      int j;
      double max;

      maxvecA.max(max, j);
      assert(max < 0 || max == A->abmaxnorm());

      maxvecC.max(max, j);
      assert(max < 0 || max == C->abmaxnorm());
   }
#endif

   if( with_sides )
   {
      bA->absminVecUpdate(minvecA);
      rhsC->absminVecUpdate(minvecC);
      lhsC->absminVecUpdate(minvecC);

      bA->absmaxVecUpdate(maxvecA);
      rhsC->absmaxVecUpdate(maxvecC);
      lhsC->absmaxVecUpdate(maxvecC);
   }

   OoqpVector* const ratiovecA = maxvecA.clone();
   OoqpVector* const ratiovecC = maxvecC.clone();

   ratiovecA->copyFrom(maxvecA);
   ratiovecC->copyFrom(maxvecC);

   ratiovecA->divideSome(minvecA, minvecA);
   ratiovecC->divideSome(minvecC, minvecC);

   int i;
   double maxratio;
   ratiovecA->max(maxratio, i);

   PIPSdebugMessage("max row ratio A: %f \n", maxratio);

   double maxvalC;
   ratiovecC->max(maxvalC, i);

   PIPSdebugMessage("max row ratio C: %f \n", maxvalC);

   if( maxvalC > maxratio )
      maxratio = maxvalC;

   delete ratiovecA;
   delete ratiovecC;

   return maxratio;
}

double QpScaler::maxColRatio(OoqpVector& maxvec, OoqpVector& minvec, const OoqpVector* rowScaleVecA, const OoqpVector* rowScaleVecC)
{
   A->getColMinMaxVec(true, true, rowScaleVecA, minvec);
   C->getColMinMaxVec(true, false, rowScaleVecC, minvec);

   A->getColMinMaxVec(false, true, rowScaleVecA, maxvec);
   C->getColMinMaxVec(false, false, rowScaleVecC, maxvec);

#ifndef NDEBUG
   if( !rowScaleVecA || !rowScaleVecC )
   {
      int j;
      double max;

      maxvec.max(max, j);

      assert(max < 0 || max == std::max(A->abmaxnorm(), C->abmaxnorm()));
   }
#endif

   OoqpVector* const ratiovec = maxvec.clone();

   ratiovec->copyFrom(maxvec);

   ratiovec->divideSome(minvec, minvec);

   int i;
   double maxratio;
   ratiovec->max(maxratio, i);
   assert(maxratio >= 0.0);

   PIPSdebugMessage("max column ratio: %f \n", maxratio);

   delete ratiovec;

   return maxratio;
}

void QpScaler::scaleObjVector(double scaling_factor)
{
   if( scaling_factor > 0.0 )
      factor_objscale = 1.0 / scaling_factor;
   else
      factor_objscale = 1.0;

   if( do_bitshifting )
   {
      int exp;
      const double mantissa = std::frexp(factor_objscale, &exp);

      if( mantissa >= 0.75 )
         factor_objscale = std::ldexp(0.5, exp + 1);
      else
         factor_objscale = std::ldexp(0.5, exp);
   }

   assert(obj);

   if( factor_objscale != 1.0 )
      obj->scalarMult(factor_objscale);

   scaling_applied = true;
}

QpScaler::~QpScaler()
{
   delete vec_rowscaleQ;
   delete vec_rowscaleA;
   delete vec_rowscaleC;
   delete vec_colscale;

}
