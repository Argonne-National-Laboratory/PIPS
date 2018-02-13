/*
 * EquiStochScaler.C
 *
 *  Created on: 20.12.2017
 *      Author: bzfrehfe
 */

//#define PIPS_DEBUG
#include "EquiStochScaler.h"
#include "StochVector.h"
#include <cmath>
#include "pipsdef.h"

EquiStochScaler::EquiStochScaler(Data* prob, bool bitshifting)
  : QpScaler(prob, bitshifting)
{

}

void EquiStochScaler::doObjScaling()
{
   assert(vec_colscale != NULL);

   obj->componentMult(*vec_colscale);

   const double absmax = obj->infnorm();
   assert(absmax >= 0);

   if( absmax > 0.0 )
      factor_objscale = 1.0 / absmax;
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

   obj->scalarMult(factor_objscale);
}

// todo scale Q
void EquiStochScaler::scale()
{

   /* We want to do the direction with lower maximal ratio first,
    * since the absolute smallest value in the scaled matrix is bounded from below by
    * the inverse of the maximum ratio of the direction that is done first */

   StochVectorHandle rowmaxA(dynamic_cast<StochVector*>(bA->clone()));
   StochVectorHandle rowminA(dynamic_cast<StochVector*>(bA->clone()));
   StochVectorHandle rowmaxC(dynamic_cast<StochVector*>(rhsC->clone()));
   StochVectorHandle rowminC(dynamic_cast<StochVector*>(rhsC->clone()));
   StochVectorHandle colmax(dynamic_cast<StochVector*>(bux->clone()));
   StochVectorHandle colmin(dynamic_cast<StochVector*>(bux->clone()));

   const double rowratio = maxRowRatio(*rowmaxA, *rowmaxC, *rowminA, *rowminC);
   const double colratio = maxColRatio(*colmax, *colmin);

   PIPSdebugMessage("rowratio %f \n", rowratio);
   PIPSdebugMessage("colratio %f \n", colratio);

   assert(vec_rowscaleA == NULL && vec_rowscaleC == NULL && vec_colscale == NULL);

   vec_colscale = colmax;
   vec_rowscaleA = rowmaxA;
   vec_rowscaleC = rowmaxC;

   // column scaling first?
   if( colratio < rowratio )
   {
      invertAndRound(do_bitshifting, *vec_colscale);

      A->getRowMinMaxVec(false, true, vec_colscale, *vec_rowscaleA);
      C->getRowMinMaxVec(false, true, vec_colscale, *vec_rowscaleC);

      invertAndRound(do_bitshifting, *vec_rowscaleA);
      invertAndRound(do_bitshifting, *vec_rowscaleC);
   }
   else // row first
   {
      invertAndRound(do_bitshifting, *vec_rowscaleA);
      invertAndRound(do_bitshifting, *vec_rowscaleC);

      A->getColMinMaxVec(false, true, vec_rowscaleA, *vec_colscale);
      C->getColMinMaxVec(false, false, vec_rowscaleC, *vec_colscale);

      invertAndRound(do_bitshifting, *vec_colscale);
   }

   PIPSdebugMessage("before scaling: \n "
         "objnorm: %f \n Anorm:  %f \n Cnorm  %f \n bAnorm %f \n rhsCnorm %f \n lhsCnorm %f \n buxnorm %f \n blxnorm %f \n  ",
        obj->infnorm(), A->abmaxnorm(), C->abmaxnorm(), bA->infnorm(), rhsC->infnorm(), lhsC->infnorm(), bux->infnorm(), blx->infnorm());

   doObjScaling();

   applyScaling();

   PIPSdebugMessage("after scaling: \n "
         "objnorm: %f \n Anorm:  %f \n Cnorm  %f \n bAnorm %f \n rhsCnorm %f \n lhsCnorm %f \n buxnorm %f \n blxnorm %f \n  ",
         obj->infnorm(), A->abmaxnorm(), C->abmaxnorm(), bA->infnorm(), rhsC->infnorm(), lhsC->infnorm(), bux->infnorm(), blx->infnorm());

#ifdef PIPS_DEBUG
   StochVectorHandle xrowmaxA(dynamic_cast<StochVector*>(bA->clone()));
   StochVectorHandle xrowminA(dynamic_cast<StochVector*>(bA->clone()));
   StochVectorHandle xrowmaxC(dynamic_cast<StochVector*>(rhsC->clone()));
   StochVectorHandle xrowminC(dynamic_cast<StochVector*>(rhsC->clone()));
   StochVectorHandle xcolmax(dynamic_cast<StochVector*>(bux->clone()));
   StochVectorHandle xcolmin(dynamic_cast<StochVector*>(bux->clone()));

   const double xrowratio = maxRowRatio(*xrowmaxA, *xrowmaxC, *xrowminA, *xrowminC);
   const double xcolratio = maxColRatio(*xcolmax, *xcolmin);

   PIPSdebugMessage("rowratio after scaling %f \n", xrowratio);
   PIPSdebugMessage("colratio after scaling %f \n", xcolratio);
#endif

   assert(A->abmaxnorm() <= 2.0 && C->abmaxnorm() <= 2.0);
}


EquiStochScaler::~EquiStochScaler()
{
}
