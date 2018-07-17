/*
 * GeoStochScaler.C
 *
 *  Created on: 16.07.2018
 *      Author: Svenja Uslu
 */

#define PIPS_DEBUG
#include "GeoStochScaler.h"
#include "StochVector.h"
#include <cmath>
#include "pipsdef.h"

GeoStochScaler::GeoStochScaler(Data* prob, bool bitshifting)
  : QpScaler(prob, bitshifting)
{
   cout<<"Creating GeoStochScaler..."<<endl;
   // todo: make parameters adjustable
   equilibrate = false;
   maxIters = 8;
   minImpr = 0.85;
   goodEnough = 1e3;
}

void GeoStochScaler::doObjScaling()
{
   assert(vec_colscale != NULL);

   obj->componentMult(*vec_colscale);

   const double absmax = obj->infnorm();
   assert(absmax >= 0);
   StochVector* objInvert = dynamic_cast<StochVector*>(obj->clone());
   objInvert->invertSave();
   double absmin = objInvert->infnorm();
   if( absmin != 0.0)
      absmin = 1.0 / absmin;
   const double scaleFactor = std::sqrt(absmax * absmin);

   if( scaleFactor > 0.0 )
      factor_objscale = 1.0 / scaleFactor;
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

void GeoStochScaler::scale()
{
   /* We want to do the direction with lower maximal ratio first,
    * since the absolute smallest value in the scaled matrix is bounded from below by
    * the inverse of the maximum ratio of the direction that is done first */

   StochVector* rowmaxA = dynamic_cast<StochVector*>(bA->clone());
   StochVector* rowminA = dynamic_cast<StochVector*>(bA->clone());
   StochVector* rowmaxC = dynamic_cast<StochVector*>(rhsC->clone());
   StochVector* rowminC = dynamic_cast<StochVector*>(rhsC->clone());
   StochVector* colmax = dynamic_cast<StochVector*>(bux->clone());
   StochVector* colmin = dynamic_cast<StochVector*>(bux->clone());

   const double rowratio = maxRowRatio(*rowmaxA, *rowmaxC, *rowminA, *rowminC);
   const double colratio = maxColRatio(*colmax, *colmin);

   PIPSdebugMessage("rowratio %f \n", rowratio);
   PIPSdebugMessage("colratio %f \n", colratio);

   double p0start, p1start;
   if( colratio < rowratio )
   {
      p0start = colratio;
      p1start = rowratio;
   }
   else
   {
      p0start = rowratio;
      p1start = colratio;
   }

   bool geoscale = p1start > goodEnough;

   if( !geoscale )
   {
      PIPSdebugMessage("No geometric scaling done, ratio good enough.\n");
      // todo: still do equilibrium scaling here or later
      return;
   }

   assert(vec_rowscaleA == NULL && vec_rowscaleC == NULL && vec_colscale == NULL);

   double p0 = 0.0;
   double p1 = 0.0;

   if( geoscale )
   {
      double p0prev = p0start;
      double p1prev = p1start;

      // todo: make a for-loop around that:

      // column scaling first?
      if(colratio < rowratio)
      {
         colmax->componentMult(*colmin);
         colmax->applySqrt();
         vec_colscale = colmax;
         invertAndRound(do_bitshifting, *vec_colscale);

         // todo: unsure here about the value for initialilzeVec
         A->getRowMinMaxVec(false, true, vec_colscale, *rowmaxA);
         C->getRowMinMaxVec(false, true, vec_colscale, *rowmaxC);
         A->getRowMinMaxVec(true, true, vec_colscale, *rowminA);
         C->getRowMinMaxVec(true, true, vec_colscale, *rowminC);

         rowmaxA->componentMult(*rowminA);
         rowmaxC->componentMult(*rowminC);
         rowmaxA->applySqrt();
         rowmaxC->applySqrt();
         vec_rowscaleA = rowmaxA;
         vec_rowscaleC = rowmaxC;

         invertAndRound(do_bitshifting, *vec_rowscaleA);
         invertAndRound(do_bitshifting, *vec_rowscaleC);
      }
      else // row first
      {
         rowmaxA->componentMult(*rowminA);
         rowmaxC->componentMult(*rowminC);
         rowmaxA->applySqrt();
         rowmaxC->applySqrt();
         vec_rowscaleA = rowmaxA;
         vec_rowscaleC = rowmaxC;
         invertAndRound(do_bitshifting, *vec_rowscaleA);
         invertAndRound(do_bitshifting, *vec_rowscaleC);

         A->getColMinMaxVec(false, true, vec_rowscaleA, *colmax);
         C->getColMinMaxVec(false, false, vec_rowscaleC, *colmax);
         A->getColMinMaxVec(true, true, vec_rowscaleA, *colmin);
         C->getColMinMaxVec(true, false, vec_rowscaleC, *colmin);

         colmax->componentMult(*colmin);
         colmax->applySqrt();
         vec_colscale = colmax;

         invertAndRound(do_bitshifting, *vec_colscale);
      }
// todo: loop around and do something with the ratios, only continue if good enough
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
}

GeoStochScaler::~GeoStochScaler()
{
}




