/*
 * GeoStochScaler.C
 *
 *  Created on: 16.07.2018
 *      Author: Svenja Uslu
 */

//#define PIPS_DEBUG
#include "GeoStochScaler.h"
#include <cmath>
#include "pipsdef.h"

GeoStochScaler::GeoStochScaler(Data* prob, bool equiScaling, bool bitshifting)
  : QpScaler(prob, bitshifting)
{
   equilibrate = equiScaling;

   // todo: adjust parameters
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

   double absmin = 0.0;
   obj->absminNonZero( absmin, pips_eps );

   if( absmin < pips_eps )
   {
      PIPSdebugMessage("Objective Scaling: absmin=%f < %f, setting objective to zero. \n", absmin, pips_eps);
      obj->setToZero();
   }
   else
   {
      const double scaleFactor = std::sqrt(absmax * absmin);
      PIPSdebugMessage("Objective Scaling: absmin=%f, absmax=%f, scaleFactor=%f \n", absmin, absmax, scaleFactor);

      assert( scaleFactor >= 0.0 );
      scaleVector(*obj, scaleFactor);

      if( equilibrate )
      {
         const double absmax = obj->infnorm();
         assert(absmax >= 0);
         scaleVector(*obj, absmax);
      }
   }
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

   const double rowratio = maxRowRatio(*rowmaxA, *rowmaxC, *rowminA, *rowminC, NULL);
   const double colratio = maxColRatio(*colmax, *colmin, NULL, NULL);

   PIPSdebugMessage("rowratio before scaling %f \n", rowratio);
   PIPSdebugMessage("colratio before scaling %f \n", colratio);

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

   assert(vec_rowscaleA == NULL && vec_rowscaleC == NULL && vec_colscale == NULL);

   vec_rowscaleA = rowmaxA;
   vec_rowscaleC = rowmaxC;
   vec_colscale = colmax;
   setScalingVecsToOne();

   bool geoscale = p1start > goodEnough;

   if( !geoscale )
   {
      PIPSdebugMessage("No geometric scaling done, ratio already good enough.\n");
      if( !equilibrate )
         return;
      PIPSdebugMessage("But will still perform equilibrium scaling.\n");
   }

   double p0 = 0.0;
   double p1 = 0.0;

   if( geoscale )
   {
      double p0prev = p0start;
      double p1prev = p1start;

      for( int i = 0; i < maxIters; i++)
      {
         // column scaling first?
         if(colratio < rowratio)
         {
            p0 = maxColRatio(*colmax, *colmin, vec_rowscaleA, vec_rowscaleC);
            applyGeoMean(*colmax, *colmin);
            vec_colscale = colmax;

            invertAndRound(do_bitshifting, *vec_colscale);

            p1 = maxRowRatio(*rowmaxA, *rowmaxC, *rowminA, *rowminC, vec_colscale);
            applyGeoMean(*rowmaxA, *rowminA);
            applyGeoMean(*rowmaxC, *rowminC);
            vec_rowscaleA = rowmaxA;
            vec_rowscaleC = rowmaxC;

            invertAndRound(do_bitshifting, *vec_rowscaleA);
            invertAndRound(do_bitshifting, *vec_rowscaleC);
            PIPSdebugMessage("Geometric Scaling round %d. colratio=%f, rowratio=%f \n", i, p0, p1);
         }
         else // row first
         {
            p0 = maxRowRatio(*rowmaxA, *rowmaxC, *rowminA, *rowminC, vec_colscale);
            applyGeoMean(*rowmaxA, *rowminA);
            applyGeoMean(*rowmaxC, *rowminC);
            vec_rowscaleA = rowmaxA;
            vec_rowscaleC = rowmaxC;
            invertAndRound(do_bitshifting, *vec_rowscaleA);
            invertAndRound(do_bitshifting, *vec_rowscaleC);

            p1 = maxColRatio(*colmax, *colmin, vec_rowscaleA, vec_rowscaleC);
            applyGeoMean(*colmax, *colmin);
            vec_colscale = colmax;

            invertAndRound(do_bitshifting, *vec_colscale);
            PIPSdebugMessage("Geometric Scaling round %d. colratio=%f, rowratio=%f \n", i, p1, p0);
         }
         // if ratio improvement is not good enough, then break:
         PIPSdebugMessage("p0=%f, p0prev=%f, p1=%f, p1prev=%f \n", p0, p0prev, p1, p1prev);
         if( p0 > minImpr * p0prev && p1 > minImpr * p1prev )
            break;

         p0prev = p0;
         p1prev = p1;
      }
      // perform geometric scaling only if there is enough (default 15%) improvement:
      geoscale = (p0 <= minImpr * p0start || p1 <= minImpr * p1start);
   }

   if( !geoscale && !equilibrate )
   {
      PIPSdebugMessage("No geometric scaling done, improvement was not good enough.\n");
   }
   else
   {
      if( equilibrate )
      {
         if( !geoscale )
            setScalingVecsToOne();

         // equiScaling using the scaling vectors from GeoScaling:
         postEquiScale();
         delete rowmaxA;
         delete rowmaxC;
         delete colmax;
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

      const double xrowratio = maxRowRatio(*xrowmaxA, *xrowmaxC, *xrowminA, *xrowminC, NULL);
      const double xcolratio = maxColRatio(*xcolmax, *xcolmin, NULL, NULL);

      PIPSdebugMessage("rowratio after scaling %f \n", xrowratio);
      PIPSdebugMessage("colratio after scaling %f \n", xcolratio);
#endif

      if( equilibrate )
         assert(A->abmaxnorm() <= 2.0 && C->abmaxnorm() <= 2.0);
   }
   delete rowminA;
   delete rowminC;
   delete colmin;
}

void GeoStochScaler::setScalingVecsToOne()
{
   assert( vec_rowscaleA && vec_rowscaleC && vec_colscale );

   vec_rowscaleA->setToConstant(1.0);
   vec_rowscaleC->setToConstant(1.0);
   vec_colscale->setToConstant(1.0);
}

/** apply an approximation to the geometric mean to Vector maxvec:
 * Multiply maxvec and minvec componentwise and take the square root of the result.
 * Return result in maxvec.
 * */
void GeoStochScaler::applyGeoMean(StochVector& maxvec, StochVector& minvec)
{
   assert( maxvec.n == minvec.n );

   maxvec.componentMult(minvec);
   maxvec.applySqrt();
}

/** apply Equilibrium Scaling after having done Geometric Scaling.
 * The scaling vectors vec_rowscaleA, vec_rowscaleC and vec_colscale should contain
 * the previously determined scaling factors.
 */
void GeoStochScaler::postEquiScale()
{
   assert(vec_rowscaleA != NULL && vec_rowscaleC != NULL && vec_colscale != NULL);

   StochVector* rowmaxA = dynamic_cast<StochVector*>(bA->clone());
   StochVector* rowminA = dynamic_cast<StochVector*>(bA->clone());
   StochVector* rowmaxC = dynamic_cast<StochVector*>(rhsC->clone());
   StochVector* rowminC = dynamic_cast<StochVector*>(rhsC->clone());
   StochVector* colmax = dynamic_cast<StochVector*>(bux->clone());
   StochVector* colmin = dynamic_cast<StochVector*>(bux->clone());

   const double rowratio = maxRowRatio(*rowmaxA, *rowmaxC, *rowminA, *rowminC, vec_colscale);
   const double colratio = maxColRatio(*colmax, *colmin, vec_rowscaleA, vec_rowscaleC);

   PIPSdebugMessage("rowratio before Post-EquiScale %f \n", rowratio);
   PIPSdebugMessage("colratio before Post-EquiScale %f \n", colratio);

   // minimum vectors are not needed here
   delete rowminA;
   delete rowminC;
   delete colmin;

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

}

GeoStochScaler::~GeoStochScaler()
{
}




