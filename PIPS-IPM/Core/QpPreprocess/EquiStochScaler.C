/*
 * EquiStochScaler.C
 *
 *  Created on: 20.12.2017
 *      Author: bzfrehfe
 */

//#define PIPS_DEBUG
#include "EquiStochScaler.h"
#include "StochVector.h"
#include "pipsdef.h"
#include "pipsport.h"

#include <cmath>

EquiStochScaler::EquiStochScaler(Data* prob, bool bitshifting)
  : StochScaler(prob, bitshifting)
{
   int myRank = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   if( myRank == 0) std::cout<<"Creating EquiStochScaler..."<<endl;
}

void EquiStochScaler::doObjScaling()
{
   assert(vec_colscale != nullptr);

   obj->componentMult(*vec_colscale);

#if 0 // note: seems to deteriorate performance and stability
   const double absmax = obj->infnorm();
   assert(absmax >= 0);

   scaleObjVector(absmax);
#endif
}

// todo scale Q
void EquiStochScaler::scale()
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

   const double rowratio = maxRowRatio(*rowmaxA, *rowmaxC, *rowminA, *rowminC, nullptr);
   const double colratio = maxColRatio(*colmax, *colmin, nullptr, nullptr);

   PIPSdebugMessage("rowratio %f \n", rowratio);
   PIPSdebugMessage("colratio %f \n", colratio);

   // minimum vectors are not needed here
   delete rowminA;
   delete rowminC;
   delete colmin;

   assert(vec_rowscaleA == nullptr && vec_rowscaleC == nullptr && vec_colscale == nullptr);

   vec_colscale = colmax;
   vec_rowscaleA = rowmaxA;
   vec_rowscaleC = rowmaxC;

   // column scaling first?
   if( colratio < rowratio && !with_sides )
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

   const double xrowratio = maxRowRatio(*xrowmaxA, *xrowmaxC, *xrowminA, *xrowminC, nullptr);
   const double xcolratio = maxColRatio(*xcolmax, *xcolmin, nullptr, nullptr);

   PIPSdebugMessage("rowratio after scaling %f \n", xrowratio);
   PIPSdebugMessage("colratio after scaling %f \n", xcolratio);
#endif

   assert(A->abmaxnorm() <= 2.0 && C->abmaxnorm() <= 2.0);
}


EquiStochScaler::~EquiStochScaler()
{
}
