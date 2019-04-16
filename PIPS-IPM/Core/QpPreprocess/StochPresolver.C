/*
 * StochPresolver.C
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */

#include "StochPresolver.h"
#include <algorithm>
#include <cassert>
#include <cmath>
//#include <math.h>
#include <iostream>
#include <utility>
#include <vector>

#include "StochVector.h"
#include "StochGenMatrix.h"
#include "SmartPointer.h"
#include "sData.h"
#include "DoubleMatrix.h"
#include "SparseGenMatrix.h"
#include "StochVectorHandle.h"
#include "OoqpVector.h"
#include "StochGenMatrix.h"
#include "sTreeCallbacks.h"
#include "DoubleMatrixTypes.h"
#include "StochPresolverSingletonRows.h"
#include "StochPresolverSingletonColumns.h"
#include "PresolveData.h"
#include "StochPresolverParallelRows.h"
#include "StochPresolverBoundStrengthening.h"
#include "StochPresolverModelCleanup.h"
#include "pipschecks.h"

StochPresolver::StochPresolver(const Data* prob)
 : QpPresolver(prob)
{
}

StochPresolver::~StochPresolver()
{
}

Data* StochPresolver::presolve()
{
   int myRank = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   if(myRank == 0)
   {
      // debugging objective
      std::cout << "============================================================" << std::endl;
      std::cout << "============================================================" << std::endl;
   }
      double absmin, absminnonzero, min, max, infnorm, onenorm, twonorm;
      int min_idx, max_idx;

      const sData* orig_prob = dynamic_cast<const sData*>(origprob);
      StochVector* copy_g = dynamic_cast<const StochVector&>(*(orig_prob->g)).cloneFull();


      /* lower bounds x */
      StochVector* copy = dynamic_cast<const StochVector&>(*(orig_prob->blx)).cloneFull();
      StochVector* copy_idx = dynamic_cast<const StochVector&>(*(orig_prob->ixlow)).cloneFull();
      copy->componentMult(*copy_idx);

      copy->absminNonZero(absminnonzero, 1e-16);
      copy->min(min, min_idx);
      copy->max(max, max_idx);
      infnorm = copy->infnorm();
      onenorm = copy->onenorm();
      twonorm = copy->twonorm();

      if(myRank == 0)
         std::cout << "lower bounds x:" << "\tabsmin: " << absmin << "\tabsminnonzero: " << absminnonzero
            << "\tmin: " << min << "\tmax: " << max << "\tinfnorm: " << infnorm << "\tonenorm: " << onenorm
            << "\ttwonorm: " << twonorm << std::endl;

      copy->componentMult(*copy_g);

      dynamic_cast<StochVector&>(*copy).absmin(absmin);
      (dynamic_cast<StochVector&>(*copy)).absminNonZero(absminnonzero, 1e-16);
      (dynamic_cast<StochVector&>(*copy)).min(min, min_idx);
      (dynamic_cast<StochVector&>(*copy)).max(max, max_idx);
      infnorm = (dynamic_cast<StochVector&>(*copy)).infnorm();
      onenorm = (dynamic_cast<StochVector&>(*copy)).onenorm();
      twonorm = (dynamic_cast<StochVector&>(*copy)).twonorm();
      if(myRank == 0)
      std::cout << "xl*g: " << "\tabsmin: " << absmin << "\tabsminnonzero: " << absminnonzero
            << "\tmin: " << min << "\tmax: " << max << "\tinfnorm: " << infnorm << "\tonenorm: " << onenorm
            << "\ttwonorm: " << twonorm << std::endl;

      /* upper bounds x */
      copy = dynamic_cast<const StochVector&>(*(orig_prob->bux)).cloneFull();
      copy_idx = dynamic_cast<const StochVector&>(*(orig_prob->ixupp)).cloneFull();
      copy->componentMult(*copy_idx);

      copy->absmin(absmin);
      copy->absminNonZero(absminnonzero, 1e-16);
      copy->min(min, min_idx);
      copy->max(max, max_idx);
      infnorm = copy->infnorm();
      onenorm = copy->onenorm();
      twonorm = copy->twonorm();

      if(myRank == 0)
         std::cout << "upper bounds x:" << "\tabsmin: " << absmin << "\tabsminnonzero: " << absminnonzero
            << "\tmin: " << min << "\tmax: " << max << "\tinfnorm: " << infnorm << "\tonenorm: " << onenorm
            << "\ttwonorm: " << twonorm << std::endl;

      copy->componentMult(*copy_g);

      (dynamic_cast<StochVector&>(*copy)).absmin(absmin);
      (dynamic_cast<StochVector&>(*copy)).absminNonZero(absminnonzero, 1e-16);
      (dynamic_cast<StochVector&>(*copy)).min(min, min_idx);
      (dynamic_cast<StochVector&>(*copy)).max(max, max_idx);
      infnorm = (dynamic_cast<StochVector&>(*copy)).infnorm();
      onenorm = (dynamic_cast<StochVector&>(*copy)).onenorm();
      twonorm = (dynamic_cast<StochVector&>(*copy)).twonorm();
      if(myRank == 0)
      std::cout << "xu*g: " << "\tabsmin: " << absmin << "\tabsminnonzero: " << absminnonzero
            << "\tmin: " << min << "\tmax: " << max << "\tinfnorm: " << infnorm << "\tonenorm: " << onenorm
            << "\ttwonorm: " << twonorm << std::endl;

      /* objective */
      copy = dynamic_cast<const StochVector&>(*(orig_prob->g)).cloneFull();

      copy->absmin(absmin);
      copy->absminNonZero(absminnonzero, 1e-16);
      copy->min(min, min_idx);
      copy->max(max, max_idx);
      infnorm = copy->infnorm();
      onenorm = copy->onenorm();
      twonorm = copy->twonorm();
      if(myRank == 0)
      {
      std::cout << "objective: " << "\tabsmin: " << absmin << "\tabsminnonzero: " << absminnonzero
            << "\tmin: " << min << "\tmax: " << max << "\tinfnorm: " << infnorm << "\tonenorm: " << onenorm
            << "\ttwonorm: " << twonorm << std::endl;

      std::cout << "============================================================" << std::endl;
      std::cout << "============================================================" << std::endl;
      }

   if( myRank == 0 )
      std::cout << "start stoch presolving" << std::endl;

   const sData* sorigprob = dynamic_cast<const sData*>(origprob);

#if 0 // todo add flag
   ofstream myfile;
   myfile.open ("before_presolving.txt");
   sorigprob->writeToStreamDense(myfile);
   myfile.close();
#endif

   /* initialize presolve data */
   PresolveData presData(sorigprob);

   assert( rootNodeInSyncSData(*sorigprob) );


   /* initialize all presolvers */
   StochPresolverBoundStrengthening presolverBS(presData);
   StochPresolverParallelRows presolverParallelRow(presData);
   StochPresolverModelCleanup presolverCleanup(presData);
   StochPresolverSingletonRows presolverSR(presData);

#ifndef NDEBUG
   if( myRank == 0 )
      std::cout <<"--- Before Presolving: " << std::endl;
   presolverSR.countRowsCols();
#endif

   // todo loop, and not exhaustive
   // some list holding all presolvers - eg one presolving run
   // some while iterating over the list over and over until either every presolver says im done or some iterlimit is reached?
   for( int i = 0; i < 1; ++i )
   {
      /* singleton rows */
//      presolverSR.applyPresolving();
//      presolverBS.applyPresolving();
//      presolverCleanup.applyPresolving();
//      presolverParallelRow.applyPresolving();
//      presolverSR.applyPresolving();
      presolverCleanup.applyPresolving();
   }


#ifndef NDEBUG
   assert( presolverSR.verifyNnzcounters() );
   if( myRank == 0 )
      std::cout << "--- After Presolving:" << std::endl;
#endif
   presolverCleanup.countRowsCols();

   //assert( rootNodeInSyncSData(*presData.presProb));

   // i assume we actually apply our changes here and then return a valid sData object to the caller
   sData* finalPresData = presData.finalize();

   assert( rootNodeInSyncSData(*finalPresData) );

#if 0
   myfile.open("after_presolving.txt");
   finalPresData->writeToStreamDense(myfile);
   myfile.close();
#endif

   // todo wrong output - n and m in sData are broken
   if( myRank==0 )
   {
      std::cout << "original problem:\t" << sorigprob->nx << " variables\t" << sorigprob->my << " equ. conss\t" << sorigprob->mz << " ineq. conss" << std::endl;
      std::cout << "presolved problem:\t" << finalPresData->nx << " variables\t" << finalPresData->my << " equ. conss\t" << finalPresData->mz << " ineq. conss" << std::endl;
   }
#ifdef TIMING
   std::cout << "sorigprob nx, my, mz" << sorigprob->nx << " " << sorigprob->my << " " << sorigprob->mz << std::endl;
   std::cout << "finalPresData nx, my, mz" << finalPresData->nx << " " << finalPresData->my << " " << finalPresData->mz << std::endl;
#endif



   if(myRank == 0)
   {
      // debugging objective
      std::cout << "============================================================" << std::endl;
      std::cout << "============================================================" << std::endl;
   }

   StochVector* copy_xl = dynamic_cast<const StochVector&>(*(finalPresData->blx)).cloneFull();
   StochVector* copy_xu = dynamic_cast<const StochVector&>(*(finalPresData->bux)).cloneFull();
   StochVector* copy_xl_idx = dynamic_cast<const StochVector&>(*(finalPresData->ixlow)).cloneFull();
   StochVector* copy_xu_idx = dynamic_cast<const StochVector&>(*(finalPresData->ixupp)).cloneFull();


   copy_xl->componentMult(*copy_xl_idx);
   copy_xu->componentMult(*copy_xu_idx);

      /* lower bounds x */
      (dynamic_cast<StochVector&>(*copy_xl)).absmin(absmin);
      (dynamic_cast<StochVector&>(*copy_xl)).absminNonZero(absminnonzero, 1e-16);
      (dynamic_cast<StochVector&>(*copy_xl)).min(min, min_idx);
      (dynamic_cast<StochVector&>(*copy_xl)).max(max, max_idx);
      infnorm = (dynamic_cast<StochVector&>(*copy_xl)).infnorm();
      onenorm = (dynamic_cast<StochVector&>(*copy_xl)).onenorm();
      twonorm = (dynamic_cast<StochVector&>(*copy_xl)).twonorm();
      if(myRank == 0)
      std::cout << "lower bounds x:" << "\tabsmin: " << absmin << "\tabsminnonzero: " << absminnonzero
            << "\tmin: " << min << "\tmax: " << max << "\tinfnorm: " << infnorm << "\tonenorm: " << onenorm
            << "\ttwonorm: " << twonorm << std::endl;

      /* upper bounds x */
      (dynamic_cast<StochVector&>(*copy_xu)).absmin(absmin);
      (dynamic_cast<StochVector&>(*copy_xu)).absminNonZero(absminnonzero, 1e-16);
      (dynamic_cast<StochVector&>(*copy_xu)).min(min, min_idx);
      (dynamic_cast<StochVector&>(*copy_xu)).max(max, max_idx);
      infnorm = (dynamic_cast<StochVector&>(*copy_xu)).infnorm();
      onenorm = (dynamic_cast<StochVector&>(*copy_xu)).onenorm();
      twonorm = (dynamic_cast<StochVector&>(*copy_xu)).twonorm();
      if(myRank == 0)
      std::cout << "upper bounds x:" << "\tabsmin: " << absmin << "\tabsminnonzero: " << absminnonzero
            << "\tmin: " << min << "\tmax: " << max << "\tinfnorm: " << infnorm << "\tonenorm: " << onenorm
            << "\ttwonorm: " << twonorm << std::endl;

      /* objective */
      (dynamic_cast<StochVector&>(*finalPresData->g)).absmin(absmin);
      (dynamic_cast<StochVector&>(*finalPresData->g)).absminNonZero(absminnonzero, 1e-16);
      (dynamic_cast<StochVector&>(*finalPresData->g)).min(min, min_idx);
      (dynamic_cast<StochVector&>(*finalPresData->g)).max(max, max_idx);
      infnorm = (dynamic_cast<StochVector&>(*finalPresData->g)).infnorm();
      onenorm = (dynamic_cast<StochVector&>(*finalPresData->g)).onenorm();
      twonorm = (dynamic_cast<StochVector&>(*finalPresData->g)).twonorm();
      if(myRank == 0)
      std::cout << "objective: " << "\tabsmin: " << absmin << "\tabsminnonzero: " << absminnonzero
            << "\tmin: " << min << "\tmax: " << max << "\tinfnorm: " << infnorm << "\tonenorm: " << onenorm
            << "\ttwonorm: " << twonorm << std::endl;
      if(myRank == 0)
      std::cout << "objective offset: " << presData.getObjOffset() << std::endl;


      copy_g = dynamic_cast<const StochVector&>(*(finalPresData->g)).cloneFull();

      copy_xl->componentMult(*copy_g);
      copy_xu->componentMult(*copy_g);

      dynamic_cast<StochVector&>(*copy_xl).absmin(absmin);
      (dynamic_cast<StochVector&>(*copy_xl)).absminNonZero(absminnonzero, 1e-16);
      (dynamic_cast<StochVector&>(*copy_xl)).min(min, min_idx);
      (dynamic_cast<StochVector&>(*copy_xl)).max(max, max_idx);
      infnorm = (dynamic_cast<StochVector&>(*copy_xl)).infnorm();
      onenorm = (dynamic_cast<StochVector&>(*copy_xl)).onenorm();
      twonorm = (dynamic_cast<StochVector&>(*copy_xl)).twonorm();
      if(myRank == 0)
      std::cout << "xl*g: " << "\tabsmin: " << absmin << "\tabsminnonzero: " << absminnonzero
            << "\tmin: " << min << "\tmax: " << max << "\tinfnorm: " << infnorm << "\tonenorm: " << onenorm
            << "\ttwonorm: " << twonorm << std::endl;

      (dynamic_cast<StochVector&>(*copy_xu)).absmin(absmin);
      (dynamic_cast<StochVector&>(*copy_xu)).absminNonZero(absminnonzero, 1e-16);
      (dynamic_cast<StochVector&>(*copy_xu)).min(min, min_idx);
      (dynamic_cast<StochVector&>(*copy_xu)).max(max, max_idx);
      infnorm = (dynamic_cast<StochVector&>(*copy_xu)).infnorm();
      onenorm = (dynamic_cast<StochVector&>(*copy_xu)).onenorm();
      twonorm = (dynamic_cast<StochVector&>(*copy_xu)).twonorm();
      if(myRank == 0)
      std::cout << "xu*g: " << "\tabsmin: " << absmin << "\tabsminnonzero: " << absminnonzero
            << "\tmin: " << min << "\tmax: " << max << "\tinfnorm: " << infnorm << "\tonenorm: " << onenorm
            << "\ttwonorm: " << twonorm << std::endl;



      if(myRank == 0)
      std::cout << "============================================================" << std::endl;
      if(myRank == 0)
      std::cout << "============================================================" << std::endl;

//   exit(1);
   return finalPresData;
}


