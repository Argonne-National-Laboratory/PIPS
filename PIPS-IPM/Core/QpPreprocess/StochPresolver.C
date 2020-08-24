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
#include <limits>
#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctype.h>

#include "StochOptions.h"
#include "PresolveData.h"
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
#include "StochPostsolver.h"
#include "StochPresolverBoundStrengthening.h"
#include "StochPresolverModelCleanup.h"
#include "StochPresolverColumnFixation.h"
#include "StochPresolverSingletonRows.h"
#include "StochPresolverSingletonColumns.h"
#include "StochPresolverParallelRows.h"
#include "pipschecks.h"
#include "pipsport.h"

StochPresolver::StochPresolver(const Data* prob, Postsolver* postsolver = nullptr)
 : QpPresolver(prob, postsolver), my_rank( PIPS_MPIgetRank(MPI_COMM_WORLD) ),
   limit_max_rounds( pips_options::getIntParameter("PRESOLVE_MAX_ROUNDS") ),
   reset_free_variables_after_presolve( pips_options::getBoolParameter("PRESOLVE_RESET_FREE_VARIABLES") ),
   print_problem( pips_options::getBoolParameter("PRESOLVE_PRINT_PROBLEM") ),
   write_presolved_problem( pips_options::getBoolParameter("PRESOLVE_WRITE_PRESOLVED_PROBLEM_MPS") ),
   verbosity( pips_options::getIntParameter("PRESOLVE_VERBOSITY") ),
   presData( new PresolveData(dynamic_cast<const sData*>(origprob), dynamic_cast<StochPostsolver*>(postsolver)) )
{
   const sData* sorigprob = dynamic_cast<const sData*>(origprob);

   if( pips_options::getBoolParameter("PRESOLVE_SINGLETON_ROWS") )
      presolvers.push_back( new StochPresolverSingletonRows(*presData, *sorigprob) );

   if( pips_options::getBoolParameter("PRESOLVE_COLUMN_FIXATION") )
      presolvers.push_back( new StochPresolverColumnFixation(*presData, *sorigprob) );

   if( pips_options::getBoolParameter("PRESOLVE_BOUND_STRENGTHENING") )
      presolvers.push_back( new StochPresolverBoundStrengthening(*presData, *sorigprob) );

   if( pips_options::getBoolParameter("PRESOLVE_PARALLEL_ROWS") )
      presolvers.push_back( new StochPresolverParallelRows(*presData, *sorigprob) );

   if( pips_options::getBoolParameter("PRESOLVE_SINGLETON_COLUMNS") )
      presolvers.push_back( new StochPresolverSingletonColumns(*presData, *sorigprob) );
}

StochPresolver::~StochPresolver()
{
   for( unsigned int i = 0; i < presolvers.size(); ++i )
      delete presolvers[i];
   delete presData;
}

Data* StochPresolver::presolve()
{
   if( my_rank == 0 )
      std::cout << "start stoch presolving" << std::endl;
   presData->printRowColStats();

   const sData* sorigprob = dynamic_cast<const sData*>(origprob);
   assert( sorigprob->isRootNodeInSync() );
   assert( presData->getPresProb().isRootNodeInSync() );

   if( print_problem )
      sorigprob->writeToStreamDense(std::cout);

   /* initialize model clean up (necessary presolver) */
   StochPresolverModelCleanup presolverCleanup(*presData, *sorigprob);

   if( my_rank == 0 && verbosity > 1 )
      std::cout <<"--- Before Presolving: " << std::endl;
   presolverCleanup.countRowsCols();

   // some while iterating over the list over and over until either every presolver says I'm done or some iterlimit is reached?
   presolverCleanup.applyPresolving();

   for( int i = 0; i < limit_max_rounds; ++i )
   {
      bool success = false;
      for( unsigned int i = 0; i < presolvers.size(); ++i )
      {
         bool presolver_success = presolvers[i]->applyPresolving();
         success = success || presolver_success;
      }
   }

   // before the finalize call fix all empty rows and columns not yet fixed
   presolverCleanup.applyPresolving();
   
   if( my_rank == 0 && verbosity > 1 )
      std::cout << "--- After Presolving:" << std::endl;
   presolverCleanup.countRowsCols();
   if( my_rank == 0 )
      std::cout << "Objective offset: " << presData->getObjOffset() << std::endl;
   assert( presData->getPresProb().isRootNodeInSync() );

   if( reset_free_variables_after_presolve )
      resetFreeVariables();

   sData* finalPresData = presData->finalize();
   assert( finalPresData );
   assert( finalPresData->isRootNodeInSync() );

   if( print_problem )
      finalPresData->writeToStreamDense(std::cout);

   if( write_presolved_problem )
   {
      std::ofstream of("presolved.mps");

      if( of.is_open() )
         finalPresData->writeMPSformat(of);
      else
         if( my_rank == 0 )
            std::cout << "Could not open presolved.mps to write out presolved problem!!" << std::endl;
   }

   if( my_rank == 0 )
      std::cout << "end stoch presolving" << std::endl;
   presData->printRowColStats();

   return finalPresData;
}

void StochPresolver::resetFreeVariables()
{
   if( my_rank == 0 )
      std::cout << "Resetting bounds found in bound strengthening" << std::endl;

   const sData* sorigprob = dynamic_cast<const sData*>(origprob);

   presData->resetOriginallyFreeVarsBounds(*sorigprob);
}
