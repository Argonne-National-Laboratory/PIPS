/*
 * StochPresolverSingletonRows.C
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

//#define PIPS_DEBUG
#include "StochPresolverSingletonRows.h"
#include <limits>
#include <cmath>

StochPresolverSingletonRows::StochPresolverSingletonRows(PresolveData& presData)
: StochPresolverBase(presData)
{
 // todo
}

StochPresolverSingletonRows::~StochPresolverSingletonRows()
{
 // todo
}


void StochPresolverSingletonRows::applyPresolving()
{

   assert(true); // todo : assert some stuff :P
   assert(presData.reductionsEmpty());

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

#ifndef NDEBUG
   if( myRank == 0 )
      std::cout << "--- Before singleton Row Presolving:" << std::endl;
   countRowsCols();
#endif

   indivObjOffset = 0.0;
   int n_singleton_equality = 0;
   int n_singleton_inequality = 0;

   // todo - i don't think the new bounds thing is a good concept - for B0 all bounds should be the same - for Bl we might get different ones but we
   // just allreduce min and max them..
   //   clearNewBoundsParent(); // rather - assert is empty...

#ifndef NDEBUG
   std::cout << "Initially found " << n_singleton_equality << " equality and " << n_singleton_inequality << " inequality singletons" << std::endl;
#endif

//   n_singleton_equality = initSingletonRows(EQUALITY_SYSTEM);

   int iter = 0;

   // main loop:
   while( n_singleton_equality + n_singleton_inequality > 0 && iter < maxIterSR )
   {
	  countSingletonRows(n_singleton_equality, n_singleton_inequality);

      /* eliminate all singleton rows in equality system */
      if( n_singleton_equality > 0 )
      {
         // main method:
         doSingletonRows(n_singleton_equality, n_singleton_inequality, EQUALITY_SYSTEM);

         synchronizeSum(n_singleton_equality, n_singleton_inequality);
      }
      else if( n_singleton_inequality > 0 )
      {
    	  assert(n_singleton_equality == 0);

         // main method:
         doSingletonRows(n_singleton_inequality, n_singleton_equality, INEQUALITY_SYSTEM);

         synchronizeSum(n_singleton_equality, n_singleton_inequality);
      }


// ???      presData.resetRedCounters();

      int newSREqLink = 0;
      int newSRIneqLink = 0;

      // todo what about inequalities?
      doSingletonLinkRows(newSREqLink, newSRIneqLink);
      synchronizeSum(newSREqLink, newSRIneqLink);

      iter++;
   }

   assert( (n_singleton_equality == 0 && n_singleton_inequality == 0) || iter >= maxIterSR );

   // Sum up individual objOffset and then add it to the global objOffset:
   sumIndivObjOffset();
   presData.addObjOffset(indivObjOffset);

#ifdef TIMING
   if( myRank == 0 )
      std::cout << "Global objOffset is now: " << presData.getObjOffset() << std::endl;
#endif

   if( myRank == 0 )
      std::cout << "Global objOffset is now: " << presData.getObjOffset() << std::endl;

#ifndef NDEBUG
   if( myRank == 0 )
      std::cout << "--- After singleton Row Presolving:" << std::endl;
   countRowsCols();
#endif
}

/** Does one round of singleton rows presolving for system A or C
 *
 * the blocks B,D,Fi,Gi (the blocks Bmat and Blmat of both A and C), the fixation and updating
 * of the columns is done. The fixed variables in one of the Amat blocks are stored in the
 * member variable colAdaptParent. Updating the blocks A,C,F0,G0 using colAdaptParent happens
 * in updateLinkingVarsBlocks() which should be called after this method.
 * Returns the number of newly found singleton rows (equality/inequality system) during adaption of B,D,Fi,Gi.
 */
void StochPresolverSingletonRows::doSingletonRows(int& n_sing_sys,
      int& n_sing_other_sys, SystemType system_type)
{
   n_sing_sys = 0;

   /* find all fixations in system and store them, update bounds whenever singleton for bonds */

   /* processes root node - finds vars to delete and updates bounds */
   procSingletonRowRoot(system_type);

   /* remove singletons from children */
   for( int node = 0; node < nChildren; node++ )
   {
         procSingletonRowChild(node, n_sing_sys, n_sing_other_sys, system_type);
   }

   /* collect deletions from all processes */

   /* apply deletions */


   /* allreduce variable bounds and check for infeasibility */

   /* update lhs and rhs of linking constraints */ // todo probably with reductions?
	updateRhsNRowLink();

   if( system_type == INEQUALITY_SYSTEM )
      combineNewBoundsParent();

   if( !presData.combineColAdaptParent() )
      abortInfeasible(MPI_COMM_WORLD );

   // update the linking variable blocks (A,C,F,G) with the fixations found in doSingletonRowsA:
   updateLinkingVarsBlocks(n_singleton_equality, n_singleton_inequality);

   // update the variable bounds for the linking variables:
   updateLinkingVarsBounds();

}

void StochPresolverSingletonRows::doSingletonLinkRows(int& newSREq, int& newSRIneq)
{
   // todo what if !hasLinking(EQUALITY_SYSTEM), but hasLinking(INEQUALITY_SYSTEM)???
   if( hasLinking(EQUALITY_SYSTEM))
   {
      setCurrentPointersToNull();
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
      presData.resetRedCounters();
      resetEqRhsAdaptionsLink();
      if( hasLinking(INEQUALITY_SYSTEM) )
         resetIneqRhsAdaptionsLink();

      for( size_t it = 0; it < presData.nRowElemsA->children.size(); it++)
      {
         if( !nodeIsDummy( it, EQUALITY_SYSTEM))
         {
            setCPBlmatsChild(presProb->A, (int)it);
            // set pointers ixlow etc and redColChild
            setCPColumnChild((int) it);
            setCPRowLinkEquality();    // set currEqRhs
            currgChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).children[it]->vec);

            std::vector<COLUMNTOADAPT> varsToBeFixed;

            for( int rowIdx = 0; rowIdx < currNnzRow->n; rowIdx++ )
            {
               if( currNnzRow->elements()[rowIdx] == 1.0 && (currBlmat->rowptr[rowIdx].start != currBlmat->rowptr[rowIdx].end) )
               {  // singleton entry is in this child:
                  assert( currBlmat->rowptr[rowIdx].start +1 == currBlmat->rowptr[rowIdx].end );

                  // collect all variables that can be fixed in varsToBeFixed:
                  int colIdx = -1;
                  double aik = 0.0;
                  getValuesForSR(*currBlmat, rowIdx, colIdx, aik);
                  assert( colIdx >= 0 && colIdx < currgChild->n );
                  assert(!PIPSisZero(aik));

                  // todo
                  //if( aik < tolerance3 )
                  // continue;

                  const double val = currEqRhsLink->elements()[rowIdx] / aik;

                  if( (currIxlowChild->elements()[colIdx] != 0.0 && PIPSisLT(val, currxlowChild->elements()[colIdx]) )
                        || (currIxuppChild->elements()[colIdx] != 0.0 && PIPSisLT(currxuppChild->elements()[colIdx], val) ) )
                  {
                     cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<endl;
                     abortInfeasible(MPI_COMM_WORLD);
                  }
                  // check if this variable was already fixed by another row:
                  bool uniqueFixVariable = true;
                  for(int i=0; i<(int)varsToBeFixed.size(); i++)
                  {
                     if( varsToBeFixed[i].colIdx == colIdx )
                     {
                        uniqueFixVariable = false;
                        if( !PIPSisEQ(varsToBeFixed[i].val, val) )
                           abortInfeasible(MPI_COMM_WORLD);
                     }
                  }
                  if( uniqueFixVariable )
                  {
                     // adapting objOffset:
                     indivObjOffset += currgChild->elements()[colIdx] * val;

                     // store in colADaptLinkBlock for G, F, B, D:
                     COLUMNTOADAPT colWithVal = {colIdx, val};
                     varsToBeFixed.push_back(colWithVal);
                  }
               }
            }
            // using colAdaptLinkBlock, go through the columns in this child:

            // go through the columns in Bmat, Blmat of the equality:
            updateCPforAdaptFixationsBChild( it, EQUALITY_SYSTEM );
            adaptOtherSystemChildB( EQUALITY_SYSTEM, varsToBeFixed, newSREq );

            // and go through the columns in Bmat, Blmat of the inequality:
            updateCPforAdaptFixationsBChild( it, INEQUALITY_SYSTEM );
            adaptOtherSystemChildB( INEQUALITY_SYSTEM, varsToBeFixed, newSRIneq );
         }
         currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
      }

      // Update nRowLink and lhs/rhs (Linking part) of both systems:
      updateRhsNRowLink();
      presData.resetRedCounters();

      // F_0 block:
      setCurrentPointersToNull();
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
      setCPBlmatsRoot(presProb->A);
      setCPColumnRoot();
      currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vecl);
      for( int rowIdx = 0; rowIdx < currNnzRow->n; rowIdx++ )
      {
         if( currNnzRow->elements()[rowIdx] == 1.0 && (currBlmat->rowptr[rowIdx].start != currBlmat->rowptr[rowIdx].end) )
         {
            assert( currBlmat->rowptr[rowIdx].start +1 == currBlmat->rowptr[rowIdx].end );

            findSingletonsB0Equality( *currBlmat, rowIdx);
         }
      }
      if( !presData.combineColAdaptParent() )
         abortInfeasible(MPI_COMM_WORLD);

      // todo delete...summed up two times! (after function call)
      updateLinkingVarsBlocks(newSREq, newSRIneq);
      presData.resetRedCounters();
   }
}

/* Finds and stores fixations in colAdapParent, deletes singletons in INEQUALITY_SYSTEM and finds and stores fixations resulting from
 * bound changes. Stores the newly found bounds in newBOundsParent.
 *
 * Later both of the stored changes have to be applied to the whole system.
 */
void StochPresolverSingletonRows::procSingletonRowRoot(SystemType system_type) {
	/* B0 node */
	processSingletonBlock(system_type, CHILD_BLOCK, -1);

	/* linking vars Bl */
	processSingletonBlock(system_type, LINKING_CONS_BLOCK, -1);
}

/* Processing the singleton rows in child it, more precisely, goes through all singleton rows in Amat and Bmat.
 * Those in Amat are stored in colAdaptParent for later processing.
 * Those in Bmat are removed and stored in colAdaptLinkBlock. Furthermore, the corresponding fixed variables (columns)
 * in Bmat and in Blmat are removed.
 * Using this colAdaptLinkBlock, the variables (columns) are removed from the inequalities Bmat, Blmat as well.
 */
void StochPresolverSingletonRows::procSingletonRowChild(int node, int& n_singleton_sys, int& n_singleton_other_sys, SystemType system_type)
{
	if( nodeIsDummy(it, system_type) )
		return;

	/* Amat */
	processSingletonBlock(system_type, LINKING_VARS_BLOCK, node);

   std::vector<COLUMNTOADAPT> linking_cols_to_adapt;
   procSingletonRowChildNonLinking(it, system_type, linking_cols_to_adapt);

   // using colAdaptLinkBlock, go through the columns in Blmat
   if( hasLinking(system_type) )
      adaptChildBlmat( linking_cols_to_adapt, system_type);

   SystemType other_sys = (system_type == EQUALITY_SYSTEM) ? INEQUALITY_SYSTEM : EQUALITY_SYSTEM;
   updateCPforAdaptFixationsBChild( it, other_sys );

   adaptOtherSystemChildB( other_sys, linking_cols_to_adapt, n_singleton_other_sys );
}

/** Finds singleton rows in the specified block and processes them
 *
 * If the processed block is from the equality system, found singleton entries will get stored for later sync and deletion.
 * If the processed block is from the inequality system, found singleton rows will get deleted, the deletion is stored (for an update of the non-zero
 * counters) and bounds are adapted accordingly.
 *
 * Method works for all types of blocks. Synchronization has to be done accordingly.
 */
void StochPresolverSingletonRows::processSingletonBlock(SystemType system_type, BlockType block_type, int node)
{
	if(block_type == LINKING_VARS_BLOCK)
		if(!hasLinking(system_type))
			return;

	int myRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	updatePointersForCurrentNode(node, system_type);

	double* ixlow = (block_type == LINKING_VARS_BLOCK) ? currIxlowParent->elements() : currIxlowChild->elements();
	double* ixupp = (block_type == LINKING_VARS_BLOCK) ? currIxuppParent->elements() : currIxuppParent->elements();
	double* xlow = (block_type == LINKING_VARS_BLOCK) ? currxlowParent->elements() : currxlowChild->elements();
	double* xupp = (block_type == LINKING_VARS_BLOCK) ? currxuppParent->elements() : currxuppChild->elements();

	SimpleVector* nnz_row = (block_type == LINKING_CONS_BLOCK) ? currNnzRowLink : currNnzRow;
	SimpleVector* nnz_col = (block_type == LINKING_VARS_BLOCK) ? currNnzColParent : currNnzColChild;

	SparseStorageDynamic* matrix;
	SparseStorageDynamic* matrix_transp;
	if(block_type == LINKING_VARS_BLOCK)
	{
		matrix = currAmat;
		matrix_transp = currAmatTrans;
	}
	else if (block_type == LINKING_CONS_BLOCK)
	{
		matrix = currBlmat;
		matrix_transp = currBlmatTrans;
	}
	else
	{
		assert(block_type == CHILD_BLOCK);
		matrix = currBmat;
		matrix_transp = currBmatTrans;
	}

	for (int i = 0; i < nnz_row->length(); ++i) {
		if (nnz_row->elements()[i] == 1.0 && matrix->rowptr[i].start + 1 == currBlmat->rowptr[i].end )
		{
			int colIdx = -1;
			double aik = 0.0;

			getValuesForSR(*matrix, i, colIdx, aik);
			assert(!PIPSisEQ(aik, 0.0));

			if ( system_type == EQUALITY_SYSTEM )
			{
				const double val = currEqRhs->elements()[i] / aik;

				if ((ixlow[colIdx] != 0.0 && PIPSisLT(val, xlow[colIdx])) || (ixupp[colIdx] != 0.0 && PIPSisLT(xupp[colIdx], val)))
				{
					std::cout << "Singleton Row Presolving detected infeasibility : fixation of variable to invalid value in EQUALITY_SYSTEM B0" << std::endl;
					std::cout << "variable index: " << colIdx << "\tvalue: " << val	<< "\tlower bound: " << xlow[colIdx] << "\tupper bound: " << xupp[colIdx] << std::endl;
					abortInfeasible(MPI_COMM_WORLD);
				}
				else
				{
					if ( myRank == 0 || node != -1 )
						storeColValInColAdaptParent(colIdx, val);
				}
			}
			else if (system_type == INEQUALITY_SYSTEM)
			{
				double newxlow = -std::numeric_limits<double>::max();
				double newxupp = std::numeric_limits<double>::max();
				double val = 0.0;

				calculateNewBoundsOnVariable(newxlow, newxupp, i, aik);

				if ( newBoundsImplyInfeasible(newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
					abortInfeasible(MPI_COMM_WORLD);
				else if ( newBoundsFixVariable(val, newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
				{
					/* let only rank 0 store the deletions in the root node - these will be synchronized afterwards */
					if ( myRank == 0 || node != -1 )
						storeColValInColAdaptParent(colIdx, val);
				}
				else
				{
					tightenBounds(newxlow, newxupp, ixlow[colIdx], xlow[colIdx], ixupp[colIdx], xupp[colIdx]);
					// set a_ik=0.0, nRow=0, nCol--
					clearRow(*matrix, i);
					// remove entry a_ik in transposed matrix as well
					removeEntryInDynamicStorage(*matrix_transp, colIdx, i, val);
					nnz_row->elements()[i]--; // todo cannot ajust them right here!
					nnz_col->elements()[colIdx]--;

					//      assert( verifyNnzcounters() ); // todo : can this handl dynamic storage too?
					assert(currNnzRow->elements()[i] == 0);
				}
			}

		}
	}
}

bool StochPresolverSingletonRows::tightenBounds(double new_xlow, double new_xupp, double& ixlow, double& old_xlow, double& ixupp, double& old_xupp) const
{
	bool tightened = false;

	if (ixlow != 0.0 && PIPSisLT(old_xlow, new_xlow))
	{
		old_xlow = new_xlow;
	}
	else if ( ixlow == 0.0 && new_xlow > -std::numeric_limits<double>::max() )
	{
		old_xlow = new_xlow;
		ixlow = 1.0;
	}

	if (ixupp != 0.0 && PIPSisLT(new_xupp, old_xupp) )
	{
		old_xupp = new_xupp;
	}
	else if (ixupp == 0.0 && new_xupp < std::numeric_limits<double>::max() )
	{
		old_xupp = new_xupp;
		ixupp = 1.0;
	}
}

void StochPresolverSingletonRows::calculateNewBoundsOnVariable(double& newxlow, double& newxupp, int rowIdx, double aik) const
{
   if( PIPSisLT(0.0, aik) )
   {
      if( currIclow->elements()[rowIdx] != 0.0 )
         newxlow = currIneqLhs->elements()[rowIdx] / aik;
      if( currIcupp->elements()[rowIdx] != 0.0 )
         newxupp = currIneqRhs->elements()[rowIdx] / aik;
   }
   else
   {
      if( currIcupp->elements()[rowIdx] != 0.0 )
         newxlow = currIneqRhs->elements()[rowIdx] / aik;
      if( currIclow->elements()[rowIdx] != 0.0 )
         newxupp = currIneqLhs->elements()[rowIdx] / aik;
   }
}

/** Should be called right after doSingletonRowsC() or another method that stores
 * information to update in newBoundsParent.
 * Updates the bounds on the linking variables.
 */
void StochPresolverSingletonRows::updateLinkingVarsBounds()
{
   setCPColumnRoot();
   double* ixlow = currIxlowParent->elements();
   double* ixupp = currIxuppParent->elements();
   double* xlow = currxlowParent->elements();
   double* xupp = currxuppParent->elements();

   // apply updated newBoundsParent to the variable bounds.
   for(int i=0; i<getNumberNewBoundsParent(); i++)
   {
      XBOUNDS newbounds = getNewBoundsParent(i);
      setNewBounds(newbounds.colIdx, newbounds.newxlow, newbounds.newxupp, ixlow, xlow, ixupp, xupp);
   }
   clearNewBoundsParent();
}

void StochPresolverSingletonRows::getValuesForSR(SparseStorageDynamic const & storage, int rowIdx, int& colIdx, double& aik) const
{
   const int indexK = storage.rowptr[rowIdx].start;
   colIdx = storage.jcolM[indexK];
   aik = storage.M[indexK];

   assert(storage.rowptr[rowIdx].start +1 == storage.rowptr[rowIdx].end);
   assert(aik != 0.0);
}
