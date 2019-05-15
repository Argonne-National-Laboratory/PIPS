/*
 * StochPresolverBoundStrengthening.C
 *
 *  Created on: 28.05.2018
 *      Author: Svenja Uslu
 */

#include "StochPresolverBoundStrengthening.h"
#include <limits>
#include <cmath>
#include "pipsdef.h"

// todo exhaustive

StochPresolverBoundStrengthening::StochPresolverBoundStrengthening(
      PresolveData& presData, const sData& origProb) :
      StochPresolverBase(presData, origProb)
{
   // todo
}

StochPresolverBoundStrengthening::~StochPresolverBoundStrengthening()
{
   // todo
}

// todo print variables bounds strengthened
// todo print removed fixed cols
void StochPresolverBoundStrengthening::applyPresolving()
{
   assert(presData.reductionsEmpty());
   assert(presData.presProb->isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(indivObjOffset == 0.0);
   assert(newBoundsParent.size() == 0);

#ifndef NDEBUG
   if( my_rank == 0 )
   {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << "--- Before Bound Strengthening Presolving:" << std::endl;
   }
   countRowsCols();
#endif

   if( my_rank == 0 )
      std::cout << "Start Bound Strengthening Presolving..." << std::endl;

   do  // todo
   {
      /* root nodes */
      strenghtenBoundsInNode( EQUALITY_SYSTEM, -1);
      strenghtenBoundsInNode( INEQUALITY_SYSTEM, -1);

      // children:
      for( int node = 0; node < nChildren; node++)
      {
         // dummy child?
         if( !presData.nodeIsDummy(node, EQUALITY_SYSTEM) )
            strenghtenBoundsInNode(EQUALITY_SYSTEM, node);

         if( !presData.nodeIsDummy(node, INEQUALITY_SYSTEM) )
            strenghtenBoundsInNode(INEQUALITY_SYSTEM, node);
      }

   }
   while( false ); // todo exhaustive

   /* update bounds on all processors */
   presData.allreduceLinkingVarBounds();

   if( my_rank == 0 )
      std::cout << "Global objOffset is now: " << presData.getObjOffset() << std::endl;

#ifndef NDEBUG
   if( my_rank == 0 )
      std::cout << "--- After bound strengthening presolving:" << std::endl;
   countRowsCols();
   if( my_rank == 0 )
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   assert(presData.reductionsEmpty());
   assert(presData.presProb->isRootNodeInSync());
   assert(verifyNnzcounters());
   assert(indivObjOffset == 0.0);
   assert(newBoundsParent.size() == 0);
}

void StochPresolverBoundStrengthening::strenghtenBoundsInNode(SystemType system_type, int node)
{
   assert( -1 <= node && node < nChildren );

   strenghtenBoundsInBlock(system_type, node, LINKING_VARS_BLOCK);
   if(hasLinking(system_type))
      strenghtenBoundsInBlock(system_type, node, LINKING_CONS_BLOCK);

   if(node != -1)
      strenghtenBoundsInBlock(system_type, node, CHILD_BLOCK);

}

/**
 * Strengthen the variable bounds in the block matrix. If childBlock==true, then a block B_i or D_i is considered.
 * partMinActivity and partMaxActivity represent the partial row activity of the respective other block.
 */
void StochPresolverBoundStrengthening::strenghtenBoundsInBlock( SystemType system_type, int node, BlockType block_type)
{
   assert( 0 <= node && node < nChildren );
   updatePointersForCurrentNode(node, system_type);

   SimpleVector& xlow = (node == -1 || block_type == LINKING_VARS_BLOCK) ? *currxlowParent : *currxlowChild;
   SimpleVector& ixlow = (node == -1 || block_type == LINKING_VARS_BLOCK) ? *currIxlowParent : *currIxlowChild;
   SimpleVector& xupp = (node == -1 || block_type == LINKING_VARS_BLOCK) ? *currxuppParent : *currxuppChild;
   SimpleVector& ixupp = (node == -1 || block_type == LINKING_VARS_BLOCK) ? *currIxuppParent : *currIxuppChild;

   SimpleVector& iclow = (block_type == LINKING_CONS_BLOCK) ? *currIclowLink : *currIclow;
   SimpleVector& clow = (block_type == LINKING_CONS_BLOCK) ? *currIneqLhsLink : *currIneqLhs;
   SimpleVector& icupp = (block_type == LINKING_CONS_BLOCK) ? *currIcuppLink : *currIcupp;
   SimpleVector& cupp = (block_type == LINKING_CONS_BLOCK) ? *currIneqRhsLink : *currIneqRhs;

   SimpleVector& rhs = (block_type == LINKING_CONS_BLOCK) ? *currEqRhsLink : *currEqRhs;

   SimpleVector& actmin = (block_type == LINKING_CONS_BLOCK) ? *currActMinLink : *currActMin;
   SimpleVector& actmax = (block_type == LINKING_CONS_BLOCK) ? *currActMaxLink : *currActMax;

   SparseStorageDynamic* mat;
   if(block_type == LINKING_CONS_BLOCK)
      mat = currBlmat;
   else if(block_type == LINKING_VARS_BLOCK)
      mat = currAmat;
   else
   {
      assert(node != -1);
      mat = currBmat;
   }
   assert(mat);

   for(int row = 0; row < mat->m; ++row)
   {
      const double actmin_row = actmin[row];
      const double actmax_row = actmax[row];

      // todo if actmin / actmax too high then the bounds deduced from them are not of precision feastol anymore
      if( (actmin_row <= -std::numeric_limits<double>::max() ) && (row_activity_max >= std::numeric_limits<double>::max()) )
         return;

      for( int j = mat->rowptr[row].start; j < mat->rowptr[row].end; j++ )
      {
         // compute the possible new bounds on variable x_colIdx:
         const int col = mat->jcolM[j];
         const double a_ik = mat->M[j];

         assert( !PIPSisZero(a_ik) );
         assert( ixlow[col] != 0.0 || ixupp[col] != 0.0 );

         /* subtract current entry from row activity */
         const double actmin_row_without_curr = ( PIPSisLE(a_ik, 0) ) ? actmin_row - a_ik * xupp[col] : actmin_row - a_ik * xlow[col];
         const double actmax_row_without_curr = ( PIPSisLE(a_ik, 0) ) ? actmax_row - a_ik * xlow[col] : actmax_row - a_ik * xupp[col];

         double lbx_new = -std::numeric_limits<double>::max();
         double ubx_new = std::numeric_limits<double>::max();

         // todo : add safguard - if actmin_row... is too big and a_ik too small the computation will not make much sense anymore
         if(system_type == EQUALITY_SYSTEM)
         {
            if( PIPSisLT(0.0, a_ik) )
            {
               if( actmin_row_without_curr > -std::numeric_limits<double>::max() )
                  ubx_new = (rhs[row] - actmin_row_without_curr) / a_ik;
               if(actmax_row_without_curr < std::numeric_limits<double>::max())
                  lbx_new = (rhs[row] - actmax_row_without_curr) / a_ik;
            }
            else
            {
               if(actmin_row_without_curr > -std::numeric_limits<double>::max())
                  lbx_new = (rhs[row] - actmin_row_without_curr) / a_ik;
               if(actmax_row_without_curr < std::numeric_limits<double>::max())
                  ubx_new = (rhs[row] - actmax_row_without_curr) / a_ik;
            }
         }
         else
         {
            if( PIPSisLT(0.0, a_ik) )
            {
               if(actmin_row_without_curr > -std::numeric_limits<double>::max() && icupp[row] != 0.0)
                  ubx_new = (cupp[row] - actmin_row_without_curr) / a_ik;
               if(actmax_row_without_curr < std::numeric_limits<double>::max() && iclow[row] != 0.0)
                  lbx_new = (clow[row] - actmax_row_without_curr) / a_ik;
            }
            else
            {
               if(actmin_row_without_curr > -std::numeric_limits<double>::max() && icupp[row] != 0.0)
                  lbx_new = (cupp[row] - actmin_row_without_curr) / a_ik;
               if(actmax_row_without_curr < std::numeric_limits<double>::max() && iclow[row] != 0.0)
                  ubx_new = (clow[row] - actmax_row_without_curr) / a_ik;
            }
         }

         presData.adjustMatrixBoundsBy(system_type, node, block_type, row, col, value);
         /* check if new bounds valid */ // todo
         if( newBoundsImplyInfeasible(newBoundLow, newBoundUpp, colIdx, ixlow.elements(), ixupp.elements(), xlow.elements(), xupp.elements()) )
            abortInfeasible(MPI_COMM_WORLD);

         /* check if new bounds imply fixation of the variable */
         double fixation_value = 0.0;
         if( newBoundsFixVariable(fixation_value, newBoundLow, newBoundUpp, colIdx, ixlow.elements(),
               ixupp.elements(), xlow.elements(), xupp.elements()) )
         {
            /* for Amat we store deletions - collect them and apply them later */
             if( node == -1 || block_type == LINKING_VARS_BLOCK )
             {
                /* only 0 process stores fixations in root node B0/Bl0 - this is to reduce MPI communications - it is not necessary */
                if( myRank == 0 && node == -1 )
                   storeColValInColAdaptParent(colIdx, fixation_value);
                else if( block_type == LINKING_VARS_BLOCK )
                   storeColValInColAdaptParent(colIdx, fixation_value);
             }
             else
             {
                /* delete variable */
                deleteNonlinkColumnFromSystem(node, colIdx, fixation_value);
             }
         }
         else
         {
            tightenBounds(newBoundLow, newBoundUpp, ixlow[colIdx], xlow[colIdx], ixupp[colIdx], xupp[colIdx]);
         }
   }
   }
}

