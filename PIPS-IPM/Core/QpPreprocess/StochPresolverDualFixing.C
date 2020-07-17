/*
 * StochPresolverDualFixing.C
 *
 *  Created on: 02.05.2018
 *      Author: bzfkempk
 */

#include "StochPresolverDualFixing.h"
#include <vector>

StochPresolverDualFixing::StochPresolverDualFixing(PresolveData& pres_data, const sData& orig_prob) :
   StochPresolverBase(pres_data, orig_prob), fixed_columns(0)
{
}

StochPresolverDualFixing::~StochPresolverDualFixing()
{
}

bool StochPresolverDualFixing::applyPresolving()
{
   assert(presData.reductionsEmpty());
   assert(presData.presDataInSync());

   int n_fixed_cols_run = 0;

#ifndef NDEBUG
   if( my_rank == 0 && verbosity > 1 )
   {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << "--- Before dual fixing:" << std::endl;
   }
   countRowsCols();
#endif

   /* check root node */
   /* allreduce findings and remove synced */
   n_fixed_cols_run = applyDualFixingNode(-1);


   /* check all other nodes if not dummy */

   for( int i = 0; i < nChildren; ++i )
   {
      if( !presData.nodeIsDummy(i) )
      {
         n_fixed_cols_run += applyDualFixingNode(i);
      }
   }


#ifndef NDEBUG
   if( my_rank == 0 && verbosity > 1 )
   {
      std::cout << "\tFixed columns in dual fixing: " << fixed_columns << std::endl;
      std::cout << "--- After dual fixing:" << std::endl;
   }
   else if( my_rank == 0 && verbosity == 1)
      std::cout << "DuFix:\t removed " << fixed_columns << " cols" << std::endl;

   countRowsCols();
   if( my_rank == 0 && verbosity > 1 )
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   assert(presData.reductionsEmpty());
   assert(presData.presDataInSync());

   if( n_fixed_cols_run != 0 )
      return true;
   else
      return false;

}

int StochPresolverDualFixing::applyDualFixingNode(int node)
{
   if( node == -1 )
      return applyDualFixingRoot();
   else
      return applyDualFixingNonRoot(node);
}

int StochPresolverDualFixing::applyDualFixingRoot()
{
   assert( currAmatTrans );
   const int n_amat = currAmat->getN();
   assert( n_amat == currAmatTrans->getM() );
   std::vector<bool> a_mat_cands(n_amat, true);


   return 0;
}

int StochPresolverDualFixing::applyDualFixingNonRoot(int node)
{
   assert( node != -1 );

   /* go through A_MAT^T and B_MAT^T (column wise in the EQUALITY_SYSTEM) and determine candidates for dual fixing */
   updatePointersForCurrentNode(node, EQUALITY_SYSTEM);
   assert( currBmatTrans );

   const int n_bmat = currBmat->getN();
   assert( n_bmat == currBmatTrans->getM() );

   std::vector<bool> b_mat_cands(n_bmat, true);

   const bool eq_linking = presData.hasLinking(EQUALITY_SYSTEM);
   for( int i = 0; i < n_bmat; ++i )
   {
      if( currBmatTrans->getRowPtr(i).start != currBmatTrans->getRowPtr(i).end )
         b_mat_cands[i] = false;

      if( eq_linking )
      {
         assert( currBlmatTrans );
         if( currBlmatTrans->getRowPtr(i).start != currBlmatTrans->getRowPtr(i).end )
            b_mat_cands[i] = false;
      }
   }

   /* go though the INEQUALITY_SYSTEM and check whether columns actually qualify for dual fixing and fix them */
   updatePointersForCurrentNode(node, INEQUALITY_SYSTEM);
   assert( currBmatTrans );
   assert( n_bmat == currBmatTrans->getM() );

   const bool ineq_linking = presData.hasLinking(INEQUALITY_SYSTEM);

   /* check the columns */
   for( int i = 0; i < n_bmat; ++i )
   {
      if( !b_mat_cands[i] )
         continue;
      const double obj_coeff = (*currgChild)[i];
      const bool pos_coeff_obj = PIPSisLT(0.0, obj_coeff);

      for( int j = currBmatTrans->getRowPtr(i).start; j < currBmatTrans->getRowPtr(i).end; ++j )
      {
         const double entry = currBmatTrans->getMat(j);
         const double row = currBmatTrans->getJcolM(j);
         /* check orientation of row */
         if( !PIPSisZero( (*currIcupp)[row] ) && !PIPSisZero( (*currIclow)[row]) )
         {
            b_mat_cands[i] = false;
            break;
         }
         assert( !PIPSisZero( (*currIcupp)[row] ) || !PIPSisZero( (*currIclow)[row] ) );

         /* convert row into <= row */
         const bool pos_coeff_row = PIPSisZero( (*currIclow)[row] ) ? pos_coeff_obj : !pos_coeff_obj;

         if( pos_coeff_row && !PIPSisLT( 0.0, entry ) )
         {
            b_mat_cands[i] = false;
            break;
         }
         else if( !PIPSisLT( entry, 0.0 ) )
         {
            b_mat_cands[i] = false;
            break;
         }
      }

      /* if we already determinde that this col was not a candidate - continue */
      if( !b_mat_cands[i] )
         continue;

      if( ineq_linking )
      {
         assert( currBlmatTrans );
         for( int j = currBlmatTrans->getRowPtr(i).start; j < currBlmatTrans->getRowPtr(i).end; ++j )
         {
            const double entry = currBlmatTrans->getMat(j);
            const double row = currBlmatTrans->getJcolM(j);
            /* check orientation of row */
            if( !PIPSisZero( (*currIcuppLink)[row] ) && !PIPSisZero( (*currIclowLink)[row] ) )
            {
               b_mat_cands[i] = false;
               break;
            }
            assert( !PIPSisZero( (*currIcuppLink)[row] ) || !PIPSisZero( (*currIclowLink)[row] ) );

            /* convert row into <= row */
            const bool pos_coeff_row = PIPSisZero( (*currIclow)[row] ) ? pos_coeff_obj : !pos_coeff_obj;

            if( pos_coeff_row && !PIPSisLT( 0.0, entry ) )
            {
               b_mat_cands[i] = false;
               break;
            }
            else if( !PIPSisLT( entry, 0.0 ) )
            {
               b_mat_cands[i] = false;
               break;
            }
         }
      }
   }

   /* after finding all candidates - fix them */
   int n_fixed = 0;
   for( int i = 0; i < n_bmat; ++i )
   {
      if( !b_mat_cands[i] )
         continue;

      const double obj_coeff = (*currgChild)[i];

      if( PIPSisLT( 0.0, obj_coeff ) )
      {
         if( PIPSisZero( (*currIxlowChild)[i] ) )
            PIPS_MPIabortInfeasible("Duals fixing found unbounded variable! Unbounded or Infeasible!",
               "StochPresolverDualFixing.C", "applyDualFixingNonRoot");
         else if( PIPSisZero(obj_coeff) )
         {
            // TODO : what is small enough here...
            // We want to remove the variable and all constraints it appears in - set variable to something very small
            presData.fixColumn( INDEX(COL, node, i), 1e-20 );
         }
         else
         {
            // TODO : maybe postsolve needs a different method here - not sure
            presData.fixColumn( INDEX(COL, node, i), (*currxlowChild)[i] );
         }
      }
      else
      {
         if( PIPSisZero( (*currIxuppChild)[i] ) )
            PIPS_MPIabortInfeasible("Duals fixing found unbounded variable! Unbounded or Infeasible!",
               "StochPresolverDualFixing.C", "applyDualFixingNonRoot");
         else if( PIPSisZero(obj_coeff) )
         {
            // TODO : what is small enough here...
            // We want to remove the variable and all constraints it appears in - set variable to something very small
            presData.fixColumn( INDEX(COL, node, i), 1e20 );
         }
         else
         {
            // TODO : maybe postsolve needs a different method here - not sure
            presData.fixColumn( INDEX(COL, node, i), (*currxuppChild)[i] );
         }
      }
      ++n_fixed;
   }

   return n_fixed;
}
