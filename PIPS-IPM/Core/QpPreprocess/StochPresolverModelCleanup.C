/*
 * StochPresolverModelCleanup.C
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#include "StochPresolverModelCleanup.h"
#include <cmath>
#include <utility>
#include <vector>
#include <string>

StochPresolverModelCleanup::StochPresolverModelCleanup(PresolveData& presData, const sData& origProb)
   : StochPresolverBase(presData, origProb), removed_entries_total(0), removed_rows_total(0)
{
}

StochPresolverModelCleanup::~StochPresolverModelCleanup()
{
}


void StochPresolverModelCleanup::applyPresolving()
{
   assert(presData.reductionsEmpty());
   assert(presData.presDataInSync());

#ifndef NDEBUG
   if( my_rank == 0 )
   {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << "--- Before model cleanup:" << std::endl;
   }
   countRowsCols();
#endif

   int n_removed_entries = 0;
   int n_removed_rows = 0;

   // removal of redundant constraints
   int n_removed_rows_eq = removeRedundantRows(EQUALITY_SYSTEM);
   int n_removed_rows_ineq = removeRedundantRows(INEQUALITY_SYSTEM);
   n_removed_rows = n_removed_rows_eq + n_removed_rows_ineq;

   /* remove entries from A and C matrices and updates transposed systems */
   int n_removed_entries_eq = removeTinyEntriesFromSystem(EQUALITY_SYSTEM);
   int n_removed_entries_ineq = removeTinyEntriesFromSystem(INEQUALITY_SYSTEM);
   n_removed_entries = n_removed_entries_eq + n_removed_entries_ineq;

   presData.allreduceAndApplyNnzChanges();

   fixEmptyColumns();

   // update all nnzCounters - set reductionStochvecs to zero afterwards
   presData.allreduceAndApplyBoundChanges();
   presData.allreduceObjOffset();
   presData.allreduceAndApplyLinkingRowActivities();

   if( distributed )
   {
      PIPS_MPIgetSumInPlace( n_removed_entries, MPI_COMM_WORLD);
      PIPS_MPIgetSumInPlace( n_removed_rows, MPI_COMM_WORLD);
   }
   removed_entries_total += n_removed_entries;
   removed_rows_total += n_removed_rows;

#ifndef NDEBUG
   if( my_rank == 0 )
   {
      std::cout << "\tRemoved redundant rows in model cleanup: " << removed_rows_total << std::endl;
      std::cout << "\tRemoved tiny entries in model cleanup: " << removed_entries_total << std::endl;
      std::cout << "--- After model cleanup:" << std::endl;
   }

   countRowsCols();
   if(my_rank == 0)
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   assert(presData.reductionsEmpty());
   assert(presData.presDataInSync());
}

/** Remove redundant rows in the constraint system. Compares the minimal and maximal row activity
 * with the row bounds (lhs and rhs). If a row is found to be redundant, it is removed.
 * If infeasiblity is detected, then Abort.
 */
int StochPresolverModelCleanup::removeRedundantRows(SystemType system_type)
{
   int nRemovedRows = 0;

   // root:
   int nRemovedRowsRoot = removeRedundantRows(system_type, -1);

   if( my_rank == 0 )
      nRemovedRows += nRemovedRowsRoot;

   // children:
   for( int node = 0; node < nChildren; node++)
      if( !presData.nodeIsDummy( node ) )
         nRemovedRows += removeRedundantRows(system_type, node);

   return nRemovedRows;
}

int StochPresolverModelCleanup::removeRedundantRows(SystemType system_type, int node)
{
   assert(!presData.nodeIsDummy(node));
   int n_removed_rows = 0;
   int n_removed_rows_link = 0;
   if(node == -1)
      n_removed_rows_link = removeRedundantRows(system_type, node, true);

   n_removed_rows += removeRedundantRows(system_type, node, false);

   return n_removed_rows + n_removed_rows_link;
}

int StochPresolverModelCleanup::removeRedundantRows(SystemType system_type, int node, bool linking)
{
   assert(-1 <= node && node < nChildren);
   assert( (linking && node == -1) || !linking );
   assert(!presData.nodeIsDummy(node));

   if(linking && !presData.hasLinking(system_type))
      return 0;
   updatePointersForCurrentNode(node, system_type);

   int n_removed_rows = 0;

   const SimpleVectorBase<int>& nnzs = (linking == false) ? *currNnzRow : *currNnzRowLink;
   const SimpleVector& rhs_eq = (linking == false) ? *currEqRhs : *currEqRhsLink;
   const SimpleVector& clow  = (linking == false) ? *currIneqLhs : *currIneqLhsLink;
   const SimpleVector& cupp = (linking == false) ? *currIneqRhs : *currIneqRhsLink;
   const SimpleVector& iclow = (linking == false) ? *currIclow : *currIclowLink;
   const SimpleVector& icupp = (linking == false) ? *currIcupp : *currIcuppLink;

   for( int row_index = 0; row_index < nnzs.n; ++row_index)
   {
      const INDEX row(ROW, node, row_index, linking, system_type);
      if( presData.wasRowRemoved( row ) )
      {
         assert(nnzs[row_index] == 0);
         continue;
      }

      double actmin_part, actmax_part;
      int actmin_ubndd, actmax_ubndd;
      presData.getRowActivities(row, actmax_part, actmin_part, actmax_ubndd, actmin_ubndd);

      if( system_type == EQUALITY_SYSTEM )
      {
         if( actmin_ubndd != 0 || actmax_ubndd != 0)
            continue;

         if( (PIPSisLT( rhs_eq[row_index], actmin_part, feastol) && actmin_ubndd == 0)  || (PIPSisLT(actmax_part, rhs_eq[row_index], feastol) && actmax_ubndd == 0))
         {
            PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Found row that cannot meet it's rhs with it's computed activities", "StochPresolverModelCleanup.C",
                  "removeRedundantRows");
         }
         else if( PIPSisLE(rhs_eq[row_index], actmin_part, feastol) && PIPSisLE(actmax_part, rhs_eq[row_index], feastol) )
         {
            presData.removeRedundantRow( row );
            n_removed_rows++;
         }
      }
      else
      {
         assert(!presData.wasRowRemoved( row ));
         assert( PIPSisLT(0.0, iclow[row_index] + icupp[row_index]) );

         if( ( !PIPSisZero(iclow[row_index]) && (actmax_ubndd == 0 && PIPSisLTFeas(actmax_part, clow[row_index])) )
               || ( !PIPSisZero(icupp[row_index]) && (actmin_ubndd == 0 && PIPSisLTFeas(cupp[row_index], actmin_part)) ) )
            PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Found row that cannot meet it's lhs or rhs with it's computed activities", "StochPresolverModelCleanup.C",
                  "removeRedundantRows");
         else if( ( PIPSisZero(iclow[row_index]) || PIPSisLE(clow[row_index], -infinity) ) &&
               ( PIPSisZero(icupp[row_index]) || PIPSisLE(infinity, cupp[row_index])) )
         {
            presData.removeRedundantRow( row );
            n_removed_rows++;
         }
         else if( ( PIPSisZero(iclow[row_index]) || (actmin_ubndd == 0 && PIPSisLEFeas(clow[row_index], actmin_part)) )
               && ( PIPSisZero(icupp[row_index]) || (actmax_ubndd == 0 && PIPSisLEFeas(actmax_part, cupp[row_index])) ) )
         {
            presData.removeRedundantRow( row );
            n_removed_rows++;
         }
      }
   }
   return n_removed_rows;
}

/* removes all small entries from the specified system
 *
 * While removing the reduction vectors in presData get set - nnzVectors are not updated an that must be
 * done using updateNnzFromReductions if needed.
 * Transposed matrices get updated in a subroutine - so after calling this method, that matrix should
 * be in a consistent state.
 */
int StochPresolverModelCleanup::removeTinyEntriesFromSystem(SystemType system_type)
{
   assert(dynamic_cast<const StochGenMatrix&>(*(presData.getPresProb().A)).children.size() == (size_t) nChildren);
   assert(dynamic_cast<const StochGenMatrix&>(*(presData.getPresProb().C)).children.size() == (size_t) nChildren);

   int n_elims = 0;

   assert(!presData.nodeIsDummy(-1));

   /* reductions in root node */
   /* process B0 and Bl0 */
   n_elims += removeTinyInnerLoop(system_type, -1, B_MAT);
   if (presData.hasLinking(system_type))
	n_elims += removeTinyInnerLoop(system_type, -1, BL_MAT);


   /* count eliminations in B0 and Bl0 only once */
   if( distributed && my_rank != 0 )
      n_elims = 0;

   // go through the children
   for( int node = 0; node < nChildren; node++ )
   {
      if( !presData.nodeIsDummy(node) )
      {
         /* Amat */
         n_elims += removeTinyInnerLoop(system_type, node, A_MAT );

         /* Bmat */
         n_elims += removeTinyInnerLoop(system_type, node, B_MAT );

         /* this has to be synchronized */
         /* Blmat */
         if( presData.hasLinking(system_type) )
            n_elims += removeTinyInnerLoop(system_type, node, BL_MAT);
         }
   }

   return n_elims;
}

/** Removes tiny entries in storage and adapts the rhs accordingly.
 *  If block_type == LINKING_VARS_BLOCK, then block Amat is considered.
 *  If block_type == CHILD_BLOCK, then block Bmat is considered. */
/* system type indicates matrix A or C, block_type indicates the block */
// todo what is the proper order for criterion 1 to 3?
// todo for criterion 3 - should ALL eliminations be considered?
int StochPresolverModelCleanup::removeTinyInnerLoop( SystemType system_type, int node, BlockType block_type)
{
   if(presData.nodeIsDummy(node))
      return 0;

   updatePointersForCurrentNode(node, system_type);

   const SparseStorageDynamic* mat = nullptr;

   const SimpleVector* x_lower = nullptr;
   const SimpleVector* x_lower_idx = nullptr;
   const SimpleVector* x_upper = nullptr;
   const SimpleVector* x_upper_idx = nullptr;
   const SimpleVectorBase<int>* nnzRow = nullptr;

   /* set matrix */
   if( block_type == B_MAT )
   {
      mat = currBmat;
      assert(currBmatTrans);
   }
   else if( block_type == A_MAT )
   {
      assert(node != -1);
      mat = currAmat;
      assert(currAmatTrans);
   }
   else if( block_type == BL_MAT)
   {
      mat = currBlmat;
      assert(currBlmatTrans);
   }

   /* set variables */
   x_lower = ( block_type == A_MAT || node == -1) ? currxlowParent : currxlowChild;
   x_lower_idx = ( block_type == A_MAT || node == -1) ? currIxlowParent : currIxlowChild;
   x_upper = ( block_type == A_MAT || node == -1) ? currxuppParent : currxuppChild;
   x_upper_idx = ( block_type == A_MAT || node == -1) ? currIxuppParent : currIxuppChild;

   /* set reduction vectors */

   /* set non-zero row vectors */
   nnzRow = (block_type == BL_MAT) ? currNnzRowLink : currNnzRow;

   int n_elims = 0;

   const SparseStorageDynamic* storage = mat;

   std::vector<std::pair<int, int> > eliminated_entries;

   /* for every row in row in matrix */
   for( int r = 0; r < storage->getM(); r++ )
   {
      double total_sum_modifications_row = 0.0;

      int start = storage->getRowPtr(r).start;
      int end = storage->getRowPtr(r).end;

      /* for every nonzero column in that row */
      for(int k = start; k < end; ++k )
      {
         const int col = storage->getJcolM(k);
         const double mat_entry = storage->getMat(k);

         /* remove all small entries */
         if( fabs( mat_entry ) < tol_matrix_entry )
         {
            presData.deleteEntry(system_type, node, block_type, r, k, end);

            std::pair<int,int> entry(r, col);
            eliminated_entries.push_back(entry);
            ++n_elims;
         }
         /* remove entries where their corresponding variables have valid lower and upper bounds, that overall do not have a real influence though */
         else if( !PIPSisZero((*x_upper_idx)[col]) && !PIPSisZero((*x_lower_idx)[col]) )
         {
            if( (fabs( mat_entry ) < tolerance1 && fabs( mat_entry ) * ( (*x_upper)[col] - (*x_lower)[col]) * (*nnzRow)[r] < tolerance2 * feastol ))
            {
               presData.deleteEntry(system_type, node, block_type, r, k, end);

               std::pair<int,int> entry(r, col);
               eliminated_entries.push_back(entry);
               ++n_elims;
            }
         }
         else if( false ) //todo if not linking constraints
         {
            // todo third criterion? for linking constraints: call extra function to know whether we have linking cons
            // that link only two blocks (not so urgent for linking)
            /* if valid lower and upper bounds */
            if( !PIPSisZero((*x_upper_idx)[col]) && !PIPSisZero((*x_lower_idx)[col]) )
            {
               if( total_sum_modifications_row + (fabs(mat_entry) * ((*x_upper)[col] - (*x_lower)[col])) < 1.0e-1 * feastol)
               {
                  total_sum_modifications_row += fabs(mat_entry) * ((*x_upper)[col] - (*x_lower)[col]);

                  presData.deleteEntry(system_type, node, block_type, r, k, end);

                  std::pair<int,int> entry(r, col);
                  eliminated_entries.push_back(entry);
                  ++n_elims;
               }
            }
         }
         /* not removed */
      }
   }

   presData.updateTransposedSubmatrix( system_type, node, block_type, eliminated_entries);

   assert(presData.elementsDeletedInTransposed());

   return n_elims;
}


/* Go through columns and fix all empty ones to the current variables lower/upper bound (depending on objective)
 * Might detect unboundedness of problem.
 */
void StochPresolverModelCleanup::fixEmptyColumns()
{

   for(int node = -1; node < nChildren; ++node)
   {
      if( presData.nodeIsDummy(node) )
         continue;
      
      updatePointersForCurrentNode(node, EQUALITY_SYSTEM);

      const SimpleVector& g = (node == -1) ? *currgParent : *currgChild;
      const SimpleVector& ixupp = (node == -1) ? *currIxuppParent : *currIxuppChild;
      const SimpleVector& ixlow = (node == -1) ? *currIxlowParent : *currIxlowChild;
      const SimpleVector& xupp = (node == -1) ? *currxuppParent : *currxuppChild;
      const SimpleVector& xlow = (node == -1) ? *currxlowParent : *currxlowChild;
      const SimpleVectorBase<int>& nnzs_col = (node == -1) ? *currNnzColParent : *currNnzColChild;

      for(int col_index = 0; col_index < nnzs_col.n; ++col_index)
      {
         const INDEX col(COL, node, col_index);
         /* column fixation candidate */
         if( nnzs_col[col_index] == 0)
         {
            /* check whether column was removed already */
            if( PIPSisZero(ixlow[col_index]) && PIPSisZero(ixupp[col_index])
               && PIPSisZero(xlow[col_index]) && PIPSisZero(xupp[col_index])
               && PIPSisZero(g[col_index]))
            {
               if( presData.wasColumnRemoved(col) )
                  continue;
            }

            if( PIPSisLT( g[col_index], 0.0) )
            {
               if( !PIPSisZero(ixupp[col_index]) )
               {
                  presData.fixEmptyColumn(col, xupp[col_index]);
               }
               else
               {
                  PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Found empty column with non-zero objective vector and no bounds in objective direction! Unbounded!", 
                     "StochPresolverModelCleanup.C", "fixEmptyColumns");
               }
            } 
            else if( PIPSisLT(0.0, g[col_index]) )
            {
               if( !PIPSisZero(ixlow[col_index]) )
               {
                  presData.fixEmptyColumn(col, xlow[col_index]);
               }
               else
               {
                  PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Found empty column with non-zero objective vector and no bounds in objective direction! Unbounded!", 
                     "StochPresolverModelCleanup.C", "fixEmptyColumns");
               }
            }
            else
            {
               assert( PIPSisEQ( g[col_index], 0.0) );
               presData.fixEmptyColumn(col, 0.0);
            }
         }
      }
   }
}


