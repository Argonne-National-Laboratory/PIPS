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
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());

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

   //todo : print specific stats
#ifndef NDEBUG
   if( my_rank == 0 )
   {
      std::cout << "Removed " << n_removed_rows << " redundant rows (" << n_removed_rows_eq <<
         " equalitiy and " << n_removed_rows_ineq << " inequality rows)" << std::endl;
      std::cout << "Removed " << n_removed_entries << " entries (" << n_removed_entries_eq <<
         " entries in equality system and " << n_removed_entries_ineq << " in inequality system)" << std::endl;
   }
#endif

#ifndef NDEBUG
   if( my_rank == 0 )
      std::cout << "--- After model cleanup:" << std::endl;
   countRowsCols();
   if(my_rank == 0)
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());
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

#ifndef NDEBUG
   if( my_rank == 0 && node == -1 )
   {
      n_removed_rows += n_removed_rows_link;
      std::cout << "\tRemoved redundant linking constraints during model cleanup: " << n_removed_rows_link << std::endl;
   }
#endif

   n_removed_rows += removeRedundantRows(system_type, node, false);

   return n_removed_rows + n_removed_rows_link;
}

// todo : deleting rows invalidates activities currently
int StochPresolverModelCleanup::removeRedundantRows(SystemType system_type, int node, bool linking)
{
   assert(-1 <= node && node < nChildren);
   assert( (linking && node == -1) || !linking );
   assert(!presData.nodeIsDummy(node));

   if(linking && !presData.hasLinking(system_type))
      return 0;
   updatePointersForCurrentNode(node, system_type);

   int n_removed_rows = 0;

   // todo
   const SimpleVectorBase<int>& nnzs = (linking == false) ? *currNnzRow : *currNnzRowLink;
   const SimpleVector& rhs_eq = (linking == false) ? *currEqRhs : *currEqRhsLink;
   const SimpleVector& clow  = (linking == false) ? *currIneqLhs : *currIneqLhsLink;
   const SimpleVector& cupp = (linking == false) ? *currIneqRhs : *currIneqRhsLink;
   const SimpleVector& iclow = (linking == false) ? *currIclow : *currIclowLink;
   const SimpleVector& icupp = (linking == false) ? *currIcupp : *currIcuppLink;

   for( int row = 0; row < nnzs.n; ++row)
   {
      if( nnzs[row] == 0 ) // empty rows might still have a rhs but should be ignored. // todo?
         continue;

      double actmin_part, actmax_part;
      int actmin_ubndd, actmax_ubndd;
      presData.getRowActivities(system_type, node, linking, row, actmax_part, actmin_part, actmax_ubndd, actmin_ubndd);

      if( system_type == EQUALITY_SYSTEM )
      {
         if( actmin_ubndd != 0 || actmax_ubndd != 0)
            continue;

         if( (PIPSisLT( rhs_eq[row], actmin_part, feastol) && actmin_ubndd == 0)  || (PIPSisLT(actmax_part, rhs_eq[row], feastol) && actmax_ubndd == 0))
         {
            PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Found row that cannot meet it's rhs with it's computed activities", "StochPresolverModelCleanup.C",
                  "removeRedundantRows");
         }
         else if( PIPSisLE(rhs_eq[row], actmin_part, feastol) && PIPSisLE(actmax_part, rhs_eq[row], feastol) )
         {
            presData.removeRedundantRow(system_type, node, row, linking);
            n_removed_rows++;
         }
      }
      else
      {
         if( ( !PIPSisZero(iclow[row]) && (actmax_ubndd == 0 && PIPSisLT(actmax_part, clow[row], feastol)) )
               || ( !PIPSisZero(icupp[row]) && (actmin_ubndd == 0 && PIPSisLT( cupp[row], actmin_part, feastol)) ) )
            PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Found row that cannot meet it's lhs or rhs with it's computed activities", "StochPresolverModelCleanup.C",
                  "removeRedundantRows");// todo infinity??
         else if( ( PIPSisZero(iclow[row]) || PIPSisLE(clow[row], -infinity) ) &&
               ( PIPSisZero(icupp[row]) || PIPSisLE(infinity, cupp[row])) )
         {
            presData.removeRedundantRow(system_type, node, row, linking);
            n_removed_rows++;
         }
         else if( ( PIPSisZero(iclow[row]) || (actmin_ubndd == 0 && PIPSisLE( clow[row], actmin_part, feastol)) )
               && ( PIPSisZero(icupp[row]) || (actmax_ubndd == 0 && PIPSisLE( actmax_part, cupp[row], feastol)) ) )
         {
            presData.removeRedundantRow(system_type, node, row, linking);
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


   // todo : traffic can be reduced by only zeroing the linking var col and the linking cons row
   // and not the zero row too
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
   for( int r = 0; r < storage->m; r++ )
   {
      double total_sum_modifications_row = 0.0;

      int start = storage->rowptr[r].start;
      int end = storage->rowptr[r].end;

      /* for every nonzero column in that row */
      for(int k = start; k < end; ++k )
      {
         const int col = storage->jcolM[k];
         const double mat_entry = storage->M[k];

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

      for(int col = 0; col < nnzs_col.n; ++col)
      {
         /* column fixation candidate */
         if( nnzs_col[col] == 0)
         {
            // todo : maybe this can also happen as an input? here we assume that the column has been fixed already
            if( PIPSisZero(ixlow[col]) && PIPSisZero(ixupp[col])
               && PIPSisZero(xlow[col]) && PIPSisZero(xupp[col]) 
               && PIPSisZero(g[col]))
               continue;

            if( PIPSisLT( g[col], 0.0) )
            {
               if( !PIPSisZero(ixupp[col]) )
               {
                  presData.fixEmptyColumn(node, col, xupp[col]);
               }
               else
               {
                  PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Found empty column with non-zero objective vector and no bounds in objective direction! Unbounded!", 
                     "StochPresolverModelCleanup.C", "fixEmptyColumns");
               }
            } 
            else if( PIPSisLT(0.0, g[col]) )
            {
               if( !PIPSisZero(ixlow[col]) )
               {
                  presData.fixEmptyColumn(node, col, xlow[col]);
               }
               else
               {
                  PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Found empty column with non-zero objective vector and no bounds in objective direction! Unbounded!", 
                     "StochPresolverModelCleanup.C", "fixEmptyColumns");
               }
            }
            else
            {
               assert( PIPSisEQ( g[col], 0.0) );
               presData.fixEmptyColumn(node, col, 0.0);
            }
         }
      }
   }
}


