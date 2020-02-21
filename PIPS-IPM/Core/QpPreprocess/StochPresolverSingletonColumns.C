/*
 * StochPresolverSingletonColumns.C
 *
 *  Created on: 08.05.2018
 *      Author: bzfuslus
 */

#include "StochPresolverSingletonColumns.h"
#include "StochVectorUtilities.h"

StochPresolverSingletonColumns::StochPresolverSingletonColumns(PresolveData &presData, const sData &origProb) :
   StochPresolverBase(presData, origProb), removed_cols(0), local_singletons(false),
   n_linking_rows_eq(dynamic_cast<const StochVector&>(*origProb.bA).vecl->length()),
   n_linking_rows_ineq(dynamic_cast<const StochVector&>(*origProb.bu).vecl->length()),
   local_linking_column_for_row_in_proc(n_linking_rows_eq + n_linking_rows_ineq, -1),
   cols( n_linking_rows_eq + n_linking_rows_ineq )
{
}

StochPresolverSingletonColumns::~StochPresolverSingletonColumns()
{
}

void StochPresolverSingletonColumns::applyPresolving()
{
   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());

   presData.putLinkingVarsSyncEvent();

#ifndef NDEBUG
   if( my_rank == 0 )
   {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << "--- Before singleton columns presolving:" << std::endl;
   }
   countRowsCols();
#endif

   if( presData.getSingletonCols().size() == 0 )
   {
#ifndef NDEBUG
      if( my_rank == 0 )
         std::cout << "No more singletons left - exiting" << std::endl;
#endif
   }

   int removed_cols_run = 0;

   // main loop:
   while( !presData.getSingletonCols().empty() )
   {
      const INDEX& col = presData.getSingletonCols().front();
      const bool removed = removeSingletonColumn(col);

      if( removed )
         ++removed_cols_run;

      presData.getSingletonCols().pop();
   }

   /* check for local singletons and communicate */
   PIPS_MPIgetLogicOrInPlace(local_singletons, MPI_COMM_WORLD);
   if(local_singletons)
   {
      /* allreduce local singleton columns and coeffs */
      std::cout << "local singleton columns!" << std::endl;
      // TODO: change this - it should be the one that provides most numerical stability..
      /* allreduce the procs that found a local singleton column in a linking row - the lowest ranking one will get to remove the column */
      PIPS_MPIminArrayInPlace(local_linking_column_for_row_in_proc, MPI_COMM_WORLD);

      /* remove local singleton columns */
      for(unsigned int i = 0; i < local_linking_column_for_row_in_proc.size(); ++i)
      {
         const int proc_that_removes = local_linking_column_for_row_in_proc[i];

         if(proc_that_removes == -1)
            continue;

         const SystemType system_type = (i < n_linking_rows_eq) ? EQUALITY_SYSTEM : INEQUALITY_SYSTEM;
         const int row_index = (i < n_linking_rows_eq) ? i : i - n_linking_rows_eq;
         const INDEX row(ROW, -1, row_index, true, system_type);

         if(my_rank == proc_that_removes)
         {
            presData.removeImpliedFreeColumnSingletonEqualityRowSynced( row, cols[i] );
            ++removed_cols_run;
         }
         else
         {
            presData.removeImpliedFreeColumnSingletonEqualityRowSynced( row, INDEX() );
         }
      }
   }

   //todo
   presData.allreduceAndApplyBoundChanges();
   presData.allreduceAndApplyNnzChanges();
   presData.allreduceAndApplyLinkingRowActivities();
   presData.allreduceLinkingVarBounds();
   presData.allreduceAndApplyObjVecChanges();
   presData.allreduceObjOffset();

#ifndef NDEBUG
   PIPS_MPIgetSumInPlace(removed_cols_run, MPI_COMM_WORLD);
   removed_cols += removed_cols_run;
   if( my_rank == 0 )
   {
      std::cout << "--- After singleton columns presolving:" << std::endl;
      std::cout << "\tRemoved columns during singleton column elimination: " << removed_cols << std::endl;
   }
   countRowsCols();
   if( my_rank == 0 )
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());
}

bool StochPresolverSingletonColumns::removeSingletonColumn(const INDEX& col)
{
   assert(col.isCol());

   const int node = col.getNode();
   const int col_index = col.getIndex();

   assert( -1 <= node && node < nChildren );
   assert( !presData.nodeIsDummy(node) );

   updatePointersForCurrentNode(node, EQUALITY_SYSTEM);
   const SimpleVectorBase<int> &nnzs_col = (node == -1) ? *currNnzColParent : *currNnzColChild;

   if( nnzs_col[col_index] == 0 )
      return false;

   assert(nnzs_col[col_index] == 1);

   bool found = false;
   /* find the associated row via checking transposed matrices */
   INDEX row = findRowForColumnSingleton(col, found);

   /* the singleton is a linking variable located on another process */
   if( !found )
   {
      assert( row.getType() == EMPTY_INDEX );
      assert( node == -1 );
      return false;
   }

   assert(row.isRow());
   assert(-1 <= row.getNode() && row.getNode() < nChildren);
   updatePointersForCurrentNode(row.getNode(), row.getSystemType());

   /* check whether col is free/implied free */
   bool lb_implied_free = false;
   bool ub_implied_free = false;

   // todo : we should not do anything if the coefficient of the column we are looking at is too small probably

   checkColImpliedFree(col, row, lb_implied_free, ub_implied_free);
   bool implied_free = lb_implied_free && ub_implied_free;

   /* if objective of variable is zero we can just remove it from the problem together with the containing row */
   double obj = (node == -1) ? (*currgParent)[col_index] : (*currgChild)[col_index];

   if( implied_free && PIPSisEQ(obj, 0.0) )
   {
      /* store local singleton cols for later communication */
      if( row.getLinking() && col.getNode() != -1)
      {
         local_singletons = true;

         const int index_row = (row.getSystemType() == EQUALITY_SYSTEM) ? row.getIndex() : n_linking_rows_eq + row.getIndex();
         local_linking_column_for_row_in_proc.at(index_row) = my_rank;
         cols.at(index_row) = INDEX(COL, col.getNode(), col.getIndex());
      }
      else
         presData.removeImpliedFreeColumnSingletonEqualityRow( row, col );
   }
   /* equality singleton variables */
   else if( row.getSystemType() == EQUALITY_SYSTEM && implied_free )
   {
      /* store local singleton cols for later */
      if( row.getLinking() && col.getNode() != -1)
      {
         local_singletons = true;

         const int index_row = (row.getSystemType() == EQUALITY_SYSTEM) ? row.getIndex() : n_linking_rows_eq + row.getIndex();
         local_linking_column_for_row_in_proc.at(index_row) = my_rank;
         cols.at(index_row) = INDEX(COL, col.getNode(), col.getIndex());
      }
      else
         presData.removeImpliedFreeColumnSingletonEqualityRow( row, col );
   }
   // TODO : singleton column located on other processes
   else if( row.getSystemType() == INEQUALITY_SYSTEM )
   {
      /* inequality singleton variables */
      const double iclow = getSimpleVecFromRowStochVec( *presData.getPresProb().iclow, row.getNode(), row.getLinking() )[row.getIndex()];
      const double icupp = getSimpleVecFromRowStochVec( *presData.getPresProb().icupp, row.getNode(), row.getLinking() )[row.getIndex()];
      const double clow = getSimpleVecFromRowStochVec( *presData.getPresProb().bl, row.getNode(), row.getLinking() )[row.getIndex()];
      const double cupp = getSimpleVecFromRowStochVec( *presData.getPresProb().bu, row.getNode(), row.getLinking() )[row.getIndex()];

      if( !PIPSisZero(clow) && !PIPSisZero(cupp) && PIPSisEQ(clow, cupp) )
         presData.removeImpliedFreeColumnSingletonEqualityRow( row, col );
      else if( !PIPSisZero(iclow) && !PIPSisZero(icupp) )
      {
         assert(!PIPSisEQ(clow, cupp));
         return false;
      }
      else
      {
         assert(!PIPSisEQ(iclow, icupp));

         const double obj_coeff = getSimpleVecFromColStochVec( *presData.getPresProb().g, col.getNode())[col.getIndex()];
         const double coeff = presData.getRowCoeff(row, col);

         assert(!PIPSisZero(coeff));

         /* convert row to less equal row */
         const double coeff_le_row = PIPSisZero(iclow) ? coeff : -coeff;

         const double lhsrhs = PIPSisZero(iclow) ? cupp : clow;

         /* both positive */
         if( PIPSisLE(0.0, coeff_le_row) && PIPSisLE(0.0, obj_coeff) )
         {
            const double ixlow = getSimpleVecFromColStochVec(*presData.getPresProb().ixlow, col.getNode())[col.getIndex()];
            const double xlow = getSimpleVecFromColStochVec( *presData.getPresProb().blx, col.getNode())[col.getIndex()];
            if( PIPSisZero(ixlow) && PIPSisLT(0.0, obj_coeff) )
               PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Found unbounded singleton column variable", "StochPresolverSingletonColumns.C", "removeSingletonColumn");

            if( PIPSisZero(ixlow) && PIPSisZero(obj_coeff) )
               presData.removeFreeColumnSingletonInequalityRow( row, col, lhsrhs, coeff); /* remove variable from whole system */

            /* fix variable to lower bound */
            presData.fixColumnInequalitySingleton(col, xlow, coeff);
         }
         /* both negative */
         else if( PIPSisLE(coeff, 0.0 && PIPSisLE(obj_coeff, 0.0)) )
         {
            const double ixupp = getSimpleVecFromColStochVec(*presData.getPresProb().ixupp, col.getNode())[col.getIndex()];
            const double xupp = getSimpleVecFromColStochVec( *presData.getPresProb().bux, col.getNode())[col.getIndex()];
            if( PIPSisZero(ixupp) && PIPSisLT(obj_coeff, 0.0) )
               PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Found unbounded singleton column variable", "StochPresolverSingletonColumns.C", "removeSingletonColumn");

            if( PIPSisZero(ixupp) && PIPSisZero(obj_coeff) )
               presData.removeFreeColumnSingletonInequalityRow( row, col, lhsrhs, coeff ); /* remove variable and associated row from whole system */

            /* fix variable to upper bound */
            presData.fixColumnInequalitySingleton(col, xupp, coeff);
         }
      }
   }
   else
      return false;

   if( row.getNode() != -1 || my_rank == 0 )
      return true;
   else
      return false;

}

INDEX StochPresolverSingletonColumns::findRowForColumnSingleton( const INDEX& col, bool& found )
{
   assert( col.isCol() );
   assert(-1 <= col.getNode() && col.getNode() < nChildren);

   if( col.getNode() == -1 )
      return findRowForLinkingSingleton(col.getIndex(), found);
   else
      return findRowForNonlinkingSingelton(col, found);
}

INDEX StochPresolverSingletonColumns::findRowForLinkingSingleton(int col, bool& found)
{
   /* go through all children and check linking variables for the singleton row */

   /* equality part */
   INDEX row = findRowForLinkingSingletonInSystem(col, EQUALITY_SYSTEM, found);
   if( found )
      return row;

   /* inequality part */
   return findRowForLinkingSingletonInSystem(col, INEQUALITY_SYSTEM, found);
}

INDEX StochPresolverSingletonColumns::findRowForNonlinkingSingelton(const INDEX& col, bool& found)
{
   assert(col.isCol());
   assert(col.getNode() != -1);

   if( presData.nodeIsDummy(col.getNode()) )
   {
      found = false;
      return INDEX();
   }

   /* check Bmat and Blmat for the singleton row */
   SystemType row_system_type = EQUALITY_SYSTEM;
   bool row_linking = false;
   const int row_node = col.getNode();
   int row_index = -1;

   /* equality part */
   updatePointersForCurrentNode(row_node, row_system_type);

   /* Bmat */
   assert(currBmatTrans);

   found = findRowForSingletonColumnInMatrix(*currBmatTrans, row_index, col.getIndex());
   if( found )
      return INDEX(ROW, row_node, row_index, row_linking, row_system_type);

   /* Blmat */
   assert(currBlmatTrans);
   found = findRowForSingletonColumnInMatrix(*currBlmatTrans, row_index, col.getIndex());
   if( found )
   {
      row_linking = true;
      return INDEX(ROW, row_node, row_index, row_linking, row_system_type);
   }

   /* inequality part */
   row_system_type = INEQUALITY_SYSTEM;
   updatePointersForCurrentNode(row_node, row_system_type);

   /* Bmat */
   assert(currBmatTrans);
   found = findRowForSingletonColumnInMatrix(*currBmatTrans, row_index, col.getIndex());
   if( found )
      return INDEX(ROW, row_node, row_index, row_linking, row_system_type);

   /* Blmat */
   assert(currBlmatTrans);
   found = findRowForSingletonColumnInMatrix(*currBlmatTrans, row_index, col.getIndex());
   if( found )
   {
      row_linking = true;
      return INDEX(ROW, row_node, row_index, row_linking, row_system_type);
   }

   return INDEX();
}

INDEX StochPresolverSingletonColumns::findRowForLinkingSingletonInSystem(int col, SystemType system_type, bool& found)
{

   for( int row_node = -1; row_node < nChildren; ++row_node )
   {
      if( !presData.nodeIsDummy(row_node) )
      {
         int row_index = -1;
         updatePointersForCurrentNode(row_node, system_type);

         /* check transposed for entry */
         if( row_node == -1 )
         {
            /* Bmat */
            assert(currBmatTrans);
            found = findRowForSingletonColumnInMatrix(*currBmatTrans, row_index, col);
            if( found )
               return INDEX(ROW, row_node, row_index, false, system_type);

            /* Blmat */
            assert(currBlmatTrans);
            found = findRowForSingletonColumnInMatrix(*currBlmatTrans, row_index, col);
            if( found )
               return INDEX(ROW, row_node, row_index, true, system_type);
         }
         else
         {
            /* Amat */
            assert(currAmatTrans);
            found = findRowForSingletonColumnInMatrix(*currAmatTrans, row_index, col);
            if( found )
               return INDEX(ROW, row_node, row_index, false, system_type);
         }
      }
   }
   assert(found == false);
   return INDEX();
}

bool StochPresolverSingletonColumns::findRowForSingletonColumnInMatrix(
   const SparseStorageDynamic &mat,
   int& row,
   const int& col)
{
   assert(col < mat.getM());
   if( mat.getRowPtr(col).start != mat.getRowPtr(col).end )
   {
      assert((mat.getRowPtr(col).end - mat.getRowPtr(col).start) == 1);
      row = mat.getJcolM(mat.getRowPtr(col).start);
      return true;
   }
   return false;
}

void StochPresolverSingletonColumns::checkColImpliedFree(const INDEX& col, const INDEX& row,
   bool &lb_implied_free,
   bool &ub_implied_free)
{
   assert(col.isCol());
   assert(row.isRow());

   ub_implied_free = lb_implied_free = false;

   updatePointersForCurrentNode(row.getNode(), row.getSystemType());

   /* check whether free */
   const double ixupp = getSimpleVecFromColStochVec(*presData.getPresProb().ixupp, col.getNode())[col.getIndex()];
   const double ixlow = getSimpleVecFromColStochVec(*presData.getPresProb().ixlow, col.getNode())[col.getIndex()];

   if( PIPSisZero(ixupp) )
      ub_implied_free = true;
   if( PIPSisZero(ixlow) )
      lb_implied_free = true;

   /* check whether bound tightening found bounds from the variables row that make it implied free */
   ub_implied_free = ub_implied_free || presData.varBoundImpliedFreeBy(true, col, row);
   lb_implied_free = lb_implied_free || presData.varBoundImpliedFreeBy(false, col, row);

}
