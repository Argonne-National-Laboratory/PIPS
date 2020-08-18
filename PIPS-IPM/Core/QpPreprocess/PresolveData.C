/*
 * PresolveData.C
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

#include "PresolveData.h"

#include "StochOptions.h"
#include "StochGenMatrix.h"
#include "DoubleMatrixTypes.h"
#include "StochMatrixUtilities.h"
#include "pipsdef.h"
#include "pipsport.h"
#include "StochVectorUtilities.h"

#include <stdexcept>
#include <limits>
#include <algorithm>
#include <string>
#include <cmath>

bool PresolveData::iTrackColumn() const
{
   return track_col && ( !nodeIsDummy(tracked_col.getNode()) && (my_rank == 0 || !tracked_col.isLinkingCol() ) );
}

bool PresolveData::iTrackRow() const
{
   return track_row && ( !nodeIsDummy(tracked_row.getNode()) && (my_rank == 0 || tracked_row.isLinkingRow() || tracked_row.getNode() != -1) );
}

PresolveData::PresolveData(const sData* sorigprob, StochPostsolver* postsolver) :
      postsolver(postsolver),
      limit_max_bound_accepted( pips_options::getDoubleParameter("PRESOLVE_MAX_BOUND_ACCEPTED") ),
      length_array_outdated_indicators(6),
      array_outdated_indicators(new bool[length_array_outdated_indicators]),
      outdated_lhsrhs(array_outdated_indicators[0]),
      outdated_nnzs(array_outdated_indicators[1]),
      outdated_linking_var_bounds(array_outdated_indicators[2]),
      outdated_activities(array_outdated_indicators[3]),
      outdated_obj_vector(array_outdated_indicators[4]),
      postsolve_linking_row_propagation_needed(array_outdated_indicators[5]),
      linking_rows_need_act_computation(0),
      nnzs_row_A(cloneStochVector<double,int>(*sorigprob->bA)),
      nnzs_row_C(cloneStochVector<double,int>(*sorigprob->icupp)),
      nnzs_col(cloneStochVector<double,int>(*sorigprob->g)),
      actmax_eq_part(cloneStochVector<int,double>(*nnzs_row_A)),
      actmin_eq_part(dynamic_cast<StochVector*>(actmax_eq_part->clone())),
      actmax_eq_ubndd(dynamic_cast<StochVectorBase<int>*>(nnzs_row_A->clone())),
      actmin_eq_ubndd(dynamic_cast<StochVectorBase<int>*>(nnzs_row_A->clone())),
      actmax_ineq_part(cloneStochVector<int,double>(*nnzs_row_C)),
      actmin_ineq_part(dynamic_cast<StochVector*>(actmax_ineq_part->clone())),
      actmax_ineq_ubndd(dynamic_cast<StochVectorBase<int>*>(nnzs_row_C->clone())),
      actmin_ineq_ubndd(dynamic_cast<StochVectorBase<int>*>(nnzs_row_C->clone())),
      my_rank(PIPS_MPIgetRank(MPI_COMM_WORLD)),
      distributed(PIPS_MPIgetDistributed(MPI_COMM_WORLD)),
      INF_NEG( -pips_options::getDoubleParameter("PRESOLVE_INFINITY") ),
      INF_POS( pips_options::getDoubleParameter("PRESOLVE_INFINITY") ),
      nChildren(nnzs_col->children.size()),
      track_row( pips_options::getBoolParameter("PRESOLVE_TRACK_ROW") ),
      track_col( pips_options::getBoolParameter("PRESOLVE_TRACK_COL") ),
      tracked_row(ROW, pips_options::getIntParameter("PRESOLVE_TRACK_ROW_NODE"),
            pips_options::getIntParameter("PRESOLVE_TRACK_ROW_INDEX"),
            pips_options::getBoolParameter("PRESOLVE_TRACK_ROW_LINKING"),
            (pips_options::getIntParameter("PRESOLVE_TRACK_ROW_SYSTEM") == 0 ? EQUALITY_SYSTEM : INEQUALITY_SYSTEM) ),
      tracked_col(COL, pips_options::getIntParameter("PRESOLVE_TRACK_COL_NODE"),
            pips_options::getIntParameter("PRESOLVE_TRACK_COL_INDEX") ),
      objOffset(0.0), obj_offset_chgs(0.0),
      objective_vec_chgs( new SimpleVector(nnzs_col->vec->length()) ),
      lower_bound_implied_by_system(dynamic_cast<StochVectorBase<int>*>(nnzs_col->clone())),
      lower_bound_implied_by_row(dynamic_cast<StochVectorBase<int>*>(nnzs_col->clone())),
      lower_bound_implied_by_node(dynamic_cast<StochVectorBase<int>*>(nnzs_col->clone())),
      upper_bound_implied_by_system(dynamic_cast<StochVectorBase<int>*>(nnzs_col->clone())),
      upper_bound_implied_by_row(dynamic_cast<StochVectorBase<int>*>(nnzs_col->clone())),
      upper_bound_implied_by_node(dynamic_cast<StochVectorBase<int>*>(nnzs_col->clone())),
      absmin_col(dynamic_cast<StochVector*>(sorigprob->g->clone())),
      absmax_col(dynamic_cast<StochVector*>(sorigprob->g->clone())),
      in_bound_tightening(false),
      store_linking_row_boundTightening_A(nnzs_row_A->vecl->length(), 0),
      store_linking_row_boundTightening_C(nnzs_row_C->vecl->length(), 0)
{
   std::memset(array_outdated_indicators, 0, length_array_outdated_indicators * sizeof(bool) );
   outdated_activities = true;

   lower_bound_implied_by_system->setToConstant(-10);
   lower_bound_implied_by_row->setToConstant(-10);
   lower_bound_implied_by_node->setToConstant(-10);
   upper_bound_implied_by_system->setToConstant(-10);
   upper_bound_implied_by_row->setToConstant(-10);
   upper_bound_implied_by_node->setToConstant(-10);

   absmin_col->setToZero();
   absmax_col->setToZero();
   objective_vec_chgs->setToZero();

   presProb = sorigprob->cloneFull(true);

   const int n_linking_vars = (nnzs_col->vec) ? nnzs_col->vec->n : 0;

   const int n_linking_A = (nnzs_row_A->vecl) ? nnzs_row_A->vecl->n : 0;
   const int n_linking_C = (nnzs_row_C->vecl) ? nnzs_row_C->vecl->n : 0;

   assert( n_linking_vars + 2 * n_linking_A + 2 * n_linking_C < std::numeric_limits<int>::max() );

   length_array_nnz_chgs = n_linking_vars + n_linking_A + n_linking_C;
   array_nnz_chgs = new int[length_array_nnz_chgs];
   std::memset(array_nnz_chgs, 0, length_array_nnz_chgs * sizeof(int));

   nnzs_col_chgs = SimpleVectorBaseHandle<int>(new SimpleVectorBase<int>(array_nnz_chgs, n_linking_vars));
   nnzs_row_A_chgs = SimpleVectorBaseHandle<int>(new SimpleVectorBase<int>(array_nnz_chgs + n_linking_vars, n_linking_A));
   nnzs_row_C_chgs = SimpleVectorBaseHandle<int>(new SimpleVectorBase<int>(array_nnz_chgs + n_linking_vars + n_linking_A, n_linking_C));

   lenght_array_act_chgs = n_linking_A * 2 + n_linking_C * 2;
   array_act_chgs = new double[lenght_array_act_chgs];
   std::memset(array_act_chgs, 0.0, lenght_array_act_chgs * sizeof(double));

   actmax_eq_chgs = SimpleVectorHandle( new SimpleVector(array_act_chgs, n_linking_A));
   actmin_eq_chgs = SimpleVectorHandle( new SimpleVector(array_act_chgs + n_linking_A, n_linking_A));
   actmax_ineq_chgs = SimpleVectorHandle( new SimpleVector(array_act_chgs + 2 * n_linking_A, n_linking_C));
   actmin_ineq_chgs = SimpleVectorHandle( new SimpleVector(array_act_chgs + 2 * n_linking_A + n_linking_C, n_linking_C));

   array_act_unbounded_chgs = new int[lenght_array_act_chgs];
   std::memset(array_act_unbounded_chgs, 0, lenght_array_act_chgs * sizeof(int));

   actmax_eq_ubndd_chgs = SimpleVectorBaseHandle<int>( new SimpleVectorBase<int>(array_act_unbounded_chgs, n_linking_A));
   actmin_eq_ubndd_chgs = SimpleVectorBaseHandle<int>( new SimpleVectorBase<int>(array_act_unbounded_chgs + n_linking_A, n_linking_A));
   actmax_ineq_ubndd_chgs = SimpleVectorBaseHandle<int>( new SimpleVectorBase<int>(array_act_unbounded_chgs + 2 * n_linking_A, n_linking_C));
   actmin_ineq_ubndd_chgs = SimpleVectorBaseHandle<int>( new SimpleVectorBase<int>(array_act_unbounded_chgs + 2 * n_linking_A + n_linking_C, n_linking_C));

   lenght_array_bound_chgs = n_linking_A + n_linking_C;
   array_bound_chgs = new double[lenght_array_bound_chgs];
   std::memset(array_bound_chgs, 0.0, lenght_array_bound_chgs * sizeof(double));

   bound_chgs_A = SimpleVectorHandle( new SimpleVector(array_bound_chgs, n_linking_A));
   bound_chgs_C = SimpleVectorHandle( new SimpleVector(array_bound_chgs + n_linking_A, n_linking_C));

   // initialize all dynamic transposed sub matrices
   getSystemMatrix(EQUALITY_SYSTEM).initTransposed(true);
   getSystemMatrix(INEQUALITY_SYSTEM).initTransposed(true);

   recomputeActivities();
   initNnzCounter( *nnzs_row_A, *nnzs_row_C, *nnzs_col);

   initSingletons();
   setUndefinedVarboundsTo(std::numeric_limits<double>::infinity());
   setUndefinedRowboundsTo(std::numeric_limits<double>::infinity());

   initAbsminAbsmaxInCols(*absmin_col, *absmax_col);

   if( iTrackColumn() )
   {
      double xupp, xlow;
      getColBounds(tracked_col, xlow, xupp);

      std::cout << "TRACKING_COLUMN: " << tracked_col << std::endl;
      std::cout << "\tbound x € [" << xlow << ", " << xupp << "]" << std::endl;
   }

   if( iTrackRow() )
   {
      std::cout << "TRACKING_ROW: " << tracked_row << std::endl;
      writeRowLocalToStreamDense(std::cout, tracked_row);
   }
}

PresolveData::~PresolveData()
{
   delete[] array_outdated_indicators;
   delete[] array_nnz_chgs;
   delete[] array_act_chgs;
   delete[] array_bound_chgs;
   delete[] array_act_unbounded_chgs;
}

/* set non existent bounds on all variables to +/- value */
void PresolveData::setUndefinedVarboundsTo(double value)
{
   StochVector& xlow = dynamic_cast<StochVector&>(*presProb->blx);
   StochVector& ixlow = dynamic_cast<StochVector&>(*presProb->ixlow);
   StochVector& xupp = dynamic_cast<StochVector&>(*presProb->bux);
   StochVector& ixupp = dynamic_cast<StochVector&>(*presProb->ixupp);

   setNotIndicatedEntriesTo(xlow, ixlow, -value);
   setNotIndicatedEntriesTo(xupp, ixupp, value);
}

void PresolveData::setUndefinedRowboundsTo(double value)
{
   StochVector& clow = dynamic_cast<StochVector&>(*presProb->bl);
   StochVector& iclow = dynamic_cast<StochVector&>(*presProb->iclow);
   StochVector& cupp = dynamic_cast<StochVector&>(*presProb->bu);
   StochVector& icupp = dynamic_cast<StochVector&>(*presProb->icupp);

   setNotIndicatedEntriesTo(clow, iclow, -value);
   setNotIndicatedEntriesTo(cupp, icupp, value);
}

// TODO: move to StochVector
void PresolveData::setNotIndicatedEntriesTo(StochVector& svec, StochVector& sivec, double value)
{
   assert(svec.children.size() == sivec.children.size());
   assert(svec.children.size() == static_cast<unsigned int>(nChildren));

   for( int node = -1; node < nChildren; ++node )
   {
      if( !nodeIsDummy(node) )
      {
         SimpleVector& vec = getSimpleVecFromColStochVec(svec, node);
         const SimpleVector& ivec = getSimpleVecFromColStochVec(sivec, node);

         assert(vec.n == ivec.n);
         for( int row = 0; row < vec.n; ++row )
         {
            if( PIPSisZero( ivec[row] ) )
               vec[row] = value;
         }
      }
   }

   if( sivec.vecl )
   {
      assert( svec.vecl );
      SimpleVector& vec = getSimpleVecFromRowStochVec(svec, -1, true );
      const SimpleVector& ivec = getSimpleVecFromRowStochVec(sivec, -1, true);

      for( int row = 0; row < vec.n; ++row )
      {
         if (PIPSisZero( ivec[row] ) )
            vec[row] = value;
      }
   }
}

void PresolveData::initAbsminAbsmaxInCols(StochVector& absmin, StochVector& absmax) const
{
   getSystemMatrix(EQUALITY_SYSTEM).getColMinMaxVec(true, true, nullptr, absmin);
   getSystemMatrix(INEQUALITY_SYSTEM).getColMinMaxVec(true, false, nullptr, absmin);

   getSystemMatrix(EQUALITY_SYSTEM).getColMinMaxVec(false, true, nullptr, absmax);
   getSystemMatrix(INEQUALITY_SYSTEM).getColMinMaxVec(false, false, nullptr, absmax);
}

sData* PresolveData::finalize()
{

#ifndef NDEBUG
   if(distributed)
      PIPS_MPIlogicOrArrayInPlace(array_outdated_indicators, length_array_outdated_indicators);
   assert(!outdated_activities && !outdated_lhsrhs && !outdated_nnzs && !outdated_linking_var_bounds && !outdated_obj_vector && !postsolve_linking_row_propagation_needed);
#endif

   /* theoretically it should not matter but there is an assert later which needs all these to be zero */
   setUndefinedVarboundsTo(0.0);
   setUndefinedRowboundsTo(0.0);

   // this removes all columns and rows that are now empty from the problem
   presProb->cleanUpPresolvedData(*nnzs_row_A, *nnzs_row_C, *nnzs_col);

   getSystemMatrix(EQUALITY_SYSTEM).deleteTransposed();
   getSystemMatrix(INEQUALITY_SYSTEM).deleteTransposed();

   return presProb;
}

int PresolveData::getNnzsRow(const INDEX& row) const
{
   if( row.inEqSys() )
   {
      if( row.isLinkingRow() )
         return getSimpleVecFromRowStochVec(*nnzs_row_A, row) + (*nnzs_row_A_chgs)[row.getIndex()];
      else
         return getSimpleVecFromRowStochVec(*nnzs_row_A, row);
   }
   else
   {
      if( row.isLinkingRow() )
         return getSimpleVecFromRowStochVec(*nnzs_row_C, row) + (*nnzs_row_C_chgs)[row.getIndex()];
      else
         return getSimpleVecFromRowStochVec(*nnzs_row_C, row);
   }
}

int PresolveData::getNnzsCol(const INDEX& col) const
{
   assert( col.isCol() );
   if( col.isLinkingCol() )
      return getSimpleVecFromColStochVec(*nnzs_col, col) + (*nnzs_col_chgs)[col.getIndex()];
   else
      return getSimpleVecFromColStochVec(*nnzs_col, col);
}

bool PresolveData::wasColumnRemoved( const INDEX& col ) const
{
   if( postsolver )
    return postsolver->wasColumnRemoved( col );
   else
      return false;
}

void PresolveData::markColumnRemoved( const INDEX& col )
{
   assert( col.isCol() );
   assert( wasColumnRemoved(col) );

   getSimpleVecFromColStochVec(*presProb->g, col) = 0.0;
   getSimpleVecFromColStochVec(*presProb->ixlow, col) = 0.0;
   getSimpleVecFromColStochVec(*presProb->ixupp, col) = 0.0;
   getSimpleVecFromColStochVec(*presProb->blx, col) = 0.0;
   getSimpleVecFromColStochVec(*presProb->bux, col) = 0.0;
}

bool PresolveData::wasRowRemoved( const INDEX& row ) const
{
   if( postsolver )
      return postsolver->wasRowRemoved( row );
   else
      return false;
}

void PresolveData::recomputeActivities(bool linking_only)
{
   recomputeActivities(linking_only, *actmax_eq_part, *actmin_eq_part, *actmax_eq_ubndd, *actmin_eq_ubndd, *actmax_ineq_part,
      *actmin_ineq_part, *actmax_ineq_ubndd, *actmin_ineq_ubndd);

   actmax_eq_chgs->setToZero();
   actmin_eq_chgs->setToZero();
   actmax_ineq_chgs->setToZero();
   actmin_ineq_chgs->setToZero();

   actmax_eq_ubndd_chgs->setToZero();
   actmin_eq_ubndd_chgs->setToZero();
   actmax_ineq_ubndd_chgs->setToZero();
   actmin_ineq_ubndd_chgs->setToZero();

   outdated_activities = false;
}

/** Recomputes the activities of all rows the process knows about. If linking_only is set to true only the linking_rows will get recomputed.
 *  Careful, recomputing linking rows requires MPI communication. Ideally all activities only have to be computed once, when creating the
 *  PresolveData.
 *  After that changes in the activities of linking rows will get stored in the SimpleVectors
 */
void PresolveData::recomputeActivities(bool linking_only, StochVector& actmax_eq_part, StochVector& actmin_eq_part, StochVectorBase<int>& actmax_eq_ubndd,
   StochVectorBase<int>& actmin_eq_ubndd, StochVector& actmax_ineq_part, StochVector& actmin_ineq_part,
   StochVectorBase<int>& actmax_ineq_ubndd, StochVectorBase<int>& actmin_ineq_ubndd) const
{
   const StochGenMatrix& mat_A = getSystemMatrix(EQUALITY_SYSTEM);
   const StochGenMatrix& mat_C = getSystemMatrix(INEQUALITY_SYSTEM);

   const StochVector& xupp = dynamic_cast<StochVector&>(*presProb->bux);
   const StochVector& ixupp = dynamic_cast<StochVector&>(*presProb->ixupp);
   const StochVector& xlow = dynamic_cast<StochVector&>(*presProb->blx);
   const StochVector& ixlow = dynamic_cast<StochVector&>(*presProb->ixlow);

   /* reset vectors keeping track of activities */
   if(!linking_only)
   {
      actmin_eq_part.setToZero();
      actmax_eq_part.setToZero();
      actmin_ineq_part.setToZero();
      actmax_ineq_part.setToZero();

      actmin_eq_ubndd.setToZero();
      actmax_eq_ubndd.setToZero();
      actmin_ineq_ubndd.setToZero();
      actmax_ineq_ubndd.setToZero();
   }
   else
   {
      actmin_eq_part.vecl->setToZero();
      actmax_eq_part.vecl->setToZero();
      actmin_ineq_part.vecl->setToZero();
      actmax_ineq_part.vecl->setToZero();

      actmin_eq_ubndd.vecl->setToZero();
      actmax_eq_ubndd.vecl->setToZero();
      actmin_ineq_ubndd.vecl->setToZero();
      actmax_ineq_ubndd.vecl->setToZero();
   }

   /* compute activities at root node */
   const SimpleVector& xupp_root = dynamic_cast<const SimpleVector&>(*xupp.vec);
   const SimpleVector& ixupp_root = dynamic_cast<const SimpleVector&>(*ixupp.vec);
   const SimpleVector& xlow_root = dynamic_cast<const SimpleVector&>(*xlow.vec);
   const SimpleVector& ixlow_root = dynamic_cast<const SimpleVector&>(*ixlow.vec);

   /* A0/B0 */
   if(!linking_only)
   {
      SimpleVector& actmin_eq_root_part = dynamic_cast<SimpleVector&>(*actmin_eq_part.vec);
      SimpleVector& actmax_eq_root_part = dynamic_cast<SimpleVector&>(*actmax_eq_part.vec);
      SimpleVector& actmin_ineq_root_part = dynamic_cast<SimpleVector&>(*actmin_ineq_part.vec);
      SimpleVector& actmax_ineq_root_part = dynamic_cast<SimpleVector&>(*actmax_ineq_part.vec);

      SimpleVectorBase<int>& actmin_eq_root_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmin_eq_ubndd.vec);
      SimpleVectorBase<int>& actmax_eq_root_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmax_eq_ubndd.vec);
      SimpleVectorBase<int>& actmin_ineq_root_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmin_ineq_ubndd.vec);
      SimpleVectorBase<int>& actmax_ineq_root_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmax_ineq_ubndd.vec);

      addActivityOfBlock(mat_A.Bmat->getStorageDynamicRef(), actmin_eq_root_part, actmin_eq_root_ubndd, actmax_eq_root_part, actmax_eq_root_ubndd,
            xlow_root, ixlow_root, xupp_root, ixupp_root);

      addActivityOfBlock(mat_C.Bmat->getStorageDynamicRef(), actmin_ineq_root_part, actmin_ineq_root_ubndd, actmax_ineq_root_part, actmax_ineq_root_ubndd,
            xlow_root, ixlow_root, xupp_root, ixupp_root);
   }

   SimpleVector& actmin_eq_link_part = dynamic_cast<SimpleVector&>(*actmin_eq_part.vecl);
   SimpleVector& actmax_eq_link_part = dynamic_cast<SimpleVector&>(*actmax_eq_part.vecl);
   SimpleVector& actmin_ineq_link_part = dynamic_cast<SimpleVector&>(*actmin_ineq_part.vecl);
   SimpleVector& actmax_ineq_link_part = dynamic_cast<SimpleVector&>(*actmax_ineq_part.vecl);

   SimpleVectorBase<int>& actmin_eq_link_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmin_eq_ubndd.vecl);
   SimpleVectorBase<int>& actmax_eq_link_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmax_eq_ubndd.vecl);
   SimpleVectorBase<int>& actmin_ineq_link_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmin_ineq_ubndd.vecl);
   SimpleVectorBase<int>& actmax_ineq_link_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmax_ineq_ubndd.vecl);

   /* Bl0 */
   if( my_rank == 0 )
   {
      addActivityOfBlock(mat_A.Blmat->getStorageDynamicRef(), actmin_eq_link_part, actmin_eq_link_ubndd, actmax_eq_link_part, actmax_eq_link_ubndd,
            xlow_root, ixlow_root, xupp_root, ixupp_root);

      addActivityOfBlock(mat_C.Blmat->getStorageDynamicRef(), actmin_ineq_link_part, actmin_ineq_link_ubndd, actmax_ineq_link_part, actmax_ineq_link_ubndd,
            xlow_root, ixlow_root, xupp_root, ixupp_root);
   }

   /* child nodes */
   for(int node = 0; node < nChildren; ++node)
   {
      if( mat_A.children[node]->isKindOf(kStochGenDummyMatrix) && mat_C.children[node]->isKindOf(kStochGenDummyMatrix) )
         continue;

      const SimpleVector& xupp_child = dynamic_cast<const SimpleVector&>(*xupp.children[node]->vec);
      const SimpleVector& ixupp_child = dynamic_cast<const SimpleVector&>(*ixupp.children[node]->vec);
      const SimpleVector& xlow_child = dynamic_cast<const SimpleVector&>(*xlow.children[node]->vec);
      const SimpleVector& ixlow_child = dynamic_cast<const SimpleVector&>(*ixlow.children[node]->vec);

      /* EQUALITY_SYSTEM */
      if( !mat_A.children[node]->isKindOf(kStochGenDummyMatrix) )
      {
         if( !linking_only )
         {
            SimpleVector& actmin_eq_child_part = dynamic_cast<SimpleVector&>(*actmin_eq_part.children[node]->vec);
            SimpleVector& actmax_eq_child_part = dynamic_cast<SimpleVector&>(*actmax_eq_part.children[node]->vec);

            SimpleVectorBase<int>& actmin_eq_child_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmin_eq_ubndd.children[node]->vec);
            SimpleVectorBase<int>& actmax_eq_child_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmax_eq_ubndd.children[node]->vec);

            /* Ai */
            addActivityOfBlock(mat_A.children[node]->Amat->getStorageDynamicRef(), actmin_eq_child_part, actmin_eq_child_ubndd, actmax_eq_child_part,
                  actmax_eq_child_ubndd, xlow_root, ixlow_root, xupp_root, ixupp_root);

            /* Bi */
            addActivityOfBlock(mat_A.children[node]->Bmat->getStorageDynamicRef(), actmin_eq_child_part, actmin_eq_child_ubndd,
                  actmax_eq_child_part, actmax_eq_child_ubndd, xlow_child, ixlow_child, xupp_child, ixupp_child);
         }
         /* Bli */
         addActivityOfBlock(mat_A.children[node]->Blmat->getStorageDynamicRef(), actmin_eq_link_part, actmin_eq_link_ubndd, actmax_eq_link_part,
               actmax_eq_link_ubndd, xlow_child, ixlow_child, xupp_child, ixupp_child);

      }
      /* INEQUALITY_SYSTEM */
      if( !mat_C.children[node]->isKindOf(kStochGenDummyMatrix) )
      {
         if( !linking_only )
         {
            SimpleVector& actmin_ineq_child_part = dynamic_cast<SimpleVector&>(*actmin_ineq_part.children[node]->vec);
            SimpleVector& actmax_ineq_child_part = dynamic_cast<SimpleVector&>(*actmax_ineq_part.children[node]->vec);

            SimpleVectorBase<int>& actmin_ineq_child_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmin_ineq_ubndd.children[node]->vec);
            SimpleVectorBase<int>& actmax_ineq_child_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmax_ineq_ubndd.children[node]->vec);

            /* Ai */
            addActivityOfBlock(mat_C.children[node]->Amat->getStorageDynamicRef(), actmin_ineq_child_part, actmin_ineq_child_ubndd, actmax_ineq_child_part,
                  actmax_ineq_child_ubndd, xlow_root, ixlow_root, xupp_root, ixupp_root);

            /* Bi */
            addActivityOfBlock(mat_C.children[node]->Bmat->getStorageDynamicRef(), actmin_ineq_child_part, actmin_ineq_child_ubndd, actmax_ineq_child_part,
                  actmax_ineq_child_ubndd, xlow_child, ixlow_child, xupp_child, ixupp_child);
         }

         /* Bli */
         addActivityOfBlock(mat_C.children[node]->Blmat->getStorageDynamicRef(), actmin_ineq_link_part, actmin_ineq_link_ubndd,
               actmax_ineq_link_part, actmax_ineq_link_ubndd, xlow_child, ixlow_child, xupp_child, ixupp_child);
      }

   }

   /* allreduce linking constraint activities */
   if( distributed )
   {
      // todo is copying and then allreducing once cheaper than allreducing 4 times ? by a lot?
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVector*>(actmin_eq_part.vecl)->elements(), actmin_eq_part.vecl->n);
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVector*>(actmax_eq_part.vecl)->elements(), actmax_eq_part.vecl->n);
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVector*>(actmin_ineq_part.vecl)->elements(), actmin_ineq_part.vecl->n);
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVector*>(actmax_ineq_part.vecl)->elements(), actmax_ineq_part.vecl->n);

      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVectorBase<int>*>(actmin_eq_ubndd.vecl)->elements(), actmin_eq_ubndd.vecl->n);
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVectorBase<int>*>(actmax_eq_ubndd.vecl)->elements(), actmax_eq_ubndd.vecl->n);
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVectorBase<int>*>(actmin_ineq_ubndd.vecl)->elements(), actmin_ineq_ubndd.vecl->n);
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVectorBase<int>*>(actmax_ineq_ubndd.vecl)->elements(), actmax_ineq_ubndd.vecl->n);
   }


   /* set activities to infinity // theoretically not necessary but for debugging */

   /* non-linking */
   if( !linking_only )
   {
      const bool linking = false;
      for(int node = -1; node < nChildren; ++node)
      {
         if(!nodeIsDummy(node))
         {
            /* EQUALITY_SYSTEM */
            SimpleVector& actmin_eq = getSimpleVecFromRowStochVec(actmin_eq_part, node, linking);
            SimpleVector& actmax_eq = getSimpleVecFromRowStochVec(actmax_eq_part, node, linking);
            const SimpleVectorBase<int>& actmin_ubndd_eq = getSimpleVecFromRowStochVec(actmin_eq_ubndd, node, linking);
            const SimpleVectorBase<int>& actmax_ubndd_eq = getSimpleVecFromRowStochVec(actmax_eq_ubndd, node, linking);

            for(int row = 0; row < actmin_eq.n; ++row)
            {
               if(actmin_ubndd_eq[row] >= 2)
                  actmin_eq[row] = INF_NEG;
               if(actmax_ubndd_eq[row] >= 2)
                  actmax_eq[row] = INF_POS;
            }

            /* INEQUALITY_SYSTEM */
            SimpleVector& actmin_ineq = getSimpleVecFromRowStochVec(actmin_ineq_part, node, linking);
            SimpleVector& actmax_ineq = getSimpleVecFromRowStochVec(actmax_ineq_part, node, linking);
            const SimpleVectorBase<int>& actmin_ubndd_ineq = getSimpleVecFromRowStochVec(actmin_ineq_ubndd, node, linking);
            const SimpleVectorBase<int>& actmax_ubndd_ineq = getSimpleVecFromRowStochVec(actmax_ineq_ubndd, node, linking);

            for(int row = 0; row < actmin_ineq.n; ++row)
            {
               if(actmin_ubndd_ineq[row] >= 2)
                  actmin_ineq[row] = INF_NEG;
               if(actmax_ubndd_ineq[row] >= 2)
                  actmax_ineq[row] = INF_POS;
            }
         }
      }
   }

   for(int row = 0; row < actmin_eq_part.vecl->n; ++row)
   {
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmin_eq_ubndd.vecl)[row] >= 2 )
         dynamic_cast<SimpleVector&>(*actmin_eq_part.vecl)[row] = INF_NEG;
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmax_eq_ubndd.vecl)[row] >= 2 )
         dynamic_cast<SimpleVector&>(*actmax_eq_part.vecl)[row] = INF_POS;
   }
   for(int row = 0; row < actmin_ineq_part.vecl->n; ++row)
   {
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmin_ineq_ubndd.vecl)[row] >= 2 )
         dynamic_cast<SimpleVector&>(*actmin_ineq_part.vecl)[row] = INF_NEG;
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmax_ineq_ubndd.vecl)[row] >= 2 )
         dynamic_cast<SimpleVector&>(*actmax_ineq_part.vecl)[row] = INF_POS;
   }

   if( iTrackRow() )
   {
      double act_min, act_max;
      int act_min_ubndd, act_max_ubndd;
      getRowActivities(tracked_row, act_max, act_min, act_max_ubndd, act_min_ubndd);

      std::cout << "TRACKING_ROW: Recomputed activity of tracked row" << std::endl;
      std::cout << "\tnew min/max activity is: " << act_min << "/" << act_max << ", min/max unbounded counters are " << act_min_ubndd << "/" << act_max_ubndd << std::endl;
   }
}

/** Computes minimal and maximal activity of all rows in given matrix. Adds activities to min/max_activities accordingly. */
void PresolveData::addActivityOfBlock( const SparseStorageDynamic& matrix, SimpleVector& min_partact, SimpleVectorBase<int>& unbounded_min,
   SimpleVector& max_partact, SimpleVectorBase<int>& unbounded_max, const SimpleVector& xlow, const SimpleVector& ixlow,
   const SimpleVector& xupp, const SimpleVector& ixupp) const
{
   assert( xlow.n == matrix.getN() && ixlow.n == matrix.getN() && xupp.n == matrix.getN() && ixupp.n == matrix.getN());
   assert( max_partact.n == matrix.getM() && min_partact.n == matrix.getM());
   assert( unbounded_min.n == matrix.getM() && unbounded_max.n == matrix.getM());

   for( int row = 0; row < matrix.getM(); ++row)
   {
      for( int j = matrix.getRowPtr(row).start; j < matrix.getRowPtr(row).end; j++)
      {
         const int col = matrix.getJcolM(j);
         const double entry = matrix.getMat(j);

         assert( !PIPSisZero(entry) );

         if( entry > 0)
         {
            if( !PIPSisZero(ixlow[col]) )
               min_partact[row] += entry * xlow[col];
            else
               ++unbounded_min[row];
            if( !PIPSisZero(ixupp[col]) )
               max_partact[row] += entry * xupp[col];
            else
               ++unbounded_max[row];
         }
         else
         {
            if( !PIPSisZero(ixupp[col]) )
               min_partact[row] += entry * xupp[col];
            else
               ++unbounded_min[row];
            if( !PIPSisZero(ixlow[col]) )
               max_partact[row] += entry * xlow[col];
            else
               ++unbounded_max[row];
         }
      }
      if(unbounded_max[row] >= 2)
         max_partact[row] = std::numeric_limits<double>::infinity();
      if(unbounded_min[row] >= 2)
         min_partact[row] = -std::numeric_limits<double>::infinity();
   }
}

void PresolveData::allreduceLinkingVarBounds()
{
   PIPS_MPIgetLogicOrInPlace( outdated_linking_var_bounds );

   if(!outdated_linking_var_bounds)
      return;

   if(distributed)
   {
      SimpleVector& xlow_new_vec = getSimpleVecFromColStochVec(*presProb->blx, -1);
      SimpleVector& xupp_new_vec = getSimpleVecFromColStochVec(*presProb->bux, -1);
      SimpleVector& ixlow_new_vec = getSimpleVecFromColStochVec(*presProb->ixlow, -1);
      SimpleVector& ixupp_new_vec = getSimpleVecFromColStochVec(*presProb->ixupp, -1);

      /* copy old values for later compairson */
      SimpleVector* xlow_old_vec = dynamic_cast<SimpleVector*>(xlow_new_vec.cloneFull());
      SimpleVector* xupp_old_vec = dynamic_cast<SimpleVector*>(xupp_new_vec.cloneFull());
      SimpleVector* ixlow_old_vec = dynamic_cast<SimpleVector*>(ixlow_new_vec.cloneFull());
      SimpleVector* ixupp_old_vec = dynamic_cast<SimpleVector*>(ixupp_new_vec.cloneFull());

      PIPS_MPImaxArrayInPlace(xlow_new_vec.elements(), xlow_new_vec.length());
      PIPS_MPImaxArrayInPlace(ixlow_new_vec.elements(), ixlow_new_vec.length());
      PIPS_MPIminArrayInPlace(xupp_new_vec.elements(), xupp_new_vec.length());
      PIPS_MPImaxArrayInPlace(ixupp_new_vec.elements(), ixupp_new_vec.length());

      // this will affect the activities of basically all rows - use with care
      for(int col = 0; col < xlow_new_vec.length(); ++col)
      {
         double xupp_old = INF_POS;
         double xlow_old = INF_NEG;
         double xupp_new = INF_POS;
         double xlow_new = INF_NEG;

         if( PIPSisEQ( (*ixupp_old_vec)[col], 1.0) )
            xupp_old = (*xupp_old_vec)[col];
         if( PIPSisEQ( (*ixlow_old_vec)[col], 1.0) )
            xlow_old = (*xlow_old_vec)[col];
         if( PIPSisEQ( ixupp_new_vec[col], 1.0) )
            xupp_new = xupp_new_vec[col];
         if( PIPSisEQ( ixlow_new_vec[col], 1.0) )
            xlow_new = xlow_new_vec[col];

         const INDEX col_idx(COL, -1, col);

         if(xupp_new != INF_POS || xlow_new != INF_NEG)
            updateRowActivities(col_idx, xlow_new, xupp_new, xlow_old, xupp_old);
      }
      
      delete xlow_old_vec;
      delete xupp_old_vec;
      delete ixlow_old_vec;
      delete ixupp_old_vec;
   }
   outdated_linking_var_bounds = false;
}

/** allreduce changes in the activities of the linking rows and update the linking row activities */
void PresolveData::allreduceAndApplyLinkingRowActivities()
{
   /* check for new linking rows that have to be computed */
   PIPS_MPIgetLogicOrInPlace(outdated_activities);
   // todo : better criterion - if there have been this and that many changes etc
   PIPS_MPIgetSumInPlace(linking_rows_need_act_computation);

   if(!outdated_activities && linking_rows_need_act_computation == 0)
      return;

   /* allreduce unbounded changes, compute rows that now have at most one unbounded entry, strore the local
    * rows in the changes array, allreduces MPI_SUM the changes array and update local activities with the global ones
    */
   if(distributed)
      PIPS_MPIsumArrayInPlace(array_act_unbounded_chgs, lenght_array_act_chgs);

   dynamic_cast<SimpleVectorBase<int>&>(*actmin_eq_ubndd->vecl).axpy( 1.0, *actmin_eq_ubndd_chgs);
   dynamic_cast<SimpleVectorBase<int>&>(*actmax_eq_ubndd->vecl).axpy( 1.0, *actmax_eq_ubndd_chgs);
   dynamic_cast<SimpleVectorBase<int>&>(*actmin_ineq_ubndd->vecl).axpy( 1.0, *actmin_ineq_ubndd_chgs);
   dynamic_cast<SimpleVectorBase<int>&>(*actmax_ineq_ubndd->vecl).axpy( 1.0, *actmax_ineq_ubndd_chgs);

   /* equality system */
   for(int row = 0; row < actmin_eq_ubndd->vecl->n; ++row)
   {
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmin_eq_ubndd->vecl)[row] < 2
            && dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl)[row] == INF_NEG)
      {
         (*actmin_eq_chgs)[row] = computeLocalLinkingRowMinOrMaxActivity( INDEX(ROW, -1, row, true, EQUALITY_SYSTEM), false);
         dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl)[row] = 0;
      }

      if( dynamic_cast<SimpleVectorBase<int>&>(*actmax_eq_ubndd->vecl)[row] < 2
            && dynamic_cast<SimpleVector&>(*actmax_eq_part->vecl)[row] == INF_POS)
      {
         (*actmax_eq_chgs)[row] = computeLocalLinkingRowMinOrMaxActivity( INDEX(ROW, -1, row, true, EQUALITY_SYSTEM), true);
         dynamic_cast<SimpleVector&>(*actmax_eq_part->vecl)[row] = 0;
      }
   }

   /* inequality system */
   for(int row = 0; row < actmin_ineq_ubndd->vecl->n; ++row)
   {
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmin_ineq_ubndd->vecl)[row] < 2
            && dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl)[row] == INF_NEG)
      {
         (*actmin_ineq_chgs)[row] = computeLocalLinkingRowMinOrMaxActivity( INDEX(ROW, -1, row, true, INEQUALITY_SYSTEM), false);
         dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl)[row] = 0;
      }

      if( dynamic_cast<SimpleVectorBase<int>&>(*actmax_ineq_ubndd->vecl)[row] < 2
            && dynamic_cast<SimpleVector&>(*actmax_ineq_part->vecl)[row] == INF_POS)
      {
         (*actmax_ineq_chgs)[row] = computeLocalLinkingRowMinOrMaxActivity( INDEX(ROW, -1, row, true, INEQUALITY_SYSTEM), true);
         dynamic_cast<SimpleVector&>(*actmax_ineq_part->vecl)[row] = 0;
      }
   }

   if( distributed )
      PIPS_MPIsumArrayInPlace(array_act_chgs, lenght_array_act_chgs);

   dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl).axpy(1.0, *actmin_eq_chgs);
   dynamic_cast<SimpleVector&>(*actmax_eq_part->vecl).axpy(1.0, *actmax_eq_chgs);
   dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl).axpy(1.0, *actmin_ineq_chgs);
   dynamic_cast<SimpleVector&>(*actmax_ineq_part->vecl).axpy(1.0, *actmax_ineq_chgs);

   actmin_eq_chgs->setToZero();
   actmax_eq_chgs->setToZero();
   actmin_ineq_chgs->setToZero();
   actmax_ineq_chgs->setToZero();

   actmin_eq_ubndd_chgs->setToZero();
   actmax_eq_ubndd_chgs->setToZero();
   actmin_ineq_ubndd_chgs->setToZero();
   actmax_ineq_ubndd_chgs->setToZero();

   outdated_activities = false;
   linking_rows_need_act_computation = 0;

#ifndef NDEBUG
   /* check if all rows with unbounded counter <= 1 are computed and the rest is set to +/- infinity */

   for(int i = 0; i < actmin_eq_ubndd->vecl->n; ++i)
   {
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmin_eq_ubndd->vecl)[i] < 2 )
         assert( std::fabs(dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl)[i]) != INF_POS );
      else
         assert( dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl)[i] == INF_NEG );
   }

   for(int i = 0; i < actmin_ineq_ubndd->vecl->n; ++i)
   {
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmin_ineq_ubndd->vecl)[i] < 2 )
         assert( std::fabs(dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl)[i]) != INF_POS );
      else
         assert( dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl)[i] == INF_NEG );
   }
#endif
}

/** allreduce changes in the nnz counters and apply them locally */
void PresolveData::allreduceAndApplyNnzChanges()
{
   PIPS_MPIgetLogicOrInPlace(outdated_nnzs);

   if(!outdated_nnzs)
      return;

   if( distributed )
      PIPS_MPIsumArrayInPlace(array_nnz_chgs, length_array_nnz_chgs);

   /* update local nnzCounters */
   nnzs_col->vec->axpy(1, *nnzs_col_chgs);
   nnzs_row_A->vecl->axpy(1, *nnzs_row_A_chgs);
   nnzs_row_C->vecl->axpy(1, *nnzs_row_C_chgs);

   // todo : this can be done more efficiently, e.g. while substracting
   // todo : this has still flaws - new singleton rows in B0 and D0 are not communicated properly for some reason - might happen else where
   for( int i = 0; i < nnzs_col_chgs->length(); ++i )
   {
      if( (*nnzs_col_chgs)[i] > 0 && dynamic_cast<SimpleVectorBase<int>&>(*nnzs_col->vec)[i] == 1 )
         singleton_cols.push( INDEX(COL, -1, i) );
   }
   
   for( int i = 0; i < nnzs_row_A_chgs->length(); ++i )
   {
      if( (*nnzs_row_A_chgs)[i] > 0 && dynamic_cast<SimpleVectorBase<int>&>(*nnzs_row_A->vec)[i] == 1 )
         singleton_rows.push( INDEX(ROW, -1, i, true, EQUALITY_SYSTEM) );
   }
   
   for( int i = 0; i < nnzs_row_C_chgs->length(); ++i )
   {
      if( (*nnzs_row_C_chgs)[i] > 0 && dynamic_cast<SimpleVectorBase<int>&>(*nnzs_row_C->vec)[i] == 1 )
         singleton_rows.push( INDEX(ROW, -1, i, true, INEQUALITY_SYSTEM) );
   }

#ifndef NDEBUG
   int minval = -1;
   int index = -1;
   nnzs_col->min(minval, index);
   assert( minval >= 0 );
   nnzs_row_A->vecl->min(minval, index);
   assert(minval >= 0);
   nnzs_row_C->vecl->min(minval, index);
   assert(minval >= 0);
#endif

   nnzs_col_chgs->setToZero();
   nnzs_row_A_chgs->setToZero();
   nnzs_row_C_chgs->setToZero();

   outdated_nnzs = false;
}

void PresolveData::allreduceAndApplyBoundChanges()
{
   PIPS_MPIgetLogicOrInPlace(outdated_lhsrhs);

   if(!outdated_lhsrhs)
      return;

   if(distributed)
      PIPS_MPIsumArrayInPlace(array_bound_chgs, lenght_array_bound_chgs);

   dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*presProb->bA).vecl)->axpy( 1.0, *bound_chgs_A);
   dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*presProb->bl).vecl)->axpy( 1.0, *bound_chgs_C);
   dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*presProb->bu).vecl)->axpy( 1.0, *bound_chgs_C);

   bound_chgs_A->setToZero();
   bound_chgs_C->setToZero();

   outdated_lhsrhs = false;
}

void PresolveData::allreduceAndApplyObjVecChanges()
{
   PIPS_MPIgetLogicOrInPlace(outdated_obj_vector);

   if(!outdated_obj_vector)
      return;

   if(distributed)
      PIPS_MPIsumArrayInPlace(objective_vec_chgs->elements(), objective_vec_chgs->length());

   dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*presProb->g).vec)->axpy(1.0, *objective_vec_chgs);

   objective_vec_chgs->setToZero();
   outdated_obj_vector = false;
}

void PresolveData::allreduceObjOffset()
{
   if(distributed)
      PIPS_MPIgetSumInPlace(obj_offset_chgs);

   objOffset += obj_offset_chgs;
   obj_offset_chgs = 0;
}

void PresolveData::initNnzCounter(StochVectorBase<int>& nnzs_row_A, StochVectorBase<int>& nnzs_row_C, StochVectorBase<int>& nnzs_col) const
{
   StochGenMatrix& A = getSystemMatrix(EQUALITY_SYSTEM);
   StochGenMatrix& C = getSystemMatrix(INEQUALITY_SYSTEM);

   StochVectorBaseHandle<int> colClone(dynamic_cast<StochVectorBase<int>*>(nnzs_col.clone()));

   A.getNnzPerRow(nnzs_row_A);
   C.getNnzPerRow(nnzs_row_C);
   A.getNnzPerCol(nnzs_col);
   C.getNnzPerCol(*colClone);

   nnzs_col.axpy(1, *colClone);
}

// todo : for columns also find row - then later no lookup is required
void PresolveData::initSingletons()
{
   /* rows of A */
   /* B0 */
   for(int i = 0; i < nnzs_row_A->vec->n; ++i)
   {
      if( dynamic_cast<SimpleVectorBase<int>&>(*nnzs_row_A->vec)[i] == 1)
         singleton_rows.push( INDEX(ROW, -1, i, false, EQUALITY_SYSTEM) );
   }

   /* Bl0 */
   if(nnzs_row_A->vecl != nullptr)
   {
      for(int i = 0; i < nnzs_row_A->vecl->n; ++i)
      {
         if( dynamic_cast<SimpleVectorBase<int>&>(*nnzs_row_A->vecl)[i] == 1)
            singleton_rows.push( INDEX(ROW, -1, i, true, EQUALITY_SYSTEM) );
      }
   }

   /* children An + Bn */
   for(int i = 0; i < nChildren; ++i)
   {
      if( !nodeIsDummy(i) )
      {
         for(int j = 0; j < nnzs_row_A->children[i]->vec->n; ++j)
         {
            if( dynamic_cast<SimpleVectorBase<int>&>(*nnzs_row_A->children[i]->vec)[j] == 1)
               singleton_rows.push( INDEX(ROW, i, j, false, EQUALITY_SYSTEM) );
         }
      }
   }

   /* rows of C */
   /* B0 */
   for(int i = 0; i < nnzs_row_C->vec->n; ++i)
      if( dynamic_cast<SimpleVectorBase<int>&>(*nnzs_row_C->vec)[i] == 1)
         singleton_rows.push( INDEX(ROW, -1, i, false, INEQUALITY_SYSTEM) );

   /* Bl0 */
   if( nnzs_row_C->vecl != nullptr)
   {
      for( int i = 0; i < nnzs_row_C->vecl->n; ++i )
      {
         if( dynamic_cast<SimpleVectorBase<int>&>(*nnzs_row_C->vecl)[i] == 1 )
            singleton_rows.push( INDEX(ROW, -1, i, true, INEQUALITY_SYSTEM) );
      }
   }

   /* children An + Bn */
   for( int i = 0; i < nChildren; ++i )
   {
      if( !nodeIsDummy(i) )
      {
         for( int j = 0; j < nnzs_row_C->children[i]->vec->n; ++j )
         {
            if( dynamic_cast<SimpleVectorBase<int>&>(*nnzs_row_C->children[i]->vec)[j] == 1 )
               singleton_rows.push( INDEX(ROW, i, j, false, INEQUALITY_SYSTEM) );
         }
      }
   }

   /* columns */
   for(int i = 0; i < nnzs_col->vec->n; ++i)
   {
      if( dynamic_cast<SimpleVectorBase<int>&>(*nnzs_col->vec)[i] == 1)
         singleton_cols.push( INDEX(COL, -1, i) );
   }

   for(int i = 0; i < nChildren; ++i)
   {
      if(nnzs_col->children[i]->vec != nullptr)
      {
         for(int j = 0; j < nnzs_col->children[i]->vec->n; ++j)
         {
            if( dynamic_cast<SimpleVectorBase<int>&>(*nnzs_col->children[i]->vec)[j] == 1)
               singleton_cols.push( INDEX(COL, i, j) );
         }
      }
   }

}

bool PresolveData::reductionsEmpty()
{
   double recv = 0;
   if(distributed)
   {
      PIPS_MPIlogicOrArrayInPlace(array_outdated_indicators, length_array_outdated_indicators);

      recv = PIPS_MPIgetSum(obj_offset_chgs);
   }
   if( outdated_obj_vector && my_rank == 0 )
      std::cout << "Error: objective vector is outdated!" << std::endl;
   if( outdated_activities && my_rank == 0 )
      std::cout << "Error: activities are outdated!" << std::endl;
   if( outdated_lhsrhs && my_rank == 0 )
      std::cout << "Error: lhs_rhs outdated!" << std::endl;
   if( outdated_linking_var_bounds && my_rank == 0 )
      std::cout << "Error: linking var bounds outdated!" << std::endl;
   if( outdated_nnzs && my_rank == 0 )
      std::cout << "Error: non-zeros outdated!" << std::endl;
   if( recv != 0.0 && my_rank == 0 )
      std::cout << "Error: objcetive offset not synced!" << std::endl;

   return !outdated_obj_vector && !outdated_activities && !outdated_lhsrhs && !outdated_linking_var_bounds && !outdated_nnzs && (recv == 0.0) && !postsolve_linking_row_propagation_needed;
}

bool PresolveData::presDataInSync() const
{
   return presProb->isRootNodeInSync() && verifyNnzcounters() && verifyActivities();
}

// todo : if small entry was removed from system no postsolve is necessary - if coefficient was removed because impact of changes in variable are small
// also rhs lhs will be adjusted - this has to be reversed later - there will also be a problem with the reduced costs in than particular row ?
// todo : postsolve if bounds adjusted because of deleted matrix entry simply reverse the adjustment - no changes in multipliers
//- row stays active / inactive ? not true..
void PresolveData::deleteEntryAtIndex( const INDEX& row, const INDEX& col, int col_index )
{
   assert( row.isRow() );
   assert( col.isCol() );
   assert( row.hasValidNode(nChildren) );
   assert( col.hasValidNode(nChildren) );

   const SparseStorageDynamic& storage = getSparseGenMatrix(row, col)->getStorageDynamicRef();
   
   const double val = storage.getMat(col_index);
   assert( col.getIndex() == storage.getJcolM(col_index) );
   assert( storage.getRowPtr(row.getIndex()).start <= col_index && col_index < storage.getRowPtr(row.getIndex()).end );

   const bool at_root = ( row.getNode() == -1 && col.getNode() == -1 );

   if( postsolver )
   {
      postsolver->notifyColModified( col );
      postsolver->notifyRowModified( row );
   }

   getSparseGenMatrix(row, col)->removeEntryAtRowColIndex(row.getIndex(), col_index);

   double& xlow = getSimpleVecFromColStochVec(*presProb->blx, col);

   if( row.isLinkingRow() || col.isLinkingCol() )
      outdated_nnzs = true;

   /* adjust rhs and lhs */
   adjustMatrixRhsLhsBy( row, -val * xlow, at_root );

   /* adjust activity */
   /* deletion of entry acts on activities like fixing it's value to zero */
   adjustRowActivityFromDeletion( row, col, val );

   /* adjust nnz counters */
   reduceNnzCounterRowBy( row, 1, at_root );
   reduceNnzCounterColumnBy( col, 1, at_root);
}

void PresolveData::resetOriginallyFreeVarsBounds(const sData& orig_prob)
{
   // todo : tell postsolver about released variables
   assert( 0 && "not yet properly implemented" );

#ifndef NDEBUG
   if(my_rank == 0)
      std::cout << "Resetting all presolved variable bounds of originally free variables" <<::endl; 
#endif

   long long n = 0;
   for( int node = -1; node < nChildren; ++node )
   {
      n += resetOriginallyFreeVarsBounds( getSimpleVecFromColStochVec(*orig_prob.ixlow, node), getSimpleVecFromColStochVec(*orig_prob.ixupp, node), node);
      if(my_rank != 0 && node == -1)
         n = 0;
   }

#ifndef NDEBUG
   PIPS_MPIgetSumInPlace(n);
   if(my_rank == 0)
      std::cout << "Reset " << n << " bounds" << std::endl;
#endif
}

// todo the -2 thing..
long PresolveData::resetOriginallyFreeVarsBounds(const SimpleVector& ixlow_orig, const SimpleVector& ixupp_orig, int node)
{
   long reset_bounds = 0;

   if( nodeIsDummy( node ) )
      return reset_bounds;

   SimpleVector& ixlow = getSimpleVecFromColStochVec(*presProb->ixlow, node);
   SimpleVector& ixupp = getSimpleVecFromColStochVec(*presProb->ixupp, node);

   SimpleVector& xlow = getSimpleVecFromColStochVec(*presProb->blx, node);
   SimpleVector& xupp = getSimpleVecFromColStochVec(*presProb->bux, node);

   SimpleVectorBase<int>& nnzs_col_vec = getSimpleVecFromColStochVec(*nnzs_col, node);

   /* check whether row that implied bound is still there - if so, and if the variable is still in that row we can remove the bound again
    * since it should still be implied 
    */

   /* todo: actually need to check whether a bound is still an implied one - if so we can reset it - this needs more mechanisms */
   /* store row that implied bound - if row still there - check if bound still implied (or even better bound implied) - if so - reset bound */
   /* print stats on how many bounds were reset */
   for(int col = 0; col < ixlow.n; ++col)
   {
      /* do not reset fixed columns */
      if( nnzs_col_vec[col] != 0 && ( PIPSisZero(ixupp_orig[col]) || PIPSisZero(ixlow_orig[col]) ) )
      {
         int sys_row_lower = getSimpleVecFromColStochVec(*lower_bound_implied_by_system, node)[col];
         int row_lower = getSimpleVecFromColStochVec(*lower_bound_implied_by_row, node)[col];
         int node_row_lower = getSimpleVecFromColStochVec(*lower_bound_implied_by_node, node)[col];
         int sys_row_upper = getSimpleVecFromColStochVec(*upper_bound_implied_by_system, node)[col];
         int row_upper = getSimpleVecFromColStochVec(*upper_bound_implied_by_row, node)[col];
         int node_row_upper = getSimpleVecFromColStochVec(*upper_bound_implied_by_node, node)[col];

         /* do not reset bounds implied by singleton rows since these rows are already removed from the problem */
         if( PIPSisZero(ixupp_orig[col]) && PIPSisEQ(ixupp[col], 1.0) )
         {
            if( !(0 <= sys_row_upper && sys_row_upper <= 1) )
               std::cout << sys_row_upper << std::endl;
            assert( 0 <= sys_row_upper && sys_row_upper <= 1);
            assert(row_upper != -10);
            assert(node_row_upper != -10);

            const StochVectorBase<int>& nnzs = (sys_row_upper == 0) ? *nnzs_row_A : *nnzs_row_C;

            if( getSimpleVecFromRowStochVec(nnzs, (node_row_upper == -2) ? -1 : node_row_upper, 
               (node_row_upper == -2) )[row_upper] != 0 )
            {
               ixupp[col] = 0;
               xupp[col] = 0.0;
               ++reset_bounds;
            }
         }

         if( PIPSisZero(ixlow_orig[col]) && PIPSisEQ(ixlow[col], 1.0) )
         {
            assert( 0 <= sys_row_lower && sys_row_lower <= 1);
            assert(row_lower != -10);
            assert(node_row_lower != -10);

            const StochVectorBase<int>& nnzs = (sys_row_lower == 0) ? *nnzs_row_A : *nnzs_row_C;

            // todo : problably should also check whether variable is still in row - not necessary for so far implemented presolvers
            if( getSimpleVecFromRowStochVec(nnzs, (node_row_lower == -2) ? -1 : node_row_lower, 
               (row_lower == -2) )[row_lower] != 0 )
            {
               ixlow[col] = 0;
               xlow[col] = 0.0;
               ++reset_bounds;
            }
         }
      }
   }
   return reset_bounds;
}

INDEX PresolveData::getRowMarkedAsImplyingColumnBound(const INDEX& col, bool upper_bound)
{
   assert(col.isCol());

   StochVectorBaseHandle<int>& by_node = upper_bound ? upper_bound_implied_by_node : lower_bound_implied_by_node;
   StochVectorBaseHandle<int>& by_system = upper_bound ? upper_bound_implied_by_system : lower_bound_implied_by_system;
   StochVectorBaseHandle<int>& by_row = upper_bound ? upper_bound_implied_by_row : lower_bound_implied_by_row;

   const int row_node_int = getSimpleVecFromColStochVec(*by_node, col);
   const int system_type_int = getSimpleVecFromColStochVec(*by_system, col);
   const int row_index = getSimpleVecFromColStochVec(*by_row, col);

   if( row_node_int == -10 )
   {
      assert( system_type_int == -10 );
      assert( row_index == -10 );
      return INDEX();
   }

   assert( 0 <= row_index );
   assert( -2 <= row_node_int && row_node_int < nChildren );

   const int row_node = ( row_node_int == -2 ) ? -1 : row_node_int;

   assert( 0 <= system_type_int && system_type_int <= 1 );

   const SystemType system_type = ( system_type_int == 0 ) ? EQUALITY_SYSTEM : INEQUALITY_SYSTEM;
   const bool linking = ( row_node_int == -2 ) ? true : false;

   return INDEX(ROW, row_node, row_index, linking, system_type);
}

void PresolveData::markRowAsImplyingColumnBound(const INDEX& col, const INDEX& row, bool upper_bound)
{
   assert( row.isRow() );
   assert( col.isCol() );

   StochVectorBaseHandle<int>& by_node = upper_bound ? upper_bound_implied_by_node : lower_bound_implied_by_node;
   StochVectorBaseHandle<int>& by_system = upper_bound ? upper_bound_implied_by_system : lower_bound_implied_by_system;
   StochVectorBaseHandle<int>& by_row = upper_bound ? upper_bound_implied_by_row : lower_bound_implied_by_row;

   getSimpleVecFromColStochVec(*by_node, col) = (row.getLinking()) ? -2 : row.getNode();
   getSimpleVecFromColStochVec(*by_system, col) = row.getSystemType();
   getSimpleVecFromColStochVec(*by_row, col) = row.getIndex();
}

/** returns whether or not the current bound on column col is implied by row */
bool PresolveData::varBoundImpliedFreeBy( bool upper, const INDEX& col, const INDEX& row)
{
   assert( col.isCol() );
   assert( row.isRow() );

   // todo : not sure whether there is an instance that ever calls this
   // todo : theoretically there might be nnzs changes in some buffers somewhere - should not happen but how to check?
   // todo : should this be only one method varBoundsImpliedFree?
   if( 0 == getNnzsRow(row) )
      return false;

   bool res = false;
   if( upper )
      res = (row == getRowMarkedAsImplyingColumnBound(col, true));
   else
      res = (row == getRowMarkedAsImplyingColumnBound(col, false));

   /* check whether bounds is actually still implied by the row -- also checks whether col is still in that row */
   if( res == true )
   {
      bool upper_implied, lower_implied;
      varboundImpliedFreeFullCheck(upper_implied, lower_implied, col, row);
      if(upper)
         return upper_implied;
      else
         return lower_implied;
   }

   return res;
}

/* uses current activities (non-updated) to check whether said column's bounds are implied by row */
void PresolveData::varboundImpliedFreeFullCheck(bool& upper_implied, bool& lower_implied, const INDEX& col, const INDEX& row) const
{
   assert( !outdated_activities );

   assert( col.isCol() );
   assert( row.isRow() );

   upper_implied = false;
   lower_implied = false;

   /* calculate implied bounds again and check whether the col bounds are actually still implied */
   /* get activities */
   double max_act, min_act;
   int max_ubndd, min_ubndd;
   getRowActivities(row, max_act, min_act, max_ubndd, min_ubndd);

   /* get matrix in order to get the coefficient of col in row */
   const SparseStorageDynamic& mat = getSparseGenMatrix(row, col)->getStorageDynamicRef();

   const int row_start = mat.getRowPtr(row.getIndex()).start;
   const int row_end = mat.getRowPtr(row.getIndex()).end;
   int col_ptr;

   /* find coefficient and column */
   for( col_ptr = row_start; col_ptr < row_end; ++col_ptr)
   {
      int col_mat = mat.getJcolM(col_ptr);
      if( col.getIndex() == col_mat )
         break;
   }

   /* if nothing was found the column has already been removed from that row and nothing is implied anymore */
   if(col_ptr == row_end)
      return;

   /* coefficient of col in row */
   const double coeff = mat.getMat(col_ptr);
   assert(!PIPSisZero(coeff));

   /* current bounds */
   double xupp, xlow;
   getColBounds(col, xlow, xupp);

   if( coeff > 0 )
   {
      min_act -= coeff * xlow;
      max_act -= coeff * xupp;
   }
   else
   {
      min_act -= coeff * xupp;
      max_act -= coeff * xlow;
   }

   const double rhs = row.inEqSys() ? getSimpleVecFromRowStochVec(*presProb->bA, row) :
      getSimpleVecFromRowStochVec(*presProb->bu, row);
   const double lhs = row.inEqSys() ? getSimpleVecFromRowStochVec(*presProb->bA, row) :
      getSimpleVecFromRowStochVec(*presProb->bl, row);

   /* check bound implied by row */

   /* calculate an check implied upper bound */
   if( max_ubndd == 0 && xupp != INF_POS )
   {
      if( 0.0 < coeff )
      {
         assert( row.inEqSys() || !PIPSisZero(getSimpleVecFromRowStochVec(*presProb->icupp, row)));
         const double implied_upperbound = (rhs - min_act) / coeff;
         upper_implied = PIPSisLE(implied_upperbound, xupp);
      }
      else
      {
         assert( row.inEqSys() || !PIPSisZero(getSimpleVecFromRowStochVec(*presProb->iclow, row)));
         const double implied_upperbound = (lhs - max_act) / coeff;
         upper_implied = PIPSisLE(implied_upperbound, xupp);
      }
   }

   /* calculate an check implied lower bound */
   if( min_ubndd == 0 && xlow != INF_NEG )
   {
      if( coeff < 0.0 )
      {
         assert( row.inEqSys() || !PIPSisZero(getSimpleVecFromRowStochVec(*presProb->iclow, row)));
         const double implied_lowerbound = (lhs - max_act) / coeff;
         lower_implied = PIPSisLE(xlow, implied_lowerbound);
      }
      else
      {
         assert( row.inEqSys() || !PIPSisZero(getSimpleVecFromRowStochVec(*presProb->icupp, row)));
         const double implied_lowerbound = (rhs - min_act) / coeff;
         lower_implied = PIPSisLE(xlow, implied_lowerbound);
      }
   }
}

void PresolveData::fixEmptyColumn(const INDEX& col, double val)
{
   assert( col.isCol() );
   assert( col.hasValidNode(nChildren) );

   if(postsolver)
   {
      const double obj_value = getSimpleVecFromColStochVec(*presProb->g, col);
      double xlow, xupp;
      getColBounds(col, xlow, xupp);

      assert(PIPSisLE(xlow, val) );
      assert(PIPSisLE(val, xupp));

      postsolver->notifyFixedEmptyColumn(col, val, obj_value, xlow, xupp);
      assert( postsolver->wasColumnRemoved(col) );
   }

   removeColumn(col, val);

   assert( getNnzsCol(col) == 0 );
}

void PresolveData::startColumnFixation()
{
   if( postsolver )
      postsolver->putLinkingRowIneqSyncEvent();
}

void PresolveData::fixColumn( const INDEX& col, double value)
{
   assert( col.isCol() );
   assert( col.hasValidNode(nChildren) );

   if( col.isLinkingCol() )
      assert( PIPS_MPIisValueEqual(col.getIndex()) );

   if( track_col && tracked_col == col )
      std::cout << "TRACKING_COLUMN: " << col << " got fixed to " << value << std::endl;
   
   /* current upper and lower bound as well as column - linking variables have to be done by all processes simultaneously because communication in postsolve is required */
   if(postsolver)
   {
      const double obj_coeff = getSimpleVecFromColStochVec(*presProb->g, col);
      postsolver->notifyFixedColumn(col, value, obj_coeff, getSystemMatrix(EQUALITY_SYSTEM), getSystemMatrix(INEQUALITY_SYSTEM));
   }

#ifndef NDEBUG
   double xlow, xupp;
   getColBounds(col, xlow, xupp);

   assert( xlow != INF_NEG );
   assert( xupp != INF_POS );
   assert( PIPSisEQ(xlow, xupp, 1e-10) );
   assert( PIPSisEQ(xlow, value, 1e-10) );
#endif

   removeColumn(col, value);

   if( !col.isLinkingCol() )
      assert( getNnzsCol(col) == 0 );
}

/* an empty col indicates that a singleton row is being removed that is not on this process for Ai and Ci call the synced method */
void PresolveData::removeSingletonRow(const INDEX& row, const INDEX& col, double xlow_new, double xupp_new, double coeff)
{
   assert( row.isRow() );
   assert( col.isCol() || col.isEmpty() );
   assert( row.hasValidNode(nChildren) );
   assert( col.hasValidNode(nChildren) );

   if( row.isLinkingRow() )
      assert( PIPS_MPIisValueEqual(row.getIndex()) );
   if( col.isLinkingCol() )
      assert( row.getNode() == -1 );

   assert( getNnzsRow(row) == 1 );

   if( !col.isEmpty() )
   {
      assert( xlow_new != INF_NEG || xupp_new != INF_POS );
      assert( (row.getSystemType() == EQUALITY_SYSTEM && xlow_new == xupp_new) ||
         (xlow_new == INF_NEG || xupp_new == INF_POS) );

      /* check for infeasibility of the newly found bounds */
      checkBoundsInfeasible(col, xlow_new, xupp_new);

      double xlow_old, xupp_old;
      getColBounds(col, xlow_old, xupp_old);

      /* adjust bounds of column - singletons columns will always be used here since we want to remove the corresponding row */
      bool tightened = updateColBounds(col, xlow_new, xupp_new);

      /* notify postsolver */
      if(tightened && postsolver)
         postsolver->notifySingletonRowBoundsTightened(row, col, xlow_old, xupp_old, xlow_new, xupp_new, coeff);
   }

   /* remove redundant row */
   removeRedundantRow( row );
}

void PresolveData::removeSingletonRowSynced(const INDEX& row, const INDEX& col, double xlow_new, double xupp_new, double coeff)
{
   assert( row.isRow() || row.isEmpty() );
   assert( col.isCol() );
   assert( row.isEmpty() || row.hasValidNode(nChildren) );
   assert( col.hasValidNode(nChildren) );
   assert( col.isLinkingCol() );
   assert( row.isEmpty() || !row.isLinkingRow() );

   assert( PIPS_MPIisValueEqual(col.getIndex()) );
   assert( PIPS_MPIisValueEqual(xlow_new) );
   assert( PIPS_MPIisValueEqual(xupp_new) );

   if( row.isRow() )
   {
      assert( getNnzsRow(row) == 1 );
      assert( (row.inEqSys() && xlow_new == xupp_new) ||
         (xlow_new == INF_NEG || xupp_new == INF_POS) );
   }
   assert( PIPS_MPIgetSum( row.isRow() ? 1 : 0 ) == 1 );

   assert( xlow_new != INF_NEG || xupp_new != INF_POS );

   /* check for infeasibility of the newly found bounds */
   checkBoundsInfeasible(col, xlow_new, xupp_new);

   double xlow_old, xupp_old;
   getColBounds(col, xlow_old, xupp_old);

   /* adjust bounds of column - singletons columns will always be used here since we want to remove the corresponding row */
   bool tightened = updateColBounds(col, xlow_new, xupp_new);

   /* notify postsolver */
   if( tightened && postsolver )
      postsolver->notifySingletonRowBoundsTightened(row, col, xlow_old, xupp_old, xlow_new, xupp_new, coeff);

   /* remove redundant row */
   if( row.isRow() )
      removeRedundantRow( row );
}


void PresolveData::startBoundTightening()
{
   if( postsolver )
      postsolver->beginBoundTightening();

   in_bound_tightening = true;
}

bool PresolveData::rowPropagatedBounds( const INDEX& row, const INDEX& col, double xlow_new, double xupp_new)
{
   assert( in_bound_tightening );
   assert( row.isRow() || row.isEmpty() );
   assert( col.isCol() );
   assert( col.hasValidNode(nChildren) );

   const bool at_root_node = row.isEmpty() ? false : col.isLinkingCol() && row.getNode() == -1;

#ifndef NDEBUG
   if( col.isLinkingCol() && !at_root_node )
   {
      assert( PIPS_MPIisValueEqual(col.getIndex()) );
      assert( PIPS_MPIisValueEqual(xlow_new) );
      assert( PIPS_MPIisValueEqual(xupp_new) );
      const int my_tightening = row.isEmpty() ? 0 : 1;
      assert( PIPS_MPIgetSum(my_tightening) == 1);
   }
   else if( col.isLinkingCol() && at_root_node )
   {
      assert( PIPS_MPIisValueEqual(col.getIndex()) );
      assert( PIPS_MPIisValueEqual(xlow_new) );
      assert( PIPS_MPIisValueEqual(xupp_new) );
   }
#endif

   if( xlow_new == INF_NEG && xupp_new == INF_POS )
      return false;


   if( row.isRow() )
      assert( row.hasValidNode(nChildren) );

   double xlow_old, xupp_old;
   getColBounds(col, xlow_old, xupp_old);

   if( track_col && tracked_col == col )
   {
      std::cout << "TRACKING_COLUMN: new bounds [" << xlow_new << ", " << xupp_new << "] propagated for column " << col
         << " from row " << row << ":" << std::endl;
      std::cout << "\tbounds were [" << xlow_old << ", " << xupp_old << "]" << std::endl;
   }

   checkBoundsInfeasible(col, xlow_new, xupp_new);

   /* adjust bounds */
   bool upper_bound_changed = false;
   bool lower_bound_changed = false;

   const double min_impact_bound_change = feastol;
   // we do not tighten bounds if impact is too low or bound is bigger than threshold_bound_tightening
   /// set upper bound
   if( PIPSisLT(xupp_new, xupp_old) )
   {
      if( PIPSisLT(std::fabs(xupp_new), limit_max_bound_accepted) && ( xupp_old == INF_POS || PIPSisLE(min_impact_bound_change, xupp_old - xupp_new) ) )
      {
         assert( xupp_new != INF_POS );
         assert( PIPSisLE(xupp_new, xupp_old) );
         if( updateColUpperBound(col, xupp_new) )
         {
            /* store node and row that implied the bound (necessary for resetting bounds later on) */
            markRowAsImplyingColumnBound(col, row, true);
            upper_bound_changed = true;
         }
      }
   }
   /// set lower bound
   // if( fabs(ubx) < 1e8 && (PIPSisZero(ixupp) || feastol * 1e3 <= fabs(xupp- ubx) ) )
   if( PIPSisLT( xlow_old, xlow_new ) )
   {
      if( PIPSisLT(std::fabs(xlow_new), limit_max_bound_accepted) && ( xlow_old == INF_NEG || PIPSisLT(min_impact_bound_change, xlow_old - xlow_new)) )
      {
         assert( xlow_new != INF_NEG );
         assert( PIPSisLE(xlow_old, xlow_new) );
         if( updateColLowerBound(col, xlow_new) )
         {
            /* store node and row that implied the bound (necessary for resetting bounds later on) */
            markRowAsImplyingColumnBound(col, row, false);
            lower_bound_changed = true;
         }
      }
   }


   /// if a linking variable was tightened (with a row that is not in B0) it's dual has to be comunicated in postsolve
   /// if a linking row was used to tighten a bound it has to be stored on all processes and used for postsolve by all processes (if the bound was tight)
   /// for linking variables we need to figure out a row that has been used for it's tightening

   /// every process should have the same root node data thus all of them should propagate their rows similarly
   if( (lower_bound_changed || upper_bound_changed) && col.isLinkingCol() )
      assert(outdated_linking_var_bounds == true);

   /// linking rows require a different postsolve event since propagating linking rows need to be stored by every process and need to be postsolved at the same time
   if( postsolver )
   {
      if( lower_bound_changed )
         postsolver->notifyRowPropagatedBound( row, col, xlow_old, xlow_new, false, getSystemMatrix(row.getSystemType()) );
      if( upper_bound_changed )
         postsolver->notifyRowPropagatedBound( row, col, xupp_old, xupp_new, true, getSystemMatrix(row.getSystemType()) );
   }

   if( (lower_bound_changed || upper_bound_changed) && row.isLinkingRow() && !at_root_node )
   {
      postsolve_linking_row_propagation_needed = true;

      if( row.inEqSys() )
         store_linking_row_boundTightening_A[row.getIndex()] = 1;
      else
         store_linking_row_boundTightening_C[row.getIndex()] = 1;
   }

   outdated_linking_var_bounds = false;
   return (lower_bound_changed || upper_bound_changed);
}

void PresolveData::endBoundTightening()
{
   syncPostsolveOfBoundsPropagatedByLinkingRows();

   if( postsolver )
      postsolver->endBoundTightening(store_linking_row_boundTightening_A, store_linking_row_boundTightening_C, getSystemMatrix(EQUALITY_SYSTEM), getSystemMatrix(INEQUALITY_SYSTEM));

   std::fill(store_linking_row_boundTightening_A.begin(), store_linking_row_boundTightening_A.end(), 0);
   std::fill(store_linking_row_boundTightening_C.begin(), store_linking_row_boundTightening_C.end(), 0);

   in_bound_tightening = false;
}

void PresolveData::checkBoundsInfeasible(const INDEX& col, double xlow_new, double xupp_new) const
{
   assert( xlow_new != INF_NEG || xupp_new != INF_POS );
   assert( col.isCol() );

   /* get current bounds */
   double xlow, xupp;
   getColBounds(col, xlow, xupp);

   if( !PIPSisLE( std::max(xlow_new, xlow), std::min(xupp_new, xupp)) )
   {
      std::cout << "[" << xlow_new << ", " << xupp_new << "] not in [" << xlow << ", " << xupp << "]" << std::endl;
      PIPS_MPIabortInfeasible("Detected infeasible new bounds!", "PresolveData.C", "checkBoundsInfeasible");
   }
}

void PresolveData::syncPostsolveOfBoundsPropagatedByLinkingRows()
{
   PIPS_MPIgetLogicOrInPlace(postsolve_linking_row_propagation_needed);

   if( !postsolve_linking_row_propagation_needed )
      return;

   // sync array of rows that propagated
   PIPS_MPIsumArrayInPlace(store_linking_row_boundTightening_A);
   PIPS_MPIsumArrayInPlace(store_linking_row_boundTightening_C);

   postsolve_linking_row_propagation_needed = false;
}

/** tightening row1 bounds with row2 where row1 = factor * row2 */
void PresolveData::tightenRowBoundsParallelRow(const INDEX& row1, const INDEX& row2, double clow_new, double cupp_new, double factor)
{
   assert( row1.isRow() && row2.isRow() );
   assert( row1.hasValidNode(nChildren) && row2.hasValidNode(nChildren) );
   assert( !row1.isLinkingRow() && !row2.isLinkingRow() );
   assert( row1.inInEqSys() && row2.inInEqSys() );
   assert( clow_new != INF_NEG || cupp_new != INF_POS );
   assert( PIPSisLE(clow_new, cupp_new) );

   if( track_row && tracked_row == row1 )
   {
      std::cout << "TRACKING_ROW: before RHS LHS adjustment" << std::endl;
      writeRowLocalToStreamDense(std::cout, row1 );
   }

   double clow, cupp;
   getRowBounds(row1, clow, cupp);

   assert( PIPSisLE( std::max(clow, clow_new), std::min(cupp, cupp_new) ) );
   if( clow_new != INF_NEG )
      assert( PIPSisLT( clow, clow_new ) );
   if( cupp_new != INF_POS )
      assert( PIPSisLT( cupp_new, cupp ) );

   if( postsolver )
      postsolver->notifyParallelRowsBoundsTightened(row1, row2, clow, cupp, clow_new, cupp_new, factor);

   if( clow_new != INF_NEG )
   {
      getSimpleVecFromRowStochVec(*presProb->iclow, row1) = 1.0;
      getSimpleVecFromRowStochVec(*presProb->bl, row1) = clow_new;
   }

   if( cupp_new != INF_POS )
   {
      getSimpleVecFromRowStochVec(*presProb->icupp, row1) = 1.0;
      getSimpleVecFromRowStochVec(*presProb->bu, row1) = cupp_new;
   }

   if( track_row && tracked_row == row1 )
   {
      std::cout << "TRACKING_ROW: after RHS LHS adjustment " << std::endl;
      writeRowLocalToStreamDense(std::cout, row1);
   }
}

/** this methods does not call any postsolve procedures but simply changes the bounds (lhs, rhs) of either A or B by value */
void PresolveData::adjustMatrixRhsLhsBy(const INDEX& row, double value, bool at_root )
{
   assert( row.isRow() );
   assert( row.hasValidNode(nChildren) );
   if( track_row && tracked_row == row )
   {
      std::cout << "TRACKING_ROW: RHS LHS of row " << row << " being adjusted by " << value << std::endl;
      writeRowLocalToStreamDense(std::cout, row );
   }

   if( PIPSisEQ(value, 0.0) )
      return;

   if( row.isLinkingRow() && !at_root )
   {
      outdated_lhsrhs = true;
      row.inEqSys() ? (*bound_chgs_A)[row.getIndex()] += value : (*bound_chgs_C)[row.getIndex()] += value;
      return;
   }

   if( row.inEqSys() )
   {
      getSimpleVecFromRowStochVec(*presProb->bA, row) += value;
   }
   else
   {
      if( PIPSisEQ(getSimpleVecFromRowStochVec(*presProb->icupp, row), 1.0) )
         getSimpleVecFromRowStochVec(*presProb->bu, row) += value;

      if( PIPSisEQ(getSimpleVecFromRowStochVec(*presProb->iclow, row), 1.0 ) )
         getSimpleVecFromRowStochVec(*presProb->bl, row) += value;
   }

   if( track_row && tracked_row == row )
   {
      std::cout << "TRACKING_ROW: after RHS LHS adjustment " << std::endl;
      writeRowLocalToStreamDense(std::cout, row);
   }
}

void PresolveData::reduceNnzCounterRowBy(const INDEX& row, int amount, bool at_root)
{
   assert( row.isRow() );
   assert( row.hasValidNode(nChildren) );
   if( at_root )
      assert( row.getNode() == -1 );
   assert( 0 <= amount );
   changeNnzCounterRow(row, -amount, at_root);
}

void PresolveData::increaseNnzCounterRowBy(const INDEX& row, int amount, bool at_root)
{
   assert( row.isRow() );
   assert( row.hasValidNode(nChildren) );
   if( at_root )
      assert( row.getNode() == -1 );
   assert( 0 <= amount );
   changeNnzCounterRow(row, amount, at_root);
}

void PresolveData::changeNnzCounterRow(const INDEX& row, int amount, bool at_root)
{
   assert( row.isRow() );
   assert( row.hasValidNode(nChildren) );

   if( amount == 0 )
      return;

   /* linking constraints get stored */
   if( row.isLinkingRow() && !at_root )
   {
      int& chgs = (row.inEqSys() ? *nnzs_row_A_chgs : *nnzs_row_C_chgs)[row.getIndex()];

      chgs += amount;
      outdated_nnzs = true;

      assert( 0 <= getNnzsRow(row) );
   }
   else
   {
      int& nnzs_row = getSimpleVecFromRowStochVec(row.inEqSys() ? *nnzs_row_A : *nnzs_row_C, row);

      nnzs_row += amount;
      assert( 0 <= nnzs_row );

      if( nnzs_row == 1 )
         singleton_rows.push( row );
   }
}

void PresolveData::reduceNnzCounterColumnBy(const INDEX& col, int amount, bool at_root)
{
   assert( col.isCol() );
   assert( col.hasValidNode(nChildren) );
   assert( 0 <= amount );
   changeNnzCounterColumn(col, -amount, at_root);
}

void PresolveData::increaseNnzCounterColumnBy(const INDEX& col, int amount, bool at_root)
{
   assert( col.isCol() );
   assert( col.hasValidNode(nChildren) );
   assert( 0 <= amount );
   changeNnzCounterColumn(col, amount, at_root);
}

void PresolveData::changeNnzCounterColumn(const INDEX& col, int amount, bool at_root)
{
   assert( col.isCol() );
   assert( col.hasValidNode(nChildren) );

   if( amount == 0 )
      return;

   /* linking constraints get stored */
   if( col.isLinkingCol() && !at_root )
   {
      int& chgs = (*nnzs_col_chgs)[col.getIndex()];
      chgs += amount;
      assert( 0 <= getNnzsCol(col) );
      outdated_nnzs = true;
   }
   else
   {
      int& nnzs = getSimpleVecFromColStochVec( *nnzs_col, col );
      nnzs += amount;
      assert(0 <= nnzs );

      if( nnzs == 1)
         singleton_cols.push( col );
   }
}

/** removes a column from the whole system A, C by fixing x to fixation
 * updates non-zero counters, rhs, lhs, objective offset and activities
 * does not call postsolve routines
 */
void PresolveData::removeColumn(const INDEX& col, double fixation)
{
   assert( col.isCol() );
   /* have to tell the postsolver first if existent */
   assert( wasColumnRemoved( col ) );
   if( col.isLinkingCol() )
      assert( PIPS_MPIisValueEqual(col.getIndex()) );
   assert( col.hasValidNode(nChildren) );
   const int dummy_index = std::numeric_limits<int>::infinity();

   if( col.isLinkingCol() )
   {
      /* B0 */
      removeColumnFromMatrix( INDEX(ROW, -1, dummy_index, false, EQUALITY_SYSTEM), col, fixation);

      /* D0 */
      removeColumnFromMatrix( INDEX(ROW, -1, dummy_index, false, INEQUALITY_SYSTEM), col, fixation);

      /* Bl0 */
      if( hasLinking(EQUALITY_SYSTEM) )
         removeColumnFromMatrix( INDEX(ROW, -1, dummy_index, true, EQUALITY_SYSTEM), col, fixation);

      /* Dl0 */
      if( hasLinking(INEQUALITY_SYSTEM) )
         removeColumnFromMatrix( INDEX(ROW, -1, dummy_index, true, INEQUALITY_SYSTEM), col, fixation);

      for(int i = 0; i < nChildren; ++i)
      {
         /* Ai, Ci */
         if( !nodeIsDummy(i) )
         {
            removeColumnFromMatrix( INDEX(ROW, i, dummy_index, false, EQUALITY_SYSTEM), col, fixation);
            removeColumnFromMatrix( INDEX(ROW, i, dummy_index, false, INEQUALITY_SYSTEM), col, fixation);
         }
      }
   }
   else
   {
      removeColumnFromMatrix( INDEX(ROW, col.getNode(), dummy_index, false, EQUALITY_SYSTEM), col, fixation);
      removeColumnFromMatrix( INDEX(ROW, col.getNode(), dummy_index, false, INEQUALITY_SYSTEM), col, fixation);

      if( hasLinking(EQUALITY_SYSTEM) )
         removeColumnFromMatrix( INDEX(ROW, -1, dummy_index, true, EQUALITY_SYSTEM), col, fixation);
      if( hasLinking(INEQUALITY_SYSTEM) )
         removeColumnFromMatrix( INDEX(ROW, -1, dummy_index, true, INEQUALITY_SYSTEM), col, fixation);
   }

   /* adjust objective function */
   if( !col.isLinkingCol() || my_rank == 0 )
   {
      double objective_factor = getSimpleVecFromColStochVec(*presProb->g, col);
      obj_offset_chgs += objective_factor * fixation;
   }

   /* mark column as removed */
   markColumnRemoved( col );
}

/** remove column - adjust lhs, rhs and activity as well as nnz_counters */
void PresolveData::removeColumnFromMatrix(const INDEX& dummy_row, const INDEX& col, double fixation)
{
   assert( dummy_row.isRow() );
   assert( dummy_row.getIndex() == std::numeric_limits<int>::infinity() );
   assert( dummy_row.hasValidNode(nChildren) );
   assert( col.isCol() );
   assert( col.hasValidNode(nChildren) );
   const bool at_root = (dummy_row.getNode() == -1 && col.isLinkingCol());

   if( at_root )
      assert( PIPS_MPIisValueEqual( col.getIndex() ) );

   SparseGenMatrix* mat = getSparseGenMatrix(dummy_row, col);
   const SparseStorageDynamic& matrix_transp = mat->getStorageDynamicTransposedRef();

   assert( col.getIndex() < matrix_transp.getM() );

   /* remove all entries in column from the sparse storage dynamic */
   for( int j = matrix_transp.getRowPtr(col.getIndex()).start; j < matrix_transp.getRowPtr(col.getIndex()).end; j++ )
   {
      const int row_index = matrix_transp.getJcolM(j);
      const double coeff = matrix_transp.getMat(j);

      const INDEX row = INDEX(ROW, dummy_row.getNode(), row_index, dummy_row.isLinkingRow(), dummy_row.getSystemType());
      assert( !PIPSisEQ(0.0, coeff) );

      if( track_row && tracked_row == row )
         std::cout << "TRACKING_ROW: fixation of " << col << " in " << row << std::endl;

      /* remove the entry, adjust activity and row counters and rhs/lhs */
      if( postsolver )
      {
         postsolver->notifyRowModified( row );
         postsolver->notifyColModified( col );
      }

      reduceNnzCounterRowBy( row, 1, at_root );

      adjustMatrixRhsLhsBy( row, - coeff * fixation, at_root );

      adjustRowActivityFromDeletion( row, col, coeff );

      if( track_row && tracked_row == row )
      {
         std::cout << "TRACKING_ROW: after removal of column" << std::endl;
         writeRowLocalToStreamDense(std::cout, row);

         double act_min, act_max;
         int act_min_ubndd, act_max_ubndd;
         getRowActivities(row, act_max, act_min, act_max_ubndd, act_min_ubndd);

         std::cout << "TRACKING_ROW: New activity of row " << row << std::endl;
         std::cout << "\tnew min/max activity is: " << act_min << "/" << act_max << ", min/max unbounded counters are " << act_min_ubndd << "/" << act_max_ubndd << std::endl;
      }
   }

   /* adjust column counters */
   reduceNnzCounterColumnBy( col, matrix_transp.getRowPtr(col.getIndex()).end - matrix_transp.getRowPtr(col.getIndex()).start, at_root);

   /* delete col in matrix and transposed */
   mat->removeCol( col.getIndex() );
}

/* tighten bounds of singleton column col1 in row implied via singleton column col2 in row2 */
void PresolveData::tightenBoundsNearlyParallelRows( const INDEX& row1, const INDEX& row2, const INDEX& col1, const INDEX& col2, double xlow_new, double xupp_new, double scalar,
      double translation, double parallel_factor )
{
   assert( row1.isRow() );
   assert( row2.isRow() );
   assert( !row1.isLinkingRow() );
   assert( !row2.isLinkingRow() );
   assert( row1.inEqSys() );

   if( row2.inInEqSys() )
   {
      assert( col1.isCol() );
      assert( col2.isEmpty() );
      assert(scalar == INF_POS);
      assert(translation == INF_POS);
   }
   else
   {
      assert( col1.isCol() || col1.isEmpty() );
      assert( col2.isCol() );
   }

   assert( !wasRowRemoved(row1) );
   assert( !wasRowRemoved(row2) );

   if( col1.isCol() )
      assert( !wasColumnRemoved(col1) );
   if( col2.isCol() )
      assert( !wasColumnRemoved(col2) );

   if( xlow_new == INF_NEG && xupp_new == INF_POS )
      return;

   if( postsolver )
   {
      const double xlow_col1 = col1.isCol() ? getSimpleVecFromColStochVec(*presProb->blx, col1) : INF_NEG;
      const double xupp_col1 = col1.isCol() ? getSimpleVecFromColStochVec(*presProb->bux, col1) : INF_POS;

      const double xlow_col2 = col2.isCol() ? getSimpleVecFromColStochVec(*presProb->blx, col2) : INF_NEG;
      const double xupp_col2 = col2.isCol() ? getSimpleVecFromColStochVec(*presProb->bux, col2) : INF_POS;

      /* one of the bounds has to be a tightening one */
      assert( PIPSisLE(xlow_col1, xlow_new) || PIPSisLE(xupp_new, xupp_col1) );

      if( col1.isCol() )
      {
         if( PIPSisZero(getSimpleVecFromColStochVec(*presProb->ixlow, col1)) )
            assert(xlow_col1 == INF_NEG);
         if( PIPSisZero(getSimpleVecFromColStochVec(*presProb->ixupp, col1)) )
            assert(xupp_col1 == INF_POS);
      }
      if( col2.isCol() )
      {
         if( PIPSisZero(getSimpleVecFromColStochVec(*presProb->ixlow, col2)) )
            assert(xlow_col2 == INF_NEG);
         if( PIPSisZero(getSimpleVecFromColStochVec(*presProb->ixupp, col2)) )
            assert(xupp_col1 == INF_POS);
      }

      const double coeff_col1 = col1.isCol() ? getRowCoeff(row1, col1) : 0.0;
      const double coeff_col2 = col2.isCol() ? getRowCoeff(row2, col2) : 0.0;

      if( col2.isCol() )
         assert( !PIPSisZero(coeff_col2) );

      const double rhs = row2.inInEqSys() ? getSimpleVecFromRowStochVec(*presProb->bA, row1) : 0.0;

      const double clow = row2.inInEqSys() ? ( PIPSisZero(getSimpleVecFromRowStochVec(*presProb->iclow, row2)) ? INF_NEG : getSimpleVecFromRowStochVec(*presProb->bl, row2) )
            : 0.0;
      const double cupp = row2.inInEqSys() ? ( PIPSisZero(getSimpleVecFromRowStochVec(*presProb->icupp, row2)) ? INF_POS : getSimpleVecFromRowStochVec(*presProb->bu, row2) )
            : 0.0;

      postsolver->notifyNearlyParallelRowBoundsTightened(row1, row2, col1, col2, xlow_col1, xupp_col1, xlow_col2, xupp_col2, coeff_col1, coeff_col2, scalar,
            translation, parallel_factor, rhs, clow, cupp);
   }

   /* adjust bounds and remove inequality row */
   updateColBounds(col1, xlow_new, xupp_new);
}

void PresolveData::startParallelRowPresolve()
{
   if( postsolver )
      postsolver->putLinkingVarsSyncEvent();
}

/* a singleton variable is substituted out of the problem and then it's original row can be removed from the problem */
void PresolveData::substituteVariableNearlyParallelRows( const INDEX& row1, const INDEX& row2, const INDEX& col1, const INDEX& col2,
   double scalar, double translation, double parallelity )
{
   assert( row1.isRow() && row2.isRow() );
   assert( col1.isCol() || (col1.isEmpty() && row1.inEqSys()) );
   assert( col2.isCol() );
   assert( row1.getNode() == row2.getNode() );

   if( row1.inInEqSys() )
      assert( PIPSisZero(translation) );

   // TODO: track row

   const double obj_col1 = col1.isCol() ? getSimpleVecFromColStochVec(*presProb->g, col1) : INF_POS;
   const double obj_col2 = getSimpleVecFromColStochVec(*presProb->g, col2);

   if( postsolver )
   {
      const double xlow_col2 = getSimpleVecFromColStochVec(*presProb->blx, col2);
      const double xupp_col2 = getSimpleVecFromColStochVec(*presProb->bux, col2);

      const double coeff_col1 = col1.isCol() ? getRowCoeff(row1, col1) : 0.0;
      const double coeff_col2 = getRowCoeff(row2, col2);

      if( PIPSisZero(getSimpleVecFromColStochVec(*presProb->ixlow, col2)) )
         assert(xlow_col2 == INF_NEG);
      if( PIPSisZero(getSimpleVecFromColStochVec(*presProb->ixupp, col2)) )
         assert(xupp_col2 == INF_POS);

      postsolver->notifyNearlyParallelRowSubstitution(row1, row2, col1, col2, scalar, translation, obj_col1, obj_col2, xlow_col2, xupp_col2, coeff_col1, coeff_col2, parallelity );
   }

   /* adapt the objective vector and the row2 */
   const double val_offset = translation * obj_col2;
   const double change_obj_var1 = scalar * obj_col2;

   bool col1_at_root = col1.isEmpty() ? false : (col1.isLinkingCol() && row1.getNode() == -1);
   bool col2_at_root = (col2.isLinkingCol() && row1.getNode() == -1);

   assert( col2.isCol() );
   const double coeff_col2 = getRowCoeff(row2, col2);

   /* add col1 with coefficient scalar * coeff_col2 to row */
   if( col1.isCol() )
   {
      const double coeff_col1_row2 = scalar * coeff_col2;
      addCoeffColToRow( coeff_col1_row2, col1, row2 );
   }

   /* remove old entry from matrix */
   const double offset = translation * coeff_col2;

   getSparseGenMatrix(row2, col2)->removeEntryAtRowCol(row2.getIndex(), col2.getIndex());
   reduceNnzCounterRowBy( row2, 1, false );
   reduceNnzCounterColumnBy( col2, 1, col2_at_root );

   /* adjust right hand side / left hand side by  -d * coeff_col2 */
   if( row2.inInEqSys() )
   {
      if( !PIPSisZero(getSimpleVecFromRowStochVec(*presProb->iclow, row2)) )
         getSimpleVecFromRowStochVec(*presProb->blx, row2) -= offset;
      if( !PIPSisZero(getSimpleVecFromRowStochVec(*presProb->icupp, row2)) )
         getSimpleVecFromRowStochVec(*presProb->bux, row2) -= offset;
   }
   else
      getSimpleVecFromRowStochVec(*presProb->bA, row2) -= offset;

   if( col1.isCol() )
   {
      if( !col1.isLinkingCol() || col1_at_root )
      {
         getSimpleVecFromColStochVec(*presProb->g, col1) += change_obj_var1;
         objOffset += val_offset;
      }
      else
      {
         /* var1 is a linking variable - store objective adaption and allreduce them */
         (*objective_vec_chgs)[col1.getIndex()] += change_obj_var1;
         outdated_obj_vector = true;
         obj_offset_chgs += val_offset;
      }
   }

   if( !col2.isLinkingCol() || col2_at_root )
      getSimpleVecFromColStochVec(*presProb->g, col2) = 0.0;
   else
   {
      (*objective_vec_chgs)[col2.getIndex()] -= getSimpleVecFromColStochVec(*presProb->g, col2);
      assert( PIPSisZero( (*objective_vec_chgs)[col2.getIndex()] + getSimpleVecFromColStochVec(*presProb->g, col2) ) );
      outdated_obj_vector = true;
   }
}

void PresolveData::removeRedundantParallelRow( const INDEX& rm_row, const INDEX& par_row )
{
   assert( rm_row.isRow() );
   assert( par_row.isRow() );
   assert( !wasRowRemoved(rm_row) );
   assert( !wasRowRemoved(par_row) );

   if(postsolver)
   {
      const double cupp_row1 = rm_row.inEqSys() ? getSimpleVecFromRowStochVec(*presProb->bA, rm_row) : getSimpleVecFromRowStochVec(*presProb->bu, rm_row);
      const double clow_row1 = rm_row.inEqSys() ? cupp_row1 : getSimpleVecFromRowStochVec(*presProb->bl, rm_row);
      const int iclow_row1 = rm_row.inEqSys() ? 1 : getSimpleVecFromRowStochVec(*presProb->iclow, rm_row);
      const int icupp_row1 = rm_row.inEqSys() ? 1 : getSimpleVecFromRowStochVec(*presProb->icupp, rm_row);

#ifndef NDEBUG
      // TODO : assert that rows are actually parallel...
#endif
      postsolver->notifyRedundantRow(rm_row, iclow_row1, icupp_row1, clow_row1, cupp_row1, getSystemMatrix( rm_row.getSystemType() ));

      assert( wasRowRemoved(rm_row) );
   }

   if( track_row && tracked_row == rm_row )
      std::cout << "TRACKING_ROW: removal of tracked row as parallel and redundant row to " << par_row << std::endl;

   removeRow( rm_row );
}

void PresolveData::removeRedundantRow( const INDEX& row )
{
   assert( row.isRow() );

   if(postsolver)
   {
      assert( !wasRowRemoved( row ) );

      const double rhs = row.inEqSys() ? getSimpleVecFromRowStochVec(*presProb->bA, row) : getSimpleVecFromRowStochVec(*presProb->bu, row);
      const double lhs = row.inEqSys() ? rhs : getSimpleVecFromRowStochVec(*presProb->bl, row);
      const int iclow = row.inEqSys() ? 1 : getSimpleVecFromRowStochVec(*presProb->iclow, row);
      const int icupp = row.inEqSys() ? 1 : getSimpleVecFromRowStochVec(*presProb->icupp, row);

#ifndef NDEBUG
      double max_act = 0;
      double min_act = 0;

      int max_ubndd = 0;
      int min_ubndd = 0;

      getRowActivities(row, max_act, min_act, max_ubndd, min_ubndd);

      if(iclow)
      {
         /// a singleton row with linking var entry or a redundant row
         assert(min_ubndd <= 1);

         if( min_ubndd == 0 )
            assert(PIPSisLEFeas(lhs, min_act));
      }
      if(icupp)
      {
         assert(max_ubndd <= 1);

         if( max_ubndd == 0 )
            assert(PIPSisLEFeas(max_act, rhs));
      }
#endif

      assert( PIPSisLE(0.0, iclow) );
      assert( PIPSisLE(0.0, icupp) );
      assert( PIPSisLT(0.0, iclow + icupp) );

      postsolver->notifyRedundantRow(row, iclow, icupp, lhs, rhs, getSystemMatrix( row.getSystemType() ));

      assert(postsolver->wasRowRemoved(row));
   }
 
   if( track_row && tracked_row == row )
      std::cout << "TRACKING_ROW: removal of tracked row as redundant row" << std::endl;

   removeRow( row );
}

void PresolveData::startSingletonColumnPresolve()
{
   if( postsolver )
   {
      postsolver->putLinkingVarsSyncEvent();
      postsolver->putLinkingRowIneqSyncEvent();
   }
}

/** dual fixing for a singleton column */
void PresolveData::fixColumnInequalitySingleton( const INDEX& col, const INDEX& row, double value, double coeff )
{
   assert( col.isCol() );

   const double ixlow = getSimpleVecFromColStochVec(*presProb->ixlow, col);
   const double ixupp = getSimpleVecFromColStochVec(*presProb->ixupp, col);
   const double xlow = getSimpleVecFromColStochVec(*presProb->blx, col);
   const double xupp = getSimpleVecFromColStochVec(*presProb->bux, col);

   if( PIPSisZero(ixlow) )
      assert(xlow == INF_NEG);
   if( PIPSisZero(ixupp) )
      assert(xupp == INF_POS);

   assert( PIPSisLE(xlow, value) );
   assert( PIPSisLE(value, xupp) );

   assert( INF_NEG < value && value < INF_POS );

   if( postsolver )
      postsolver->notifyFixedSingletonFromInequalityColumn(col, row, value, coeff, xlow, xupp);

   updateColBounds(col, value, value);
}

void PresolveData::removeFreeColumnSingletonInequalityRow( const INDEX& row, const INDEX& col, double coeff )
{
   assert( row.isRow() );
   assert( col.isCol() );
   assert( row.inInEqSys() );
   assert( !wasColumnRemoved(col) );
   assert( !wasRowRemoved(row) );

   if( row.isLinkingRow() )
      assert( PIPS_MPIisValueEqual(row.getIndex(), MPI_COMM_WORLD) );

   double clow, cupp;
   getRowBounds(row, clow, cupp);
   assert( clow != INF_NEG || cupp != INF_POS );

   if( clow != INF_NEG )
      assert( cupp == INF_POS );
   if( cupp == INF_POS )
      assert( clow == INF_NEG );

   const double rhs = (clow != INF_NEG ) ? -clow : cupp;

   if( track_col && tracked_col == col )
     std::cout << "TRACKING_COLUMN: tracked column removed as free column singleton" << std::endl;
   if( track_row && tracked_row == row )
      std::cout << "TRACKING_ROW: removal of " << row << " since it contained a free column singleton" << std::endl;

   double& ixlow = getSimpleVecFromColStochVec(*(presProb->ixlow), col);
   double& xlow = getSimpleVecFromColStochVec(*(presProb->blx), col);
   double& ixupp = getSimpleVecFromColStochVec(*(presProb->ixupp), col);
   double& xupp = getSimpleVecFromColStochVec(*(presProb->bux), col);

   if( PIPSisZero(ixlow) )
      assert( xlow == INF_NEG );
   if( PIPSisZero(ixupp) )
      assert( xupp == INF_POS );

   if( postsolver )
      postsolver->notifyFreeColumnSingletonInequalityRow( row, col, rhs, coeff, xlow, xupp, getSystemMatrix(row.getSystemType()) );

   removeRow(row);

   /* remove col and mark it for the fix empty columns presolver */
   ixlow = ixupp = 0.0;
   xlow = INF_NEG;
   xupp = INF_POS;
   if( col.getNode() == -1 )
      outdated_linking_var_bounds = true;

   assert( getNnzsCol(col) == 0 );
   if( col.isLinkingCol() && my_rank == 0 )
      assert( PIPSisZero(getSimpleVecFromColStochVec(*presProb->g, col) + (*objective_vec_chgs)[col.getIndex()]) );
   else if( !col.isLinkingCol() )
      assert( PIPSisZero(getSimpleVecFromColStochVec( *presProb->g, col)) );
}

void PresolveData::removeFreeColumnSingletonInequalityRowSynced( const INDEX& row, const INDEX& col, double coeff )
{
   assert( row.isRow() );
#ifndef NDEBUG
   const int my_col = col.isCol();
   assert( 1 == PIPS_MPIgetSum( my_col, MPI_COMM_WORLD ) );
#endif
   MPI_Barrier(MPI_COMM_WORLD);
   assert( PIPS_MPIisValueEqual(row.getIndex(), MPI_COMM_WORLD) );

   if( col.isCol() )
      assert( !wasColumnRemoved(col) );
   assert( !wasRowRemoved(row) );
   assert( row.isLinkingRow() );
   if( col.isEmpty() )
      assert( coeff == 0.0 );

   const double clow = getSimpleVecFromRowStochVec(*presProb->bl, row);
   const double cupp = getSimpleVecFromRowStochVec(*presProb->bu, row);

   assert( clow != INF_NEG || cupp != INF_POS );
   if( clow != INF_NEG )
      assert( cupp == INF_POS);
   if( cupp != INF_POS)
      assert( clow == INF_NEG );

   const double rhs = (clow != INF_NEG ) ? -clow : cupp;

   if( col.isCol() )
   {
      if( track_col && tracked_col == col )
         std::cout << "TRACKING_COLUMN: tracked column removed as free column singleton" << std::endl;
      if( track_row && tracked_row == row )
         std::cout << "TRACKING_ROW: removal of " << row << " since it contained a free column singleton" << std::endl;
   }

   PIPS_MPIgetSumInPlace(coeff, MPI_COMM_WORLD);

   if( col.isCol() )
   {
      double& ixlow = getSimpleVecFromColStochVec(*(presProb->ixlow), col);
      double& xlow = getSimpleVecFromColStochVec(*(presProb->blx), col);
      double& ixupp = getSimpleVecFromColStochVec(*(presProb->ixupp), col);
      double& xupp = getSimpleVecFromColStochVec(*(presProb->bux), col);

      if( PIPSisZero(ixlow) )
         assert( xlow == INF_NEG );
      if( PIPSisZero(ixupp) )
         assert( xupp == INF_POS );
      if( postsolver )
         postsolver->notifyFreeColumnSingletonInequalityRow( row, col, rhs, coeff, xlow, xupp, getSystemMatrix(row.getSystemType()) );

      removeRow(row);

      /* remove col and mark it for the fix empty columns presolver */
      ixlow = xlow = xupp = ixupp = 0;

      assert( PIPSisZero(getSimpleVecFromColStochVec( *presProb->g, col)) );
      assert( getSimpleVecFromColStochVec(*nnzs_col, col) == 0 );

      if( col.isLinkingCol() )
         outdated_linking_var_bounds = true;
   }
   else
   {
      if( postsolver )
         postsolver->notifyFreeColumnSingletonInequalityRow( row, col, rhs, coeff, INF_NEG, INF_POS, getSystemMatrix(row.getSystemType()) );

      removeRow(row);
   }
}

/* for linking rows with non-linking singleton columns this presolver must be called simultaneously - if col is on another process col should an empty index */
void PresolveData::removeImpliedFreeColumnSingletonEqualityRowSynced( const INDEX& row, const INDEX& col )
{
   assert( row.isRow() );
   assert( row.getLinking() );
   assert( col.isCol() || col.isEmpty() );
   assert( !wasRowRemoved(row) );

   if( col.isCol() )
   {
      assert( !col.isLinkingCol() );
      assert( !wasColumnRemoved(col) );
   }

#ifndef NDEBUG
   /* only one process has a local linking column */
   const int my_col = col.isCol();
   assert( PIPS_MPIgetSum(my_col) == 1 );
#endif
   assert( PIPS_MPIisValueEqual(row.getIndex()) );

   if( col.isCol() )
   {
      assert( getSimpleVecFromColStochVec(*nnzs_col, col) == 1);

      if( track_col && tracked_col == col )
        std::cout << "TRACKING_COLUMN: tracked column removed as (implied) free column singleton" << std::endl;
   }

   if( track_row && tracked_row == row )
   {
      std::cout << "TRACKING_ROW: removal of " << row << " since it contained an (implied) free column singleton" << std::endl;
      writeRowLocalToStreamDense(std::cout, row);
   }

   if( row.inInEqSys() )
   {
      assert( PIPSisEQ(getSimpleVecFromRowStochVec( *presProb->bl, row), getSimpleVecFromRowStochVec( *presProb->bu, row)) );
      assert( !PIPSisZero(getSimpleVecFromRowStochVec( *presProb->iclow, row)) );
      assert( !PIPSisZero(getSimpleVecFromRowStochVec( *presProb->icupp, row)) );
   }

   const double rhs = row.inEqSys() ? getSimpleVecFromRowStochVec( *presProb->bA, row) :
      getSimpleVecFromRowStochVec( *presProb->bl, row);

   double obj_coeff = 0.0;
   double col_coeff = 0.0;

   if( col.isCol() )
   {
      assert( !nodeIsDummy(col.getNode()) );
      obj_coeff = col.isLinkingCol() ? getSimpleVecFromColStochVec( *presProb->g, col) + (*objective_vec_chgs)[col.getIndex()] :
            getSimpleVecFromColStochVec(*presProb->g, col);
      col_coeff = getRowCoeff(row, col);
   }

   // TODO: make more efficient
   PIPS_MPIgetSumInPlace(col_coeff);
   PIPS_MPIgetSumInPlace(obj_coeff);

   if( col.isCol() )
   {
      const double& ixlow = getSimpleVecFromColStochVec(*(presProb->ixlow), col);
      const double& xlow = getSimpleVecFromColStochVec(*(presProb->blx), col);
      const double& ixupp = getSimpleVecFromColStochVec(*(presProb->ixupp), col);
      const double& xupp = getSimpleVecFromColStochVec(*(presProb->bux), col);

      if( PIPSisZero(ixlow) )
         assert( xlow == INF_NEG );
      if( PIPSisZero(ixupp) )
         assert( xupp == INF_POS );

      postsolver->notifyFreeColumnSingletonEquality( row, col, rhs, obj_coeff, col_coeff, xlow, xupp, getSystemMatrix(row.getSystemType()) );
   }
   else
      postsolver->notifyFreeColumnSingletonEquality( row, col, rhs, obj_coeff, col_coeff, INF_NEG, INF_POS, getSystemMatrix(row.getSystemType()) );

   /* adapt objective from substitution */
   adaptObjectiveSubstitutedRow( row, col, obj_coeff, col_coeff);

   /* remove row and mark column as empty - will be removed in model cleanup on all processes */
   removeRow( row );

   if( col.isCol() )
   {
      double& ixlow = getSimpleVecFromColStochVec(*(presProb->ixlow), col);
      double& xlow = getSimpleVecFromColStochVec(*(presProb->blx), col);
      double& ixupp = getSimpleVecFromColStochVec(*(presProb->ixupp), col);
      double& xupp = getSimpleVecFromColStochVec(*(presProb->bux), col);

      /* remove col and mark it for the fix empty columns presolver */
      ixlow = xlow = xupp = ixupp = 0;

      assert( PIPSisZero(getSimpleVecFromColStochVec( *presProb->g, col)) );
      assert( getSimpleVecFromColStochVec(*nnzs_col, col) == 0 );
   }
}

void PresolveData::removeImpliedFreeColumnSingletonEqualityRow( const INDEX& row, const INDEX& col )
{
   /* removing multiple linking variables in one run is possible since they get communicated only when the objective changes get communicated too
    * thus no process will remove a linking variable as singleton column with outdated objective vector information
    */
   assert( row.isRow() );
   assert( col.isCol() );
   assert( !nodeIsDummy(row.getNode()) );
   assert( !wasRowRemoved(row) );
   assert( !wasColumnRemoved(col) );

   if( track_col && tracked_col == col )
      std::cout << "TRACKING_COLUMN: tracked column removed as (implied) free column singleton" << std::endl;

   if( track_row && tracked_row == row )
   {
      std::cout << "TRACKING_ROW: removal of " << row << " since it contained an (implied) free column singleton" << std::endl;
      writeRowLocalToStreamDense(std::cout, row);
   }

   if( row.inInEqSys() )
   {
      assert( PIPSisEQ(getSimpleVecFromRowStochVec( *presProb->bl, row), getSimpleVecFromRowStochVec( *presProb->bu, row)) );
      assert( !PIPSisZero(getSimpleVecFromRowStochVec( *presProb->iclow, row)) );
      assert( !PIPSisZero(getSimpleVecFromRowStochVec( *presProb->icupp, row)) );
   }

   const double rhs = row.inEqSys() ? getSimpleVecFromRowStochVec( *presProb->bA, row) : getSimpleVecFromRowStochVec( *presProb->bl, row);
   const double obj_coeff = col.isLinkingCol() ? getSimpleVecFromColStochVec( *presProb->g, col) + (*objective_vec_chgs)[col.getIndex()] :
      getSimpleVecFromColStochVec(*presProb->g, col);
   const double col_coeff = getRowCoeff(row, col);

   double& ixlow = getSimpleVecFromColStochVec(*(presProb->ixlow), col);
   double& xlow = getSimpleVecFromColStochVec(*(presProb->blx), col);
   double& ixupp = getSimpleVecFromColStochVec(*(presProb->ixupp), col);
   double& xupp = getSimpleVecFromColStochVec(*(presProb->bux), col);

   if( PIPSisZero(ixlow) )
      assert( xlow == INF_NEG );
   if( PIPSisZero(ixupp) )
      assert( xupp == INF_POS );

   if(postsolver)
      postsolver->notifyFreeColumnSingletonEquality( row, col, rhs, obj_coeff, col_coeff, xlow, xupp, getSystemMatrix(row.getSystemType()) );

   /* adapt objective from substitution */
   adaptObjectiveSubstitutedRow( row, col, obj_coeff, col_coeff );

   /* remove row and mark column as empty - will be removed in model cleanup on all processes */
   removeRow( row );

   /* remove col and mark it for the fix empty columns presolver */
   ixlow = ixupp = 0.0;
   xlow = INF_NEG;
   xupp = INF_POS;
   if( col.getNode() == -1 )
      outdated_linking_var_bounds = true;

   assert( getNnzsCol(col) == 0 );
   if( col.isLinkingCol() && my_rank == 0 )
      assert( PIPSisZero(getSimpleVecFromColStochVec(*presProb->g, col) + (*objective_vec_chgs)[col.getIndex()]) );
   else if( !col.isLinkingCol() )
      assert( PIPSisZero(getSimpleVecFromColStochVec( *presProb->g, col)) );
}

/* column col getting substituted with row - must be called simultaneously for linking rows */
void PresolveData::adaptObjectiveSubstitutedRow( const INDEX& row, const INDEX& col, double obj_coeff, double col_coeff )
{
   assert(row.isRow());

   if( !col.isEmpty() )
   {
      assert( col.isCol() );
      assert( col.hasValidNode(nChildren) );

      if( col.isLinkingCol() && row.getNode() == -1 && !row.isLinkingRow() )
         assert(PIPS_MPIisValueEqual(row.getIndex()));
   }

   if( col.isEmpty() )
      assert( row.isLinkingRow() );
   if( row.isLinkingRow() )
      assert(PIPS_MPIisValueEqual(row.getIndex()));

   assert( row.hasValidNode(nChildren) );

   const BlockType block_type = row.isLinkingRow() ? BL_MAT : row.getBlockOfColInRow(col);

   if( col.isCol() )
   {
#ifndef NDEBUG
      const SparseStorageDynamic& col_mat_tp = getSparseGenMatrix(row, col)->getStorageDynamicTransposedRef();
#endif
      assert( (col_mat_tp.getRowPtr(col.getIndex()).end - col_mat_tp.getRowPtr(col.getIndex()).start) == 1 );
      assert( row.getIndex() == col_mat_tp.getJcolM(col_mat_tp.getRowPtr(col.getIndex()).start) );
      assert( col_coeff == col_mat_tp.getMat(col_mat_tp.getRowPtr(col.getIndex()).start) );
      assert( obj_coeff == getSimpleVecFromColStochVec( *presProb->g, col) );
   }

   assert( !PIPSisZero(col_coeff) );
   
   const int dummy_index = std::numeric_limits<int>::infinity();

   if( !row.isLinkingRow() )
   {
      assert( block_type == B_MAT || block_type == A_MAT );

      const INDEX col_b_mat(COL, row.getNode(), dummy_index);

      /* Bmat */
      const SparseStorageDynamic& b_mat = getSparseGenMatrix(row, col_b_mat)->getStorageDynamicRef();

      for(int i = b_mat.getRowPtr(row.getIndex()).start ; i < b_mat.getRowPtr(row.getIndex()).end; ++i)
      {
         const int col_idx = b_mat.getJcolM(i);
         if(col_idx != col.getIndex() || block_type == A_MAT)
            getSimpleVecFromColStochVec( *presProb->g, row.getNode())[col_idx] -= obj_coeff * b_mat.getMat(i) / col_coeff;
      }

      /* Amat */
      if( row.getNode() != -1 )
      {
         const INDEX col_a_mat(COL, -1, dummy_index);

         const SparseStorageDynamic& a_mat = getSparseGenMatrix( row, col_a_mat )->getStorageDynamicRef();
         
         for(int i = a_mat.getRowPtr(row.getIndex()).start; i < a_mat.getRowPtr(row.getIndex()).end; ++i)
         {
            const int col_idx = a_mat.getJcolM(i);
            if( col_idx != col.getIndex() || block_type == B_MAT )
            {
               (*objective_vec_chgs)[col_idx] -= obj_coeff * a_mat.getMat(i) / col_coeff;
               outdated_obj_vector = true;
            }
         }
      }
   }
   else
   {
      assert( block_type == BL_MAT );
      const INDEX col_b0_mat(COL, -1, dummy_index);
      /* Bl0 */
      const SparseStorageDynamic& bl0_mat = getSparseGenMatrix( row, col_b0_mat )->getStorageDynamicRef();
      
      for(int i = bl0_mat.getRowPtr(row.getIndex()).start; i < bl0_mat.getRowPtr(row.getIndex()).end; ++i)
      {
         const int col_ptr = bl0_mat.getJcolM(i);
         const double aij = bl0_mat.getMat(i);

         if( col.isEmpty() )
            getSimpleVecFromColStochVec( *presProb->g, -1)[col_ptr] -= obj_coeff * aij / col_coeff;
         else if( col_ptr != col.getIndex() || !col.isLinkingCol() )
            getSimpleVecFromColStochVec( *presProb->g, -1)[col_ptr] -= obj_coeff * aij / col_coeff;
      }

      /* Bl_i */
      for( int node = 0; node < nChildren; ++node )
      {
         if(!nodeIsDummy(node))
         {
            const INDEX col_bli_mat(COL, node, dummy_index);
            const SparseStorageDynamic& bli_mat = getSparseGenMatrix( row, col_bli_mat )->getStorageDynamicRef();

            for(int i = bli_mat.getRowPtr(row.getIndex()).start ; i < bli_mat.getRowPtr(row.getIndex()).end; ++i)
            {
               const int col_ptr = bli_mat.getJcolM(i);
               const double aij = bli_mat.getMat(i);

               if( col.isEmpty() )
                  getSimpleVecFromColStochVec( *presProb->g, node)[col_ptr] -= obj_coeff * aij / col_coeff;
               else if(col_ptr != col.getIndex() || col.getNode() != node)
                  getSimpleVecFromColStochVec( *presProb->g, node)[col_ptr] -= obj_coeff * aij / col_coeff;
            }
         }
      }
   }

   /* rhs/clow == cupp */
   if( row.inInEqSys() )
   {
      assert( PIPSisEQ(getSimpleVecFromRowStochVec( *presProb->bl, row), getSimpleVecFromRowStochVec( *presProb->bu, row)) );
      assert( !PIPSisZero(getSimpleVecFromRowStochVec( *presProb->iclow, row)) );
      assert( !PIPSisZero(getSimpleVecFromRowStochVec( *presProb->icupp, row)) );
   }

   const double rhs = row.inEqSys() ? getSimpleVecFromRowStochVec( *presProb->bA, row ) :
      getSimpleVecFromRowStochVec( *presProb->bu, row );

   /* this is the case where everyone removes the same row in parallel */
   if( col.isCol() )
   {
      if( ( col.isLinkingCol() && row.getNode() == -1) || row.isLinkingRow() )
      {
         if(my_rank == 0)
            obj_offset_chgs += obj_coeff * rhs / col_coeff;
      }
      else
         obj_offset_chgs += obj_coeff * rhs / col_coeff;

   /* if a liking column in a local row gets removed we have to communicate the changes */
   if( col.isCol() && col.isLinkingCol() && row.getNode() != -1 )
   {
      (*objective_vec_chgs)[col.getIndex()] -= getSimpleVecFromColStochVec( *presProb->g, col);
      outdated_obj_vector = true;
   }
   else if( col.isCol() )
      getSimpleVecFromColStochVec( *presProb->g, col) = 0.0;
   }
}

/** adds coeff * col to row */
void PresolveData::addCoeffColToRow( double coeff, const INDEX& col, const INDEX& row )
{
   assert( !PIPSisZero(coeff) );

   assert( row.isRow() );
   assert( col.isCol() );

   assert( !row.isLinkingRow() );

   BlockType block_type = row.getBlockOfColInRow(col);
   const int node = row.isLinkingRow() ? col.getNode() : row.getNode();

   SparseGenMatrix* sparse_mat = getSparseGenMatrixFromStochMat( getSystemMatrix( row.getSystemType() ), node, block_type );
   assert( sparse_mat );

   const SparseStorageDynamic& mat = sparse_mat->getStorageDynamicRef();

   if( postsolver )
   {
      postsolver->notifyColModified(col);
      postsolver->notifyRowModified(row);
   }

   const int row_start = mat.getRowPtr(row.getIndex()).start;
   const int row_end = mat.getRowPtr(row.getIndex()).end;
   const bool new_coeff_in_row = (std::find( mat.getJcolM() + row_start, mat.getJcolM() + row_end, col.getIndex() ) == mat.getJcolM() + row_end);

   /* only increase non-zero counters if coeff has not been in the row before */
   if( new_coeff_in_row )
   {
      const bool at_root = row.getNode() == -1 && col.isLinkingCol();
      increaseNnzCounterColumnBy(col, 1, at_root);
      increaseNnzCounterRowBy(row, 1, at_root);
   }

   sparse_mat->addColToRow( coeff, col.getIndex(), row.getIndex() );
}

/* removes row from local system - sets rhs lhs and activities to zero */
void PresolveData::removeRow( const INDEX& row )
{
   assert( row.isRow() );
   assert( row.hasValidNode(nChildren) );
   assert( !nodeIsDummy( row.getNode() ) );

   const int dummy_index = std::numeric_limits<int>::infinity();
   /* the postsolver must have been notified before the actual removal of the row! */
   if( postsolver )
      assert(wasRowRemoved(row));

   if( row.isLinkingRow() )
   {
      assert( row.getNode() == -1 );

      /* Bl0 */
      const INDEX col_bl0(COL, -1, dummy_index);
      removeRowFromMatrix(row, col_bl0);

      /* linking rows Bli */
      for(int child = 0; child < nChildren; ++child)
      {
         if( !nodeIsDummy(child) )
         {
            const INDEX col_bli(COL, child, dummy_index);
            removeRowFromMatrix(row, col_bli);
         }
      }
   }
   else
   {
      /* Bmat */
      const INDEX col_b_mat(COL, row.getNode(), dummy_index);
      removeRowFromMatrix(row, col_b_mat );

      /* Amat */
      if(row.getNode() != -1)
      {
         const INDEX col_a_mat(COL, -1, dummy_index);
         removeRowFromMatrix(row, col_a_mat);
      }
   }


   /* set lhs rhs to zero */
   if( row.inEqSys() )
      getSimpleVecFromRowStochVec(*presProb->bA, row) = 0.0;
   else
   {
      getSimpleVecFromRowStochVec(*presProb->bl, row) = 0.0;
      getSimpleVecFromRowStochVec(*presProb->bu, row) = 0.0;
   }

   /* set activities and unbounded counters to zero */
   if( row.inEqSys() )
   {
      getSimpleVecFromRowStochVec(*actmax_eq_part, row) = 0.0;
      getSimpleVecFromRowStochVec(*actmin_eq_part, row) = 0.0;
      getSimpleVecFromRowStochVec(*actmax_eq_ubndd, row) = 0;
      getSimpleVecFromRowStochVec(*actmin_eq_ubndd, row) = 0;

      if( row.isLinkingRow() )
      {
         const int row_index = row.getIndex();
         (*actmax_eq_chgs)[row_index] = 0.0;
         (*actmin_eq_chgs)[row_index] = 0.0;
         (*actmax_eq_ubndd_chgs)[row_index] = 0;
         (*actmin_eq_ubndd_chgs)[row_index] = 0;
      }
   }
   else
   {
      getSimpleVecFromRowStochVec(*actmax_ineq_part, row) = 0.0;
      getSimpleVecFromRowStochVec(*actmin_ineq_part, row) = 0.0;
      getSimpleVecFromRowStochVec(*actmax_ineq_ubndd, row) = 0;
      getSimpleVecFromRowStochVec(*actmin_ineq_ubndd, row) = 0;

      if( row.isLinkingRow() )
      {
         const int row_index = row.getIndex();
         (*actmax_ineq_chgs)[row_index] = 0.0;
         (*actmin_ineq_chgs)[row_index] = 0.0;
         (*actmax_ineq_ubndd_chgs)[row_index] = 0;
         (*actmin_ineq_ubndd_chgs)[row_index] = 0;
      }
   }


#ifndef NDEBUG
   /* assert non-zero counters of row are zero - only works for non-linking rows */
   if( row.inEqSys() )
   {
      if( !row.isLinkingRow() )
         assert( getSimpleVecFromRowStochVec(*nnzs_row_A, row) == 0 );
   }
   else
   {
      if( !row.isLinkingRow() )
         assert( getSimpleVecFromRowStochVec(*nnzs_row_C, row) == 0 );
   }  
#endif
}

void PresolveData::removeRowFromMatrix( const INDEX& row, const INDEX& col )
{
   assert( row.isRow() );
   assert( col.isCol() );
   assert( !nodeIsDummy(row.getNode()) );
   assert( !nodeIsDummy(col.getNode()) );

   assert( row.hasValidNode(nChildren) );
   assert( col.hasValidNode(nChildren) );

   assert( col.getIndex() == std::numeric_limits<int>::infinity() );

   SparseGenMatrix* mat = getSparseGenMatrix(row, col);

   assert(mat);
   assert(mat->hasDynamicStorage());

   const SparseStorageDynamic& mat_storage = mat->getStorageDynamicRef();

   /* update non-zero counters and tell postsolver which rows/cols changed */
   const int row_start = mat_storage.getRowPtr(row.getIndex()).start;
   const int row_end = mat_storage.getRowPtr(row.getIndex()).end;

   const bool at_root = ( row.getNode() == -1 && col.isLinkingCol() );
   reduceNnzCounterRowBy(row, row_end - row_start, at_root);

   for(int k = row_start; k < row_end; k++)
   {
      const int col_idx = mat_storage.getJcolM(k);

      const INDEX col_modified(COL, col.getNode(), col_idx);
      if(postsolver)
      {
         postsolver->notifyColModified( col_modified );
         postsolver->notifyRowModified( row );
      }
      reduceNnzCounterColumnBy( col_modified, 1, at_root );
   }

   mat->removeRow( row.getIndex() );
}

bool PresolveData::verifyActivities() const
{
   assert(!outdated_activities && linking_rows_need_act_computation == 0);

   bool activities_correct = true;

   StochVectorHandle actmax_eq_part_new(dynamic_cast<StochVector*>(actmax_eq_part->clone()));
   StochVectorHandle actmin_eq_part_new(dynamic_cast<StochVector*>(actmin_eq_part->clone()));

   StochVectorBaseHandle<int> actmax_eq_ubndd_new(dynamic_cast<StochVectorBase<int>*>(actmax_eq_ubndd->clone()));
   StochVectorBaseHandle<int> actmin_eq_ubndd_new(dynamic_cast<StochVectorBase<int>*>(actmin_eq_ubndd->clone()));

   StochVectorHandle actmax_ineq_part_new(dynamic_cast<StochVector*>(actmax_ineq_part->clone()));
   StochVectorHandle actmin_ineq_part_new(dynamic_cast<StochVector*>(actmin_ineq_part->clone()));

   StochVectorBaseHandle<int> actmax_ineq_ubndd_new(dynamic_cast<StochVectorBase<int>*>(actmax_ineq_ubndd->clone()));
   StochVectorBaseHandle<int> actmin_ineq_ubndd_new(dynamic_cast<StochVectorBase<int>*>(actmin_ineq_ubndd->clone()));

   actmax_eq_part_new->setToZero();
   actmin_eq_part_new->setToZero();

   actmax_eq_ubndd_new->setToZero();
   actmin_eq_ubndd_new->setToZero();

   actmax_ineq_part_new->setToZero();
   actmin_ineq_part_new->setToZero();

   actmax_ineq_ubndd_new->setToZero();
   actmin_ineq_ubndd_new->setToZero();

   recomputeActivities(false, *actmax_eq_part_new, *actmin_eq_part_new, *actmax_eq_ubndd_new, *actmin_eq_ubndd_new, *actmax_ineq_part_new, *actmin_ineq_part_new,
      *actmax_ineq_ubndd_new, *actmin_ineq_ubndd_new);

   const int tracked_rank = 0;
   if( !actmax_eq_part_new->componentEqual(*actmax_eq_part, feastol))
   {
      if(my_rank == tracked_rank)
         std::cout << "on rank " << my_rank << " found actmax_eq_part not correct" << std::endl;
      activities_correct = false;
   }

   if( !actmin_eq_part_new->componentEqual(*actmin_eq_part, feastol))
   {
      if(my_rank == tracked_rank)
         std::cout << "on rank " << my_rank << " found actmin_eq_part not correct" << std::endl;
      activities_correct = false;
   }

   if( !actmax_eq_ubndd_new->componentEqual(*actmax_eq_ubndd, feastol))
   {
      if(my_rank == tracked_rank)
         std::cout << "on rank " << my_rank << " found actmax_eq_ubndd not correct" << std::endl;
      activities_correct = false;
   }

   if( !actmin_eq_ubndd_new->componentEqual(*actmin_eq_ubndd, feastol))
   {
      if(my_rank == tracked_rank)
         std::cout << "on rank " << my_rank << " found actmin_eq_ubndd not correct" << std::endl;
      activities_correct = false;
   }

   if( !actmax_ineq_part_new->componentEqual(*actmax_ineq_part, feastol))
   {
      if(my_rank == tracked_rank)
         std::cout << "on rank " << my_rank << " found actmax_ineq_part not correct" << std::endl;
      activities_correct = false;
   }

   if( !actmin_ineq_part_new->componentEqual(*actmin_ineq_part, feastol))
   {
      if(my_rank == tracked_rank)
         std::cout << "on rank " << my_rank << " found actmin_ineq_part not correct" << std::endl;
      activities_correct = false;
   }

   if( !actmax_ineq_ubndd_new->componentEqual(*actmax_ineq_ubndd, feastol))
   {
      if(my_rank == tracked_rank)
         std::cout << "on rank " << my_rank << " found actmax_ineq_ubndd not correct" << std::endl;
      activities_correct = false;
   }

   if( !actmin_ineq_ubndd_new->componentEqual(*actmin_ineq_ubndd, feastol))
   {
      if(my_rank == tracked_rank)
         std::cout << "on rank " << my_rank << " found actmin_ineq_ubndd not correct" << std::endl;
      activities_correct = false;
   }

   return activities_correct;
}



/** Verifies if the nnzCounters are still correct. */
bool PresolveData::verifyNnzcounters() const
{
   assert(!outdated_nnzs);

   bool nnzCorrect = true;
   StochVectorBaseHandle<int> nnzs_col_new(dynamic_cast<StochVectorBase<int>*>(nnzs_col->cloneFull()));
   StochVectorBaseHandle<int> nnzs_row_A_new(dynamic_cast<StochVectorBase<int>*>(nnzs_row_A->cloneFull()));
   StochVectorBaseHandle<int> nnzs_row_C_new(dynamic_cast<StochVectorBase<int>*>(nnzs_row_C->cloneFull()));

   nnzs_col_new->setToZero();
   nnzs_row_A_new->setToZero();
   nnzs_row_C_new->setToZero();

   initNnzCounter(*nnzs_row_A_new, *nnzs_row_C_new, *nnzs_col_new);

   // linking variables:
   SimpleVectorBase<int>* nColOrigSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_col_new->vec);
   SimpleVectorBase<int>* nColUpdatedSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_col->vec);
   assert( nColUpdatedSimple->n == nColOrigSimple->n );
   for( int i = 0; i < nColUpdatedSimple->n; i++)
   {
      if( (*nColUpdatedSimple)[i] != (*nColOrigSimple)[i])
      {
         std::cout << "Nnz Counter linking column " << i << " not correct: "
               << (*nColUpdatedSimple)[i] << " vs. " << (*nColOrigSimple)[i] << std::endl;
         nnzCorrect = false;
         break;
      }
   }
   // non linking variables:
   for( int it = 0; it < nChildren; it++)
   {
      nColOrigSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_col_new->children[it]->vec);
      nColUpdatedSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_col->children[it]->vec);
      assert( nColUpdatedSimple->n == nColOrigSimple->n );
      for( int i = 0; i < nColUpdatedSimple->n; i++)
      {
         if( (*nColUpdatedSimple)[i] != (*nColOrigSimple)[i])
         {
            std::cout << "Nnz Counter non-linking column " << i << " of child " << it << " not correct: "
                  << (*nColUpdatedSimple)[i] << " vs. " << (*nColOrigSimple)[i] << std::endl;
            nnzCorrect = false;
            break;
         }
      }
   }

   // rows A:
   SimpleVectorBase<int>* nRowAOrigSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_row_A_new->vec);
   SimpleVectorBase<int>* nRowAUpdatedSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_row_A->vec);
   assert( nRowAUpdatedSimple->n == nRowAOrigSimple->n );
   for( int i = 0; i < nRowAUpdatedSimple->n; i++)
   {
      if( (*nRowAUpdatedSimple)[i] != (*nRowAOrigSimple)[i])
      {
         std::cout << "Nnz Counter root A row " << i << " not correct: " << (*nRowAUpdatedSimple)[i] << " vs. "
               << (*nRowAOrigSimple)[i] << std::endl;
         nnzCorrect = false;
         break;
      }
   }
   // child rows:
   for( int it = 0; it < nChildren; it++)
   {
      nRowAOrigSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_row_A_new->children[it]->vec);
      nRowAUpdatedSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_row_A->children[it]->vec);
      assert( nRowAUpdatedSimple->n == nRowAOrigSimple->n );
      for( int i = 0; i < nRowAUpdatedSimple->n; i++)
      {
         if( (*nRowAUpdatedSimple)[i] != (*nRowAOrigSimple)[i])
         {
            std::cout << "Nnz Counter non-linking A row " << i << " of child " << it << " not correct: "
                  << (*nRowAUpdatedSimple)[i] << " vs. " << (*nRowAOrigSimple)[i] << std::endl;
            nnzCorrect = false;
            break;
         }
      }
   }
   if(nnzs_row_A_new->vecl) // linking rows:
   {
      nRowAOrigSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_row_A_new->vecl);
      nRowAUpdatedSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_row_A->vecl);
      assert( nRowAUpdatedSimple->n == nRowAOrigSimple->n );
      for( int i=0; i<nRowAUpdatedSimple->n; i++)
      {
         if( (*nRowAUpdatedSimple)[i] != (*nRowAOrigSimple)[i])
         {
            std::cout << "Nnz Counter linking row of A " << i << " not correct: " << (*nRowAUpdatedSimple)[i] << " vs. "
                  << (*nRowAOrigSimple)[i] << std::endl;
            nnzCorrect = false;
            break;
         }
      }
   }
   // rows C:
   SimpleVectorBase<int>* nRowCOrigSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_row_C_new->vec);
   SimpleVectorBase<int>* nRowCUpdatedSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_row_C->vec);
   assert( nRowCUpdatedSimple->n == nRowCOrigSimple->n );
   for( int i = 0; i < nRowCUpdatedSimple->n; i++)
   {
      if( (*nRowCUpdatedSimple)[i] != (*nRowCOrigSimple)[i])
      {
         std::cout << "Nnz Counter root C row " << i << " not correct: " << (*nRowCUpdatedSimple)[i] << " vs. "
               << (*nRowCOrigSimple)[i] << std::endl;
         nnzCorrect = false;
         break;
      }
   }
   // child rows:
   for( int it = 0; it < nChildren; it++)
   {
      nRowCOrigSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_row_C_new->children[it]->vec);
      nRowCUpdatedSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_row_C->children[it]->vec);
      assert( nRowCUpdatedSimple->n == nRowCOrigSimple->n );
      for( int i = 0; i < nRowCUpdatedSimple->n; i++)
      {
         if( (*nRowCUpdatedSimple)[i] != (*nRowCOrigSimple)[i])
         {
            std::cout << "Nnz Counter non-linking C row " << i << " of child "<< it <<" not correct: "
                  << (*nRowCUpdatedSimple)[i] << " vs. " << (*nRowCOrigSimple)[i] << std::endl;
            nnzCorrect = false;
            break;
         }
      }
   }
   if(nnzs_row_C_new->vecl) // linking rows:
   {
      nRowCOrigSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_row_C_new->vecl);
      nRowCUpdatedSimple = dynamic_cast<SimpleVectorBase<int>*>(nnzs_row_C->vecl);
      assert( nRowCUpdatedSimple->n == nRowCOrigSimple->n );
      for( int i = 0; i < nRowCUpdatedSimple->n; i++)
      {
         if( (*nRowCUpdatedSimple)[i] != (*nRowCOrigSimple)[i])
         {
            std::cout << "Nnz Counter linking row of C " << i << " not correct: " << (*nRowCUpdatedSimple)[i]
                  << " vs. " << (*nRowCOrigSimple)[i] << std::endl;
            nnzCorrect = false;
            break;
         }
      }
   }
   return nnzCorrect;
}

bool PresolveData::nodeIsDummy(int node) const
{
   assert( -1 <= node && node < nChildren );

   if( node == -1 )
      return false;

   assert( getSystemMatrix(EQUALITY_SYSTEM).children[node]->isKindOf(kStochGenDummyMatrix)
		   == getSystemMatrix(INEQUALITY_SYSTEM).children[node]->isKindOf(kStochGenDummyMatrix));

   StochGenMatrix& matrix = getSystemMatrix(EQUALITY_SYSTEM);
   // todo : asserts
   if( matrix.children[node]->isKindOf(kStochGenDummyMatrix))
   {
      assert( dynamic_cast<StochVector&>(*(presProb->bux)).children[node]->isKindOf(kStochDummy) );
      assert( dynamic_cast<StochVector&>(*(presProb->blx)).children[node]->isKindOf(kStochDummy) );

      assert( dynamic_cast<StochVector&>(*(presProb->bA)).children[node]->isKindOf(kStochDummy) );
      assert( dynamic_cast<StochVector&>(*(presProb->bux)).children[node]->isKindOf(kStochDummy) );
      assert( dynamic_cast<StochVector&>(*(presProb->blx)).children[node]->isKindOf(kStochDummy) );
      assert( nnzs_row_A->children[node]->isKindOf(kStochDummy) );
      assert( dynamic_cast<StochVector&>(*(presProb->bu)).children[node]->isKindOf(kStochDummy) );
      assert( dynamic_cast<StochVector&>(*(presProb->bl)).children[node]->isKindOf(kStochDummy) );
      assert( dynamic_cast<StochVector&>(*(presProb->icupp)).children[node]->isKindOf(kStochDummy) );
      assert( dynamic_cast<StochVector&>(*(presProb->iclow)).children[node]->isKindOf(kStochDummy) );
      assert( nnzs_row_C->children[node]->isKindOf(kStochDummy) );

      return true;
   }

   return false;
}

bool PresolveData::hasLinking(SystemType system_type) const
{
   int mlink, nlink;
   const INDEX dummy_link_row(ROW, -1, std::numeric_limits<int>::infinity(), true, system_type);
   const INDEX dummy_link_col(COL, -1, std::numeric_limits<int>::infinity());
   const SparseGenMatrix* mat = getSparseGenMatrix(dummy_link_row, dummy_link_col);
   mat->getSize(mlink, nlink);
   if( mlink > 0 )
   {
      // todo: assert that all vectors and matrices have linking part
      return true;
   }
   return false;
}

/** adjusts unbounded counters of row as well as activity (if applicable) 
 *  assumes col has been fixed to coeff
 */ // todo refactor
void PresolveData::adjustRowActivityFromDeletion( const INDEX& row, const INDEX& col, double coeff )
{
   assert( row.isRow() );
   assert( col.isCol() );
   assert( row.hasValidNode(nChildren) );
   assert( col.hasValidNode(nChildren) );

   assert( !PIPSisZero(coeff) );
   assert( !nodeIsDummy(row.getNode()) );
   assert( !nodeIsDummy(col.getNode()) );

   const bool local_linking_row = (row.isLinkingRow() && !col.isLinkingCol());

   /* get upper and lower bound on variable */
   const double ixlow = getSimpleVecFromColStochVec(*(presProb->ixlow), col);
   const double xlow = getSimpleVecFromColStochVec(*(presProb->blx), col);
   const double ixupp = getSimpleVecFromColStochVec(*(presProb->ixupp), col);
   const double xupp = getSimpleVecFromColStochVec(*(presProb->bux), col);

   /* get unbounded counters */
   int* actmax_ubndd = row.inEqSys() ? &getSimpleVecFromRowStochVec(*actmax_eq_ubndd, row)
         : &getSimpleVecFromRowStochVec(*actmax_ineq_ubndd, row);
   int* actmin_ubndd = row.inEqSys() ? &getSimpleVecFromRowStochVec(*actmin_eq_ubndd, row)
         : &getSimpleVecFromRowStochVec(*actmin_ineq_ubndd, row);

   /* sum of reduction and counters */
   int actmax_ubndd_val = *actmax_ubndd;
   int actmin_ubndd_val = *actmin_ubndd;

   if( local_linking_row )
   {
      actmax_ubndd = row.inEqSys() ? &(*actmax_eq_ubndd_chgs)[row.getIndex()]
            : &(*actmax_ineq_ubndd_chgs)[row.getIndex()];
      actmin_ubndd = row.inEqSys() ? &(*actmin_eq_ubndd_chgs)[row.getIndex()]
            : &(*actmin_ineq_ubndd_chgs)[row.getIndex()];

      actmax_ubndd_val += *actmax_ubndd;
      actmin_ubndd_val += *actmin_ubndd;
   }

   assert(0 <= actmax_ubndd_val);
   assert(0 <= actmin_ubndd_val);

   /* get partial activities */
   double* actmax_part = row.inEqSys() ? &getSimpleVecFromRowStochVec(*actmax_eq_part, row)
         : &getSimpleVecFromRowStochVec(*actmax_ineq_part, row);
   double* actmin_part = row.inEqSys() ? &getSimpleVecFromRowStochVec(*actmin_eq_part, row)
         : &getSimpleVecFromRowStochVec(*actmin_ineq_part, row);

   if( local_linking_row )
   {
      actmax_part = row.inEqSys() ? &(*actmax_eq_chgs)[row.getIndex()] : &(*actmax_ineq_chgs)[row.getIndex()];
      actmin_part = row.inEqSys() ? &(*actmin_eq_chgs)[row.getIndex()] : &(*actmin_ineq_chgs)[row.getIndex()];
   }

   assert(actmin_part);
   assert(actmax_part);
   assert(actmin_ubndd);
   assert(actmax_ubndd);

   /* depending on the sign of coeff add / substract lbx * coeff/ ubx * coeff from actmin and actmax */
   if( PIPSisLT(0.0, coeff) )
   {
      if( PIPSisEQ(ixupp, 1.0) )
         (*actmax_part) -= coeff * xupp;
      else
      {
         assert(1 <= actmax_ubndd_val);
         --(*actmax_ubndd);
         if( actmax_ubndd_val == 2)
            computeRowMinOrMaxActivity(row, true);
      }

      if( PIPSisEQ(ixlow, 1.0) )
         (*actmin_part) -= coeff * xlow;
      else
      {
         assert(1 <= actmin_ubndd_val);
         --(*actmin_ubndd);
         if( actmin_ubndd_val == 2)
            computeRowMinOrMaxActivity(row, false);
      }
   }
   else
   {
      if( PIPSisEQ(ixlow, 1.0) )
         (*actmax_part) -= coeff * xlow;
      else
      {
         assert(1 <= actmax_ubndd_val);
         --(*actmax_ubndd);
         if( actmax_ubndd_val == 2)
            computeRowMinOrMaxActivity(row, true);
      }

      if( PIPSisEQ(ixupp, 1.0) )
         (*actmin_part) -= coeff * xupp;
      else
      {
         assert(1 <= actmin_ubndd_val);
         --(*actmin_ubndd);
         if( actmin_ubndd_val == 2)
            computeRowMinOrMaxActivity(row, false);
      }
   }

   if( local_linking_row )
      outdated_activities = true;
}

/// computes min or max activity of linking row regardless of unbounded counters
double PresolveData::computeLocalLinkingRowMinOrMaxActivity(const INDEX& row, bool upper) const
{
   assert( row.isRow() );
   assert( row.isLinkingRow() );
   double act_part = 0.0;

   const int dummy_index = std::numeric_limits<int>::infinity();
   for( int node = -1; node < nChildren; ++node )
   {
      /* Bl0 block is added only if my_rank == 0 */
      if(node == -1 && my_rank != 0)
         continue;

      if( nodeIsDummy(node) )
         continue;

      /* get upper, lower bounds */
      const SimpleVector& ixlow = getSimpleVecFromColStochVec(*(presProb->ixlow), node);
      const SimpleVector& xlow = getSimpleVecFromColStochVec(*(presProb->blx), node);
      const SimpleVector& ixupp = getSimpleVecFromColStochVec(*(presProb->ixupp), node);
      const SimpleVector& xupp = getSimpleVecFromColStochVec(*(presProb->bux), node);

      /* get matrix */
      const INDEX dummy_col(COL, node, dummy_index);
      const SparseStorageDynamic& mat = getSparseGenMatrix(row, dummy_col)->getStorageDynamicRef();

      for( int j = mat.getRowPtr(row.getIndex()).start; j < mat.getRowPtr(row.getIndex()).end; j++ )
      {
         const int col = mat.getJcolM(j);
         const double entry = mat.getMat(j);

         assert( 0 <= col && col < ixlow.n);
         assert( !PIPSisZero(entry));

         if( PIPSisLT(0.0, entry) )
         {
            if( !upper && !PIPSisZero(ixlow[col]) )
               act_part += entry * xlow[col];
            if( upper && !PIPSisZero(ixupp[col]) )
               act_part += entry * xupp[col];
         }
         else
         {
            if( !upper && !PIPSisZero(ixupp[col]) )
               act_part += entry * xupp[col];
            if( upper && !PIPSisZero(ixlow[col]) )
               act_part += entry * xlow[col];
         }
      }
   }

   return act_part;
}

/// recomputes and updates the activity of a locally available row or increases the recompute counter for linking constraints
/// does not compute activities for linking rows
void PresolveData::computeRowMinOrMaxActivity(const INDEX& row, bool upper)
{
   assert( row.isRow() );
   assert( row.hasValidNode(nChildren) );

   const int dummy_index = std::numeric_limits<int>::infinity();
   
   if( track_row && tracked_row == row )
   {
      if( row.isLinkingRow() )
         std::cout << "TRACKING_ROW: " << row << " is linking constraint and needs " << (upper ? "upper" : "lower") << " activity recomputation " << std::endl;
      else
         std::cout << "TRACKING_ROW: computing " << (upper ? "upper" : "lower") << " activity of " << row << std::endl;
   }

   /* single linking rows that have to be computed for the first time */
   if( row.isLinkingRow() )
   {
      ++linking_rows_need_act_computation;
      return;
   }

   /* upper lower bounds linking vars */
   const SimpleVector& ixlow_root = getSimpleVecFromColStochVec(*(presProb->ixlow), -1);
   const SimpleVector& xlow_root = getSimpleVecFromColStochVec(*(presProb->blx), -1);
   const SimpleVector& ixupp_root = getSimpleVecFromColStochVec(*(presProb->ixupp), -1);
   const SimpleVector& xupp_root = getSimpleVecFromColStochVec(*(presProb->bux), -1);

   /* get upper, lower bounds */
   const SimpleVector& ixlow = getSimpleVecFromColStochVec(*(presProb->ixlow), row.getNode());
   const SimpleVector& xlow = getSimpleVecFromColStochVec(*(presProb->blx), row.getNode());
   const SimpleVector& ixupp = getSimpleVecFromColStochVec(*(presProb->ixupp), row.getNode());
   const SimpleVector& xupp = getSimpleVecFromColStochVec(*(presProb->bux), row.getNode());

   /* get matrix */
   const SparseStorageDynamic& Bmat = getSparseGenMatrix(row, INDEX(COL, row.getNode(), dummy_index) )->getStorageDynamicRef();

   /* get activity vector */
   SimpleVector* act_vec;
   if( row.inEqSys() )
   {
      act_vec = (!upper) ? &getSimpleVecFromRowStochVec(*actmin_eq_part, row.getNode(), row.getLinking()) :
            &getSimpleVecFromRowStochVec(*actmax_eq_part, row.getNode(), row.getLinking());
   }
   else
   {
      act_vec = (!upper) ? &getSimpleVecFromRowStochVec(*actmin_ineq_part, row.getNode(), row.getLinking()) :
            &getSimpleVecFromRowStochVec(*actmax_ineq_part, row.getNode(), row.getLinking());
   }

   double& act_part = (*act_vec)[row.getIndex()];
   act_part = 0;

   /* Bmat */
   for( int j = Bmat.getRowPtr(row.getIndex()).start; j < Bmat.getRowPtr(row.getIndex()).end; j++)
   {
      const int col = Bmat.getJcolM(j);
      const double entry = Bmat.getMat(j);

      assert( 0 <= col && col < ixlow.n );
      assert( !PIPSisZero(entry) );

      if( PIPSisLT(0.0, entry) )
      {
         if( !upper && !PIPSisZero(ixlow[col]) )
            act_part += entry * xlow[col];
         if( upper && !PIPSisZero(ixupp[col]) )
            act_part += entry * xupp[col];
      }
      else
      {
         if( !upper && !PIPSisZero(ixupp[col]) )
            act_part += entry * xupp[col];
         if( upper && !PIPSisZero(ixlow[col]) )
            act_part += entry * xlow[col];
      }
   }

   /* Amat */
   if( row.getNode() != -1 )
   {
      const SparseStorageDynamic& Amat = getSparseGenMatrix(row, INDEX(COL, -1, dummy_index))->getStorageDynamicRef();
      for( int j = Amat.getRowPtr(row.getIndex()).start; j < Amat.getRowPtr(row.getIndex()).end; j++ )
      {
         const int col = Amat.getJcolM(j);
         const double entry = Amat.getMat(j);

         assert( 0 <= col && col < ixlow_root.n);
         assert( !PIPSisZero(entry));

         if( PIPSisLT(0.0, entry) )
         {
            if( !upper && !PIPSisZero(ixlow_root[col]) )
               act_part += entry * xlow_root[col];
            if( upper && !PIPSisZero(ixupp_root[col]) )
               act_part += entry * xupp_root[col];
         }
         else
         {
            if( !upper && !PIPSisZero(ixupp_root[col]) )
               act_part += entry * xupp_root[col];
            if( upper && !PIPSisZero(ixlow_root[col]) )
               act_part += entry * xlow_root[col];
         }
      }
   }

   if( track_row && tracked_row == row )
         std::cout << "TRACKING_ROW: new " << (upper ? "upper" : "lower") << " activity is " << act_part << std::endl;
}

/** goes through given column and adjusts activities by a_ij * (old_ubx - ubx)
 *  If the variable previously was unbounded in upper/lower direction unbounded_counters
 *  get adjusted and the row activity is computed for the first time if necessary
 */
void PresolveData::updateRowActivities(const INDEX& col, double xlow_new, double xupp_new, double xlow_old, double xupp_old)
{
   assert( col.isCol() );
   assert( col.hasValidNode(nChildren) );
   assert(xlow_new != INF_NEG || xupp_new != INF_POS);

   const int dummy_index = std::numeric_limits<int>::infinity();

   if( track_col && tracked_col == col )
   {
      std::cout << "TRACKING_COLUMN: " << col << " gets it's activities updated" << std::endl;
      std::cout << "\t bounds changed from [" << xlow_old << ", " << xupp_old << "] to [" << xlow_new << ", " << xupp_new << "]" << std::endl;
   }

   /* if node == -1 go through all linking var blocks of both systems */
   if( col.isLinkingCol() )
   {
      if( track_col && tracked_col == col )
         std::cout << "TRACKING_COLUMN: the column is linking (root node)" << std::endl;
      /* B0 and D0block */
      updateRowActivitiesBlock(INDEX(ROW, -1, dummy_index, false, EQUALITY_SYSTEM), col, xlow_new, xupp_new, xlow_old, xupp_old);
      updateRowActivitiesBlock(INDEX(ROW, -1, dummy_index, false, INEQUALITY_SYSTEM), col, xlow_new, xupp_new, xlow_old, xupp_old);

      /* Bl0 and Dl0 */
      updateRowActivitiesBlock(INDEX(ROW, -1, dummy_index, true, EQUALITY_SYSTEM), col, xlow_new, xupp_new, xlow_old, xupp_old);
      updateRowActivitiesBlock(INDEX(ROW, -1, dummy_index, true, INEQUALITY_SYSTEM), col, xlow_new, xupp_new, xlow_old, xupp_old);

      for(int child = 0; child < nChildren; ++child)
      {
         if( !nodeIsDummy(child) )
         {
            /* Ai and Ci */
            updateRowActivitiesBlock( INDEX(ROW, child, dummy_index, false, EQUALITY_SYSTEM), col, xlow_new, xupp_new, xlow_old, xupp_old);
            updateRowActivitiesBlock( INDEX(ROW, child, dummy_index, false, INEQUALITY_SYSTEM), col, xlow_new, xupp_new, xlow_old, xupp_old);
         }
      }
   }
   else
   {
      if( track_col && tracked_col == col )
         std::cout << "TRACKING_COLUMN: the column is non-linking (non-root)" << std::endl;
      /* Bmat, Blmat */
      /* Bmat */
      updateRowActivitiesBlock(INDEX(ROW, col.getNode(), dummy_index, false, EQUALITY_SYSTEM), col, xlow_new, xupp_new, xlow_old, xupp_old);

      /* Blmat */
      updateRowActivitiesBlock(INDEX(ROW, -1, dummy_index, true, EQUALITY_SYSTEM), col, xlow_new, xupp_new, xlow_old, xupp_old);

      /* Dmat Dlmat */

      /* Dmat */
      updateRowActivitiesBlock(INDEX(ROW, col.getNode(), dummy_index, false, INEQUALITY_SYSTEM), col, xlow_new, xupp_new, xlow_old, xupp_old);

      /* Dlmat */
      updateRowActivitiesBlock(INDEX(ROW, -1, dummy_index, true, INEQUALITY_SYSTEM), col, xlow_new, xupp_new, xlow_old, xupp_old);
   }
}


void PresolveData::updateRowActivitiesBlock( const INDEX& row, const INDEX& col, double xlow_new, double xupp_new, double xlow_old, double xupp_old)
{
   assert( row.isRow() );
   assert( col.isCol() );
   assert( row.getIndex() == std::numeric_limits<int>::infinity() );

   if( xlow_new != INF_NEG && PIPSisLT( xlow_old, xlow_new) )
      updateRowActivitiesBlock(row, col, xlow_new, xlow_old, false);
	if( xupp_new != INF_POS && PIPSisLT( xupp_new, xupp_old) )
	   updateRowActivitiesBlock(row, col, xupp_new, xupp_old, true);
}


void PresolveData::updateRowActivitiesBlock(const INDEX& row, const INDEX& col, double bound, double old_bound, bool upper)
{
   assert( row.isRow() );
   assert( col.isCol() );

   assert( row.getIndex() == std::numeric_limits<int>::infinity() );

   assert( row.hasValidNode(nChildren) );
   assert( col.hasValidNode(nChildren) );

   /* dummies do not adjust anything */
   if( nodeIsDummy(row.getNode()) )
      return;

   /* we do not have to adjust activities if no new bound was found */
   if( bound == INF_NEG || bound == INF_POS )
      return;

   /* no changes if a worse bound has been found */
   if( (upper && !PIPSisLT( bound, old_bound)) || (!upper && PIPSisLT(bound, old_bound)) )
      return;

   const SparseStorageDynamic& mat_transp = getSparseGenMatrix(row, col)->getStorageDynamicTransposedRef();

   assert(col.getIndex() < mat_transp.getM() );

   SimpleVectorBase<int>& actmax_ubndd = row.inEqSys() ? getSimpleVecFromRowStochVec(*actmax_eq_ubndd, row.getNode(), row.isLinkingRow()) :
         getSimpleVecFromRowStochVec(*actmax_ineq_ubndd, row.getNode(), row.isLinkingRow());
   SimpleVectorBase<int>& actmin_ubndd = row.inEqSys() ? getSimpleVecFromRowStochVec(*actmin_eq_ubndd, row.getNode(), row.isLinkingRow()) :
         getSimpleVecFromRowStochVec(*actmin_ineq_ubndd, row.getNode(), row.isLinkingRow());
   SimpleVector& actmax_part = row.inEqSys() ? getSimpleVecFromRowStochVec(*actmax_eq_part, row.getNode(), row.isLinkingRow()) :
         getSimpleVecFromRowStochVec(*actmax_ineq_part, row.getNode(), row.isLinkingRow());
   SimpleVector& actmin_part = row.inEqSys() ? getSimpleVecFromRowStochVec(*actmin_eq_part, row.getNode(), row.isLinkingRow()) :
         getSimpleVecFromRowStochVec(*actmin_ineq_part, row.getNode(), row.isLinkingRow());

   /* we always set these variables but only use them in case an actual linking row gets it's activities updated */
   SimpleVectorBase<int>& actmax_ubndd_chgs = row.inEqSys() ? *actmax_eq_ubndd_chgs : *actmax_ineq_ubndd_chgs;
   SimpleVectorBase<int>& actmin_ubndd_chgs = row.inEqSys() ? *actmin_eq_ubndd_chgs : *actmin_ineq_ubndd_chgs;
   SimpleVector& actmax_part_chgs = row.inEqSys() ? *actmax_eq_chgs : *actmax_ineq_chgs;
   SimpleVector& actmin_part_chgs = row.inEqSys() ? *actmin_eq_chgs : *actmin_ineq_chgs;

   if( track_col && tracked_col == col )
      std::cout << "TRACKING_COLUMN: updating activities of " << col << " in " << row <<
         " with new " << ( upper ? "upper" : "lower" ) << " bound " << bound << " from " << old_bound << std::endl;

   for( int j = mat_transp.getRowPtr(col.getIndex()).start; j < mat_transp.getRowPtr(col.getIndex()).end; ++j )
   {
      const int row_index = mat_transp.getJcolM(j);
      const INDEX row_INDEX(ROW, row.getNode(), row_index, row.getLinking(), row.getSystemType());
      const double entry = mat_transp.getMat(j);

      assert( !PIPSisZero(entry) );

      /* get affected partial activity and act_ubndd */
      bool switch_upperlower = upper;
      if( PIPSisLT(entry, 0.0) )
         switch_upperlower = !upper;

      SimpleVector& act_part = (switch_upperlower) ? actmax_part : actmin_part;
      SimpleVector& act_part_chgs = (switch_upperlower) ? actmax_part_chgs : actmin_part_chgs;
      SimpleVectorBase<int>& act_ubndd = (switch_upperlower) ? actmax_ubndd : actmin_ubndd;
      SimpleVectorBase<int>& act_ubndd_chgs = (switch_upperlower) ? actmax_ubndd_chgs : actmin_ubndd_chgs;

      if( track_row && tracked_row == row_INDEX )
      {
         std::cout << "TRACKING_ROW: activities of " << row_INDEX << " are getting updated. Current state:" << std::endl;
         if( row_INDEX.isLinkingRow() )
         {
            std::cout << "\twas linking row with " << ( switch_upperlower ? "upper" : "lower" ) << " activity " << act_part[row_index] << " + " << act_part_chgs[row_index] << "(changes) and " <<
                  (switch_upperlower ? "upper" : "lower" ) << " unbounded counters " << act_ubndd[row_index] << " + " << act_ubndd_chgs[row_index] << "(changes)" << std::endl;
         }
         else
         {
            std::cout << "\twas non-linking row with " << ( switch_upperlower ? "upper" : "lower" ) << " activity " << act_part[row_index] << " and " << (switch_upperlower ? "upper" : "lower" ) <<
                  " unbounded counters " << act_ubndd[row_index] << std::endl;
         }
      }

      /* if the old bound was not set we have to modify the unbounded counters */
      if( old_bound == INF_NEG || old_bound == INF_POS )
      {
         if( track_row && tracked_row == row_INDEX )
               std::cout << "TRACKING_ROW: ubndd counters are being changed" << std::endl;

         /* every process works the root node - even in the linking row case */
         if( !col.isLinkingCol() && row_INDEX.isLinkingRow() )
         {
            assert( act_ubndd[row_index] + act_ubndd_chgs[row_index] > 0 );
            --act_ubndd_chgs[row_index];
            outdated_activities = true;
         }
         else
         {
            if( row_INDEX.isLinkingRow() )
               assert( act_ubndd[row_index] + act_ubndd_chgs[row_index] > 0 );
            else
               assert( act_ubndd[row_index] > 0 );

            --act_ubndd[row_index];
         }

         /* from now on activity of row has to be computed */
         if( ( row_INDEX.isLinkingRow()  && (act_ubndd[row_index] + act_ubndd_chgs[row_index] == 1) )
               || (act_ubndd[row_index] == 1 && !row_INDEX.isLinkingRow()) )
         {

            if( track_row && tracked_row == row_INDEX )
            {
               if( row_INDEX.isLinkingRow()  )
                  std::cout << "TRACKING_ROW: " << ( switch_upperlower ? "upper" : "lower" ) << " now needs act computation " << std::endl;
               else
                  std::cout << "TRACKING_ROW: first time computation of " << ( switch_upperlower ? "upper" : "lower" ) << " activity " << std::endl;
            }

            if( row_INDEX.isLinkingRow() )
               ++linking_rows_need_act_computation;
            else
               computeRowMinOrMaxActivity( row_INDEX, switch_upperlower);

         }
         else if( (row_INDEX.isLinkingRow() && act_ubndd[row_index] + act_ubndd_chgs[row_index] == 0) ||
               (!row_INDEX.isLinkingRow() && act_ubndd[row_index] == 0) )
         {
            if( track_row && tracked_row == row_INDEX )
                  std::cout << "TRACKING_ROW: adjusting activities " << std::endl;

            if( !col.isLinkingCol() && row_INDEX.isLinkingRow() )
            {
               outdated_activities = true;
               act_part_chgs[row_index] += bound * entry;
            }
            else
               act_part[row_index] += bound * entry;
         }
      }
      else
      {
         assert( old_bound != INF_NEG && old_bound != INF_POS );

         if( track_row && tracked_row == row_INDEX )
            std::cout << "TRACKING_ROW: " << " no changes in ubndd counters - only adjust activities" << std::endl;

         /* better bound was found */
         if( (row_INDEX.isLinkingRow() && act_ubndd[row_index] + act_ubndd_chgs[row_index] <= 1) ||
               (!row_INDEX.isLinkingRow() && act_ubndd[row_index] <= 1) )
         {

            if( !col.isLinkingCol() && row_INDEX.isLinkingRow() )
            {
               outdated_activities = true;
               act_part_chgs[row_index] += (bound - old_bound) * entry;
            }
            else
               act_part[row_index] += (bound - old_bound) * entry;
         }
      }

      if( track_row && tracked_row == row_INDEX )
      {
         std::cout << "TRACKING_ROW: activities of tracked row have been updated: " << std::endl;
         if( row_INDEX.isLinkingRow() )
         {
            std::cout << "\tnow linking row with " << ( switch_upperlower ? "upper" : "lower" ) <<
               " activity " << act_part[row_index] << " + " << act_part_chgs[row_index] << "(changes) and " <<
               (switch_upperlower ? "upper" : "lower" ) << " unbounded counters " << 
               act_ubndd[row_index] << " + " << act_ubndd_chgs[row_index] << "(changes)" << std::endl;
         }
         else
         {
            std::cout << "\tnow non-linking row with " << ( switch_upperlower ? "upper" : "lower" ) <<
            " activity " << act_part[row_index] << " and " << (switch_upperlower ? "upper" : "lower" ) <<
            " unbounded counters " << act_ubndd[row_index] << std::endl;
         }
      }
   }
}

/// returns known activity and unbounded counters for specified row
/// +/- infinity if there are two or more unbouded entries in a row 
void PresolveData::getRowActivities( const INDEX& row, double& max_act,
      double& min_act, int& max_ubndd, int& min_ubndd) const
{
   assert(row.isRow());

   assert( !nodeIsDummy( row.getNode() ) );

   max_ubndd = row.inEqSys() ? getSimpleVecFromRowStochVec(*actmax_eq_ubndd, row)
         : getSimpleVecFromRowStochVec(*actmax_ineq_ubndd, row);
   min_ubndd = row.inEqSys() ? getSimpleVecFromRowStochVec(*actmin_eq_ubndd, row)
         : getSimpleVecFromRowStochVec(*actmin_ineq_ubndd, row);

   if( row.isLinkingRow() )
   {
      max_ubndd += row.inEqSys() ? (*actmax_eq_ubndd_chgs)[row.getIndex()]
            : (*actmax_ineq_ubndd_chgs)[row.getIndex()];
      min_ubndd += row.inEqSys() ? (*actmin_eq_ubndd_chgs)[row.getIndex()]
            : (*actmin_ineq_ubndd_chgs)[row.getIndex()];
   }

   max_act = row.inEqSys() ? getSimpleVecFromRowStochVec(*actmax_eq_part, row)
         : getSimpleVecFromRowStochVec(*actmax_ineq_part, row);
   min_act = row.inEqSys() ? getSimpleVecFromRowStochVec(*actmin_eq_part, row)
         : getSimpleVecFromRowStochVec(*actmin_ineq_part, row);

   if( row.isLinkingRow() )
   {
      max_act += row.inEqSys() ? (*actmax_eq_chgs)[row.getIndex()] : (*actmax_ineq_chgs)[row.getIndex()];
      min_act += row.inEqSys() ? (*actmin_eq_chgs)[row.getIndex()] : (*actmin_ineq_chgs)[row.getIndex()];
   }

   if( max_ubndd >= 2)
      assert( max_act == INF_POS );
   else
   {
      if( !row.isLinkingRow() )
      {
         assert(max_act < INF_POS);
      }
   }

   if( min_ubndd >= 2)
      assert( min_act == INF_NEG );
   else
   {
      if( !row.isLinkingRow() )
      {
         assert( INF_NEG < min_act);
      }
   }
}

void PresolveData::getRowBounds( const INDEX& row, double& lhs, double& rhs) const
{
   assert( row.isRow() );

   if( row.inEqSys() )
   {
      rhs = getSimpleVecFromRowStochVec(presProb->bA, row);
      lhs = getSimpleVecFromRowStochVec(presProb->bA, row);

      if( row.isLinkingRow() && outdated_lhsrhs )
      {
         rhs += (*bound_chgs_A)[row.getIndex()];
         lhs += (*bound_chgs_A)[row.getIndex()];
      }

      assert( std::fabs(rhs) != INF_POS );
   }
   else
   {
      rhs = getSimpleVecFromRowStochVec(presProb->bu, row);
      lhs = getSimpleVecFromRowStochVec(presProb->bl, row);

      if( row.isLinkingRow() && outdated_lhsrhs )
      {
         rhs += (*bound_chgs_C)[row.getIndex()];
         lhs += (*bound_chgs_C)[row.getIndex()];
      }

#ifndef NDEBUG
      if( PIPSisZero(getSimpleVecFromRowStochVec(presProb->iclow, row)) )
         assert( lhs == INF_NEG );
      else
         assert( std::fabs(lhs) != INF_POS );

      if( PIPSisZero(getSimpleVecFromRowStochVec(presProb->icupp, row)) )
         assert( rhs == INF_POS );
      else
         assert( std::fabs(rhs) != INF_POS );
#endif
   }
}

void PresolveData::setRowBounds( const INDEX& row, double clow, double cupp )
{
   assert( row.isRow() );
   assert( PIPSisLE( clow, cupp ) );
   assert( clow != INF_NEG || cupp != INF_POS );

   if( row.inEqSys() )
   {
      assert( clow == cupp );
      if( row.isLinkingRow() )
      {
         const int row_index = row.getIndex();
         const double old_clow = getSimpleVecFromRowStochVec(presProb->bA, row);
         double& clow_chgs = (*bound_chgs_A)[row_index];

         clow_chgs = clow - old_clow;

         outdated_lhsrhs = true;
         assert( PIPSisEQ(getSimpleVecFromRowStochVec(presProb->bA, row) + clow_chgs, clow) );
      }
      else
         getSimpleVecFromRowStochVec(presProb->bA, row) = clow;
   }
   else
   {
      if( row.isLinkingRow() )
      {
         const int row_index = row.getIndex();
         // TODO : this needs tow arrays for synchronization if it is supposed to work..
         if( clow != INF_NEG )
         {
            const double old_clow = getSimpleVecFromRowStochVec(presProb->bl, row);
            double& clow_chgs = (*bound_chgs_C)[row_index];

            if( old_clow == INF_NEG )
            {
               getSimpleVecFromRowStochVec(presProb->iclow, row) = 1;
               clow_chgs = clow;
            }
            else
               clow_chgs = clow - old_clow;

            outdated_lhsrhs = true;
         }

         if( cupp != INF_POS )
         {
            const double old_cupp = getSimpleVecFromRowStochVec(presProb->bu, row);
            double& cupp_chgs = (*bound_chgs_C)[row_index];

            if( old_cupp == INF_NEG )
            {
               getSimpleVecFromRowStochVec(presProb->icupp, row) = 1;
               cupp_chgs = cupp;
            }
            else
               cupp_chgs = cupp - old_cupp;

            outdated_lhsrhs = true;
         }

         if( getSimpleVecFromRowStochVec(presProb->bl, row) != INF_NEG )
            assert( PIPSisEQ(getSimpleVecFromRowStochVec(presProb->bl, row) + (*bound_chgs_C)[row_index], clow) );
         else
            assert( PIPSisEQ((*bound_chgs_C)[row_index], clow) );

         if( getSimpleVecFromRowStochVec(presProb->bu, row) != INF_POS )
            assert( PIPSisEQ(getSimpleVecFromRowStochVec(presProb->bl, row) + (*bound_chgs_C)[row_index], clow) );
         else
            assert( PIPSisEQ((*bound_chgs_C)[row_index], cupp) );
      }
      else
      {
         if( clow != INF_NEG )
         {
            getSimpleVecFromRowStochVec(presProb->iclow, row) = 1;
            getSimpleVecFromRowStochVec(presProb->bl, row) = clow;
         }

         if( cupp != INF_POS )
         {
            getSimpleVecFromRowStochVec(presProb->icupp, row) = 1;
            getSimpleVecFromRowStochVec(presProb->bu, row) = cupp;
         }
      }
   }
}

void PresolveData::getColBounds( const INDEX& col, double& xlow, double& xupp) const
{
   assert( col.isCol() );

   xlow = getSimpleVecFromColStochVec(presProb->blx, col);
   xupp = getSimpleVecFromColStochVec(presProb->bux, col);

#ifndef NDBEUG
   if( PIPSisZero(getSimpleVecFromColStochVec(presProb->ixupp, col)) )
      assert( xupp == INF_POS );
   else
      assert( std::fabs(xupp) != INF_POS );
   if( PIPSisZero(getSimpleVecFromColStochVec(presProb->ixlow, col)) )
      assert( xlow == INF_NEG );
   else
      assert( std::fabs(xlow) != INF_POS );
#endif
}

/** updates the bounds on a variable as well as activities */
bool PresolveData::updateColBounds( const INDEX& col, double xlow, double xupp)
{
   assert( xlow != INF_NEG || xupp != INF_POS );
   assert( col.isCol() );

   double xlow_old, xupp_old;
   getColBounds( col, xlow_old, xupp_old );

   if( track_col && tracked_col == col )
      std::cout << "TRACKING_COLUMN: updating column bounds from [" << xlow << ", " << xupp << "] for " << col << " with [" <<
         xlow << ", " << xupp << "]" << std::endl;

   checkBoundsInfeasible(col, xlow, xupp);

   bool updated = false;
   if( xupp < INF_POS && PIPSisLT(xupp, xupp_old) )
   {
      updated = true;

      getSimpleVecFromColStochVec(presProb->bux, col) = xupp;
      getSimpleVecFromColStochVec(presProb->ixupp, col) = 1.0;
   }

   if( INF_NEG < xlow && PIPSisLT( xlow_old, xlow ) )
   {
      updated = true;

      getSimpleVecFromColStochVec(presProb->blx, col) = xlow;
      getSimpleVecFromColStochVec(presProb->ixlow, col) = 1.0;
   }

   if( updated )
   {
      if(  track_col && tracked_col == col )
      {
         std::cout << "TRACKING_COLUMN: bounds are now [" << xlow << ", " << xupp << "]" << std::endl;
         std::cout << "TRACKING_COLUMN: moving on to update activities" << std::endl;
      }
      updateRowActivities(col, xlow, xupp, xlow_old, xupp_old);
   }
   else
   {
      if( track_col && tracked_col == col )
         std::cout << "TRACKING_COLUMN: " << col << " was not updated" << std::endl;
   }

#ifndef NDEBUG
   double xlow_new, xupp_new;
   getColBounds(col, xlow_new, xupp_new);

   if( xlow != INF_NEG )
      assert( xlow_new != INF_NEG );
   if( xupp != INF_POS )
      assert( xupp_new != INF_POS );
#endif

   if( updated && col.isLinkingCol() )
      outdated_linking_var_bounds = true;

   return updated;
}


double PresolveData::getRowCoeff( const INDEX& row, const INDEX& col ) const
{
   assert( row.isRow() );
   assert( col.isCol() );
   assert( !nodeIsDummy(col.getNode()) );

   const SparseStorageDynamic& mat = getSparseGenMatrix(row, col)->getStorageDynamicRef();

   const int row_start = mat.getRowPtr(row.getIndex()).start;
   const int row_end = mat.getRowPtr(row.getIndex()).end;

   assert(row_start != row_end);

   for(int i = row_start; i < row_end; ++i)
   {
      if(mat.getJcolM(i) == col.getIndex())
         return mat.getMat(i);
   }
   assert(false && "could not find coefficient");
   return 0.0;
}


StochGenMatrix& PresolveData::getSystemMatrix(SystemType system_type) const
{
   if( system_type == EQUALITY_SYSTEM )
      return dynamic_cast<StochGenMatrix&>(*presProb->A);
   else
   {
      assert(system_type == INEQUALITY_SYSTEM);
      return dynamic_cast<StochGenMatrix&>(*presProb->C);
   }
}

SparseGenMatrix* PresolveData::getSparseGenMatrix(const INDEX& row, const INDEX& col) const
{
   assert( row.isRow() );
   assert( col.isCol() );
   assert( row.hasValidNode(nChildren) );
   assert( col.hasValidNode(nChildren) );


   const int mat_node = row.isLinkingRow() ? col.getNode() : row.getNode();
   const BlockType block_type = row.getBlockOfColInRow(col);
   assert(!nodeIsDummy(mat_node));

   const StochGenMatrix& sMat = getSystemMatrix(row.getSystemType());
   SparseGenMatrix* res = getSparseGenMatrixFromStochMat(sMat, mat_node, block_type);

   return res;
}

/* only prints the part of a linking constraint the current process knows about */
// todo so far only prints linking constraints on MPI proc 1
void PresolveData::writeRowLocalToStreamDense(std::ostream& out, const INDEX& row) const
{
   assert( row.isRow() );

   const int node = row.getNode();
   const int row_index = row.getIndex();
   const SystemType system_type = row.getSystemType();
   const bool linking = row.getLinking();

   if(nodeIsDummy(node))
      return;

   if( node == -1 && my_rank != 0)
      return;

   out << "SystemType: " << system_type << "\tnode: " << node << "\tLinkingCons: " << linking << "\trow: " << row_index << std::endl;

   if( row.inEqSys() )
      out << getSimpleVecFromRowStochVec(*presProb->bA, row) << " = ";
   else
   {
      double iclow = getSimpleVecFromRowStochVec(*presProb->iclow, row);
      double clow = PIPSisEQ(iclow, 1.0) ? getSimpleVecFromRowStochVec(*presProb->bl, row) : INF_NEG;

      out << clow << " <= ";
   }

   const int dummy_index = std::numeric_limits<int>::infinity();
   if( !row.isLinkingRow() )
   {
      if(node != -1)
      {
         writeMatrixRowToStreamDense(out, *getSparseGenMatrix(row, INDEX(COL, -1, dummy_index)), -1, row_index, getSimpleVecFromColStochVec(*presProb->ixupp, -1),
            getSimpleVecFromColStochVec(*presProb->bux, -1), getSimpleVecFromColStochVec(*presProb->ixlow, -1), getSimpleVecFromColStochVec(*presProb->blx, -1));
      }

      writeMatrixRowToStreamDense(out, *getSparseGenMatrix(row, INDEX(COL, row.getNode(), dummy_index)), node, row_index, getSimpleVecFromColStochVec(*presProb->ixupp, node),
         getSimpleVecFromColStochVec(*presProb->bux, node),getSimpleVecFromColStochVec(*presProb->ixlow, node),getSimpleVecFromColStochVec(*presProb->blx, node));
   }
   else if( row.isLinkingRow() )
   {
      assert(node == -1);
      writeMatrixRowToStreamDense(out, *getSparseGenMatrix(row, INDEX(COL, row.getNode(), dummy_index)), node, row_index, getSimpleVecFromColStochVec(*presProb->ixupp, node),
         getSimpleVecFromColStochVec(*presProb->bux, node),getSimpleVecFromColStochVec(*presProb->ixlow, node),getSimpleVecFromColStochVec(*presProb->blx, node));

      for(int child = 0; child < nChildren; ++child)
      {
         if(!nodeIsDummy(child))
            writeMatrixRowToStreamDense(out, *getSparseGenMatrix(row, INDEX(COL, row.getNode(), dummy_index)), child, row_index, getSimpleVecFromColStochVec(*presProb->ixupp, child),
                  getSimpleVecFromColStochVec(*presProb->bux, child),getSimpleVecFromColStochVec(*presProb->ixlow, child),getSimpleVecFromColStochVec(*presProb->blx, child));
      }
   }

   if(system_type == INEQUALITY_SYSTEM)
   {
      double icupp = getSimpleVecFromRowStochVec(*presProb->icupp, row);
      double cupp = PIPSisEQ(icupp, 1.0) ? getSimpleVecFromRowStochVec(*presProb->bu, row) : INF_POS;

      out << " <= " << cupp;
   }
   out << std::endl;
}

void PresolveData::writeMatrixRowToStreamDense(std::ostream& out, const SparseGenMatrix& mat, int node, int row, const SimpleVector& ixupp, const SimpleVector& xupp,
      const SimpleVector& ixlow, const SimpleVector& xlow) const
{
   const SparseStorageDynamic& storage = mat.getStorageDynamicRef();

   const int start = storage.getRowPtr(row).start;
   const int end = storage.getRowPtr(row).end;

   for(int k = start; k < end; ++k)
   {
      const double val = storage.getMat(k);
      const int col = storage.getJcolM(k);

      out << " + " << val << " * x_" << node << "_" << col << " € [" << ( PIPSisZero(ixlow[col]) ? -std::numeric_limits<double>::infinity() : xlow[col]) << ", "
            << ( PIPSisZero(ixupp[col]) ? std::numeric_limits<double>::infinity() : xupp[col]) << "]";
   }
}

void PresolveData::printVarBoundStatistics(std::ostream& out) const
{
   const StochVector& xupp = dynamic_cast<const StochVector&>(*presProb->bux);
   const StochVector& xlow = dynamic_cast<const StochVector&>(*presProb->blx);
   const StochVector& ixupp = dynamic_cast<const StochVector&>(*presProb->ixupp);
   const StochVector& ixlow = dynamic_cast<const StochVector&>(*presProb->ixlow);

   StochVectorHandle xlow_def = StochVectorHandle(dynamic_cast<StochVector*>(xlow.cloneFull()));
   StochVectorHandle xupp_def = StochVectorHandle(dynamic_cast<StochVector*>(xupp.cloneFull()));

   xlow_def->componentMult(ixlow);
   xupp_def->componentMult(ixupp);

   double nr_bounds_upper = ixupp.dotProductSelf(1.0);
   double nr_bounds_lower = ixlow.dotProductSelf(1.0);

   double xupp_twonorm = xupp_def->twonorm();
   double xupp_infnorm = xupp_def->infnorm();
   double xupp_onenorm = xupp_def->onenorm();

   double xlow_twonorm = xlow_def->twonorm();
   double xlow_infnorm = xlow_def->infnorm();
   double xlow_onenorm = xlow_def->onenorm();

   double min_low;
   double max_low;
   double min_upp;
   double max_upp;
   int a;

   xupp_def->min(min_upp, a);
   xupp_def->max(max_upp, a);
   xlow_def->min(min_low, a);
   xlow_def->max(max_low, a);

   if(distributed && my_rank == 0)
   {
      std::cout << "_____________________________________________________" << std::endl;
      std::cout << "____________________Stats__Bounds____________________" << std::endl;
      std::cout << "_____________________________________________________" << std::endl;

      std::cout << "xlow:" << std::endl;
      std::cout << "nr_bounds\t" << nr_bounds_lower << std::endl;
      std::cout << "twonorm  \t" << xlow_twonorm << std::endl;
      std::cout << "onenorm  \t" << xlow_onenorm << std::endl;
      std::cout << "infnorm  \t" << xlow_infnorm << std::endl;
      std::cout << "min      \t" << max_low << std::endl; 
      std::cout << "max      \t" << max_low << std::endl;
      std::cout << std::endl;
      std::cout << "xupp:" << std::endl;
      std::cout << "nr_bounds\t" << nr_bounds_upper << std::endl;
      std::cout << "twonorm  \t" << xupp_twonorm << std::endl;
      std::cout << "onenorm  \t" << xupp_onenorm << std::endl;
      std::cout << "infnorm  \t" << xupp_infnorm << std::endl;
      std::cout << "min      \t" << min_upp << std::endl;
      std::cout << "max      \t" << max_upp << std::endl;

      std::cout << "_____________________________________________________" << std::endl;
      std::cout << "_____________________________________________________" << std::endl;
   }
}

void PresolveData::printRowColStats() const
{
   const int nnzs_cols = nnzs_col->onenorm();

#ifndef NDEBUG
   const int nnzs_A = nnzs_row_A->onenorm();
   const int nnzs_C = nnzs_row_C->onenorm();

   assert( nnzs_cols == nnzs_A + nnzs_C );
#endif

   const int cols = nnzs_col->getNnzs();
   const int rows_A = nnzs_row_A->getNnzs();
   const int rows_C = nnzs_row_C->getNnzs();

   if( my_rank == 0 )
   {
      std::cout << "Problem Dimensions:" << std::endl;
      std::cout << nnzs_cols << " non-zeros" << std::endl;
      std::cout << cols << " cols\t" << rows_A << " rows A\t" << rows_C << " rows C " << std::endl;
   }
}
