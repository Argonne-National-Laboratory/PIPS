/*
 * PresolveData.C
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

// todo : preallocate singleton rows and singleton columns?

#include "PresolveData.h"
#include "StochGenMatrix.h"
#include "DoubleMatrixTypes.h"
#include "StochVectorUtilities.h"
#include "pipsdef.h"

#include <stdexcept>
#include <limits>
#include <algorithm>
#include <string>
#include <cmath>

/* can specify a cloumn here which bound changes we wanna track maybe */

// #ifndef NDEBUG
//   #define TRACK_C
//   #define COLUMN 100
//   #define COL_NODE 1
// #endif

// #ifndef NDEBUG
//    #define TRACK_R
//    #define ROW 42
//    #define ROW_NODE 2
//    #define ROW_BLOCK BL_MAT
//    #define ROW_IS_LINK false
//    #define ROW_SYS EQUALITY_SYSTEM
// #endif

#ifdef TRACK_C
#define TRACK_COLUMN(node, col)                    \
   (COL_NODE == node && PIPSisEQ(COLUMN, col)      \
      && (my_rank == 0 || node != -1)              \
      && !nodeIsDummy(COL_NODE))
#else
#define TRACK_COLUMN(node, col) false
#endif

#ifdef TRACK_R
#define TRACK_ROW(node, row, sys, link)            \
   (PIPSisEQ(row, ROW) && node == ROW_NODE         \
      && sys == ROW_SYS                            \
      && !nodeIsDummy(ROW_NODE)           \
      && (my_rank == 0 || link ||                  \
         ROW_NODE != -1) )
#else
#define TRACK_ROW(node, row, sys, link) false
#endif

#ifdef TRACK_R
#define I_TRACK_ROW                                \
   (!nodeIsDummy(ROW_NODE)        \
      && (my_rank == 0 || ROW_IS_LINK              \
         || ROW_NODE != -1))
#else
#define I_TRACK_ROW false
#endif

#ifdef TRACK_C
#define I_TRACK_COLUMN                             \
   !nodeIsDummy(COL_NODE)         \
      && (my_rank == 0 || COL_NODE != -1)
#else
#define I_TRACK_COLUMN false
#endif

PresolveData::PresolveData(const sData* sorigprob, StochPostsolver* postsolver) :
      postsolver(postsolver),
      length_array_outdated_indicators(5),
      array_outdated_indicators(new bool[length_array_outdated_indicators]),
      outdated_lhsrhs(array_outdated_indicators[0]),
      outdated_nnzs(array_outdated_indicators[1]),
      outdated_linking_var_bounds(array_outdated_indicators[2]),
      outdated_activities(array_outdated_indicators[3]),
      outdated_obj_vector(array_outdated_indicators[4]),
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
      nChildren(nnzs_col->children.size()),
      objOffset(0.0), obj_offset_chgs(0.0),
      objective_vec_chgs(dynamic_cast<SimpleVector*>(nnzs_col->vec->cloneFull())),
      lower_bound_implied_by_system(dynamic_cast<StochVectorBase<int>*>(nnzs_col->clone())),
      lower_bound_implied_by_row(dynamic_cast<StochVectorBase<int>*>(nnzs_col->clone())),
      lower_bound_implied_by_node(dynamic_cast<StochVectorBase<int>*>(nnzs_col->clone())),
      upper_bound_implied_by_system(dynamic_cast<StochVectorBase<int>*>(nnzs_col->clone())),
      upper_bound_implied_by_row(dynamic_cast<StochVectorBase<int>*>(nnzs_col->clone())),
      upper_bound_implied_by_node(dynamic_cast<StochVectorBase<int>*>(nnzs_col->clone())),
      elements_deleted(0), elements_deleted_transposed(0)
{
   std::memset(array_outdated_indicators, 0, length_array_outdated_indicators * sizeof(bool) );
   outdated_activities = true;

   lower_bound_implied_by_system->setToConstant(-10);
   lower_bound_implied_by_row->setToConstant(-10);
   lower_bound_implied_by_node->setToConstant(-10);
   upper_bound_implied_by_system->setToConstant(-10);
   upper_bound_implied_by_row->setToConstant(-10);
   upper_bound_implied_by_node->setToConstant(-10);

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
   dynamic_cast<StochGenMatrix&>(*presProb->A).initTransposed(true);
   dynamic_cast<StochGenMatrix&>(*presProb->C).initTransposed(true);

   recomputeActivities();
   initNnzCounter( *nnzs_row_A, *nnzs_row_C, *nnzs_col);

   initSingletons();
   setUndefinedVarboundsTo(std::numeric_limits<double>::max());

#ifdef TRACK_C
   if( I_TRACK_COLUMN )
   {
      std::cout << "TRACKING_COLUMN: colum " << COLUMN << " node " << COL_NODE << std::endl;
   }
#endif

#ifdef TRACK_R
   if( I_TRACK_ROW )
   {
      std::cout << "TRACKING_ROW: row " << ROW << " node " << ROW_NODE << " in BlockType " << ROW_BLOCK << " SystemType " << ROW_SYS << std::endl;
      writeRowLocalToStreamDense(std::cout, ROW_SYS, ROW_NODE, ROW_BLOCK, ROW);
   }
#endif
}

PresolveData::~PresolveData()
{
   delete[] array_outdated_indicators;
   delete[] array_nnz_chgs;
   delete[] array_act_chgs;
   delete[] array_bound_chgs;
   delete[] array_act_unbounded_chgs;
}

/* set non existent bounds on linking variables to +/- double max */
void PresolveData::setUndefinedVarboundsTo(double value)
{
   StochVector& xlow = dynamic_cast<StochVector&>(*presProb->blx);
   StochVector& ixlow = dynamic_cast<StochVector&>(*presProb->ixlow);
   StochVector& xupp = dynamic_cast<StochVector&>(*presProb->bux);
   StochVector& ixupp = dynamic_cast<StochVector&>(*presProb->ixupp);

   SimpleVector& vec_xlow = dynamic_cast<SimpleVector&>(*xlow.vec);
   SimpleVector& vec_ixlow = dynamic_cast<SimpleVector&>(*ixlow.vec);
   SimpleVector& vec_xupp = dynamic_cast<SimpleVector&>(*xupp.vec);
   SimpleVector& vec_ixupp = dynamic_cast<SimpleVector&>(*ixupp.vec);
   for(int i = 0; i < vec_xlow.length(); ++i)
   {
      if( PIPSisZero(vec_ixlow[i]) )
         vec_xlow[i] = -value;

      if( PIPSisZero(vec_ixupp[i]) )
         vec_xupp[i] = value;
   }
}

sData* PresolveData::finalize()
{

#ifndef NDEBUG
   if(distributed)
      PIPS_MPIlogicOrArrayInPlace(array_outdated_indicators, length_array_outdated_indicators, MPI_COMM_WORLD);
   assert(!outdated_activities && !outdated_lhsrhs && !outdated_nnzs && !outdated_linking_var_bounds && !outdated_obj_vector);
#endif

   /* theoretically it should not matter but there is an assert later which needs all these to be zero */
   setUndefinedVarboundsTo(0.0);

   // this removes all columns and rows that are now empty from the problem
   presProb->cleanUpPresolvedData(*nnzs_row_A, *nnzs_row_C, *nnzs_col);

   dynamic_cast<StochGenMatrix&>(*presProb->A).deleteTransposed();
   dynamic_cast<StochGenMatrix&>(*presProb->C).deleteTransposed();

   return presProb;
}

/** Recomputes the activities of all rows the process knows about. If linking_only is set to true only the linking_rows will get recomputed.
 *  Careful, recomputing linking rows requires MPI communication. Ideally all activities only have to be computed once, when creating the
 *  PresolveData.
 *  After that changes in the activities of linking rows will get stored in the SimpleVectors
 */
void PresolveData::recomputeActivities(bool linking_only)
{
   const StochGenMatrix& mat_A = dynamic_cast<const StochGenMatrix&>(*presProb->A);
   const StochGenMatrix& mat_C = dynamic_cast<const StochGenMatrix&>(*presProb->C);

   const StochVector& xupp = dynamic_cast<StochVector&>(*presProb->bux);
   const StochVector& ixupp = dynamic_cast<StochVector&>(*presProb->ixupp);
   const StochVector& xlow = dynamic_cast<StochVector&>(*presProb->blx);
   const StochVector& ixlow = dynamic_cast<StochVector&>(*presProb->ixlow);

   /* reset vectors keeping track of activities */
   if(!linking_only)
   {
      actmin_eq_part->setToZero();
      actmax_eq_part->setToZero();
      actmin_ineq_part->setToZero();
      actmax_ineq_part->setToZero();

      actmin_eq_ubndd->setToZero();
      actmax_eq_ubndd->setToZero();
      actmin_ineq_ubndd->setToZero();
      actmax_ineq_ubndd->setToZero();
   }
   else
   {
      actmin_eq_part->vecl->setToZero();
      actmax_eq_part->vecl->setToZero();
      actmin_ineq_part->vecl->setToZero();
      actmax_ineq_part->vecl->setToZero();

      actmin_eq_ubndd->vecl->setToZero();
      actmax_eq_ubndd->vecl->setToZero();
      actmin_ineq_ubndd->vecl->setToZero();
      actmax_ineq_ubndd->vecl->setToZero();
   }

   /* compute activities at root node */
   const SimpleVector& xupp_root = dynamic_cast<const SimpleVector&>(*xupp.vec);
   const SimpleVector& ixupp_root = dynamic_cast<const SimpleVector&>(*ixupp.vec);
   const SimpleVector& xlow_root = dynamic_cast<const SimpleVector&>(*xlow.vec);
   const SimpleVector& ixlow_root = dynamic_cast<const SimpleVector&>(*ixlow.vec);

   /* A0/B0 */
   if(!linking_only)
   {
      SimpleVector& actmin_eq_root_part = dynamic_cast<SimpleVector&>(*actmin_eq_part->vec);
      SimpleVector& actmax_eq_root_part = dynamic_cast<SimpleVector&>(*actmax_eq_part->vec);
      SimpleVector& actmin_ineq_root_part = dynamic_cast<SimpleVector&>(*actmin_ineq_part->vec);
      SimpleVector& actmax_ineq_root_part = dynamic_cast<SimpleVector&>(*actmax_ineq_part->vec);

      SimpleVectorBase<int>& actmin_eq_root_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmin_eq_ubndd->vec);
      SimpleVectorBase<int>& actmax_eq_root_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmax_eq_ubndd->vec);
      SimpleVectorBase<int>& actmin_ineq_root_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmin_ineq_ubndd->vec);
      SimpleVectorBase<int>& actmax_ineq_root_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmax_ineq_ubndd->vec);

      addActivityOfBlock(mat_A.Bmat->getStorageDynamicRef(), actmin_eq_root_part, actmin_eq_root_ubndd, actmax_eq_root_part, actmax_eq_root_ubndd,
            xlow_root, ixlow_root, xupp_root, ixupp_root);

      addActivityOfBlock(mat_C.Bmat->getStorageDynamicRef(), actmin_ineq_root_part, actmin_ineq_root_ubndd, actmax_ineq_root_part, actmax_ineq_root_ubndd,
            xlow_root, ixlow_root, xupp_root, ixupp_root);
   }

   SimpleVector& actmin_eq_link_part = dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl);
   SimpleVector& actmax_eq_link_part = dynamic_cast<SimpleVector&>(*actmax_eq_part->vecl);
   SimpleVector& actmin_ineq_link_part = dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl);
   SimpleVector& actmax_ineq_link_part = dynamic_cast<SimpleVector&>(*actmax_ineq_part->vecl);

   SimpleVectorBase<int>& actmin_eq_link_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmin_eq_ubndd->vecl);
   SimpleVectorBase<int>& actmax_eq_link_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmax_eq_ubndd->vecl);
   SimpleVectorBase<int>& actmin_ineq_link_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmin_ineq_ubndd->vecl);
   SimpleVectorBase<int>& actmax_ineq_link_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmax_ineq_ubndd->vecl);

   /* Bl0 */
   if(my_rank == 0)
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
            SimpleVector& actmin_eq_child_part = dynamic_cast<SimpleVector&>(*actmin_eq_part->children[node]->vec);
            SimpleVector& actmax_eq_child_part = dynamic_cast<SimpleVector&>(*actmax_eq_part->children[node]->vec);

            SimpleVectorBase<int>& actmin_eq_child_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmin_eq_ubndd->children[node]->vec);
            SimpleVectorBase<int>& actmax_eq_child_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmax_eq_ubndd->children[node]->vec);

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
            SimpleVector& actmin_ineq_child_part = dynamic_cast<SimpleVector&>(*actmin_ineq_part->children[node]->vec);
            SimpleVector& actmax_ineq_child_part = dynamic_cast<SimpleVector&>(*actmax_ineq_part->children[node]->vec);

            SimpleVectorBase<int>& actmin_ineq_child_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmin_ineq_ubndd->children[node]->vec);
            SimpleVectorBase<int>& actmax_ineq_child_ubndd = dynamic_cast<SimpleVectorBase<int>&>(*actmax_ineq_ubndd->children[node]->vec);

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
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVector*>(actmin_eq_part->vecl)->elements(), actmin_eq_part->vecl->n, MPI_COMM_WORLD);
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVector*>(actmax_eq_part->vecl)->elements(), actmax_eq_part->vecl->n, MPI_COMM_WORLD);
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVector*>(actmin_ineq_part->vecl)->elements(), actmin_ineq_part->vecl->n, MPI_COMM_WORLD);
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVector*>(actmax_ineq_part->vecl)->elements(), actmax_ineq_part->vecl->n, MPI_COMM_WORLD);

      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVectorBase<int>*>(actmin_eq_ubndd->vecl)->elements(), actmin_eq_ubndd->vecl->n, MPI_COMM_WORLD);
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVectorBase<int>*>(actmax_eq_ubndd->vecl)->elements(), actmax_eq_ubndd->vecl->n, MPI_COMM_WORLD);
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVectorBase<int>*>(actmin_ineq_ubndd->vecl)->elements(), actmin_ineq_ubndd->vecl->n, MPI_COMM_WORLD);
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVectorBase<int>*>(actmax_ineq_ubndd->vecl)->elements(), actmax_ineq_ubndd->vecl->n, MPI_COMM_WORLD);
   }

   /* set activities to infinity // theoretically not necessary but for debugging */
   for(int row = 0; row < actmin_eq_part->vecl->n; ++row)
   {
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmin_eq_ubndd->vecl)[row] >= 2 )
         dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl)[row] = -std::numeric_limits<double>::infinity();
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmax_eq_ubndd->vecl)[row] >= 2 )
         dynamic_cast<SimpleVector&>(*actmax_eq_part->vecl)[row] = std::numeric_limits<double>::infinity();
   }
   for(int row = 0; row < actmin_ineq_part->vecl->n; ++row)
   {
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmin_ineq_ubndd->vecl)[row] >= 2 )
         dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl)[row] = -std::numeric_limits<double>::infinity();
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmax_ineq_ubndd->vecl)[row] >= 2 )
         dynamic_cast<SimpleVector&>(*actmax_ineq_part->vecl)[row] = std::numeric_limits<double>::infinity();
   }

   actmax_eq_chgs->setToZero();
   actmin_eq_chgs->setToZero();
   actmax_ineq_chgs->setToZero();
   actmin_ineq_chgs->setToZero();

   actmax_eq_ubndd_chgs->setToZero();
   actmin_eq_ubndd_chgs->setToZero();
   actmax_ineq_ubndd_chgs->setToZero();
   actmin_ineq_ubndd_chgs->setToZero();

   outdated_activities = false;

#ifdef TRACK_R
   if( I_TRACK_ROW )
   {
   double act_min = (ROW_SYS == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmin_eq_part, ROW_NODE, ROW_IS_LINK )[ROW] :
         getSimpleVecFromRowStochVec(*actmin_ineq_part, ROW_NODE, ROW_IS_LINK )[ROW];
   double act_max = (ROW_SYS == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmax_eq_part, ROW_NODE, ROW_IS_LINK )[ROW] :
         getSimpleVecFromRowStochVec(*actmax_ineq_part, ROW_NODE, ROW_IS_LINK )[ROW];

   int act_min_ubndd = (ROW_SYS == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmin_eq_ubndd, ROW_NODE, ROW_IS_LINK )[ROW] :
         getSimpleVecFromRowStochVec(*actmin_ineq_ubndd, ROW_NODE, ROW_IS_LINK )[ROW];
   int act_max_ubndd = (ROW_SYS == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmax_eq_ubndd, ROW_NODE, ROW_IS_LINK )[ROW] :
         getSimpleVecFromRowStochVec(*actmax_ineq_ubndd, ROW_NODE, ROW_IS_LINK )[ROW];


   std::cout << "TRACKING_ROW: Recomputed activity of tracked row" << std::endl;
   std::cout << "\tnew min/max activity is: " << act_min << "/" << act_max << ", min/max unbounded counters are " << act_min_ubndd << "/" << act_max_ubndd << std::endl;
   }
#endif

}


/** Computes minimal and maximal activity of all rows in given matrix. Adds activities to min/max_activities accordingly. */
void PresolveData::addActivityOfBlock( const SparseStorageDynamic& matrix, SimpleVector& min_partact, SimpleVectorBase<int>& unbounded_min,
   SimpleVector& max_partact, SimpleVectorBase<int>& unbounded_max, const SimpleVector& xlow, const SimpleVector& ixlow,
   const SimpleVector& xupp, const SimpleVector& ixupp) const
{
   assert( xlow.n == matrix.n && ixlow.n == matrix.n && xupp.n == matrix.n && ixupp.n == matrix.n );
   assert( max_partact.n == matrix.m && min_partact.n == matrix.m);
   assert( unbounded_min.n == matrix.m && unbounded_max.n == matrix.m);

   for( int row = 0; row < matrix.m; ++row)
   {
      for( int j = matrix.rowptr[row].start; j < matrix.rowptr[row].end; j++)
      {
         const int col = matrix.jcolM[j];
         const double entry = matrix.M[j];

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
   PIPS_MPIgetLogicOrInPlace( outdated_linking_var_bounds, MPI_COMM_WORLD );

   if(!outdated_linking_var_bounds)
      return;

   if(distributed)
   {
      SimpleVector& xlow = getSimpleVecFromColStochVec(*presProb->blx, -1);
      SimpleVector& xupp = getSimpleVecFromColStochVec(*presProb->bux, -1);
      SimpleVector& ixlow = getSimpleVecFromColStochVec(*presProb->ixlow, -1);
      SimpleVector& ixupp = getSimpleVecFromColStochVec(*presProb->ixupp, -1);

      /* copy old values for later compairson */
      SimpleVector* xlow_old = dynamic_cast<SimpleVector*>(xlow.cloneFull());
      SimpleVector* xupp_old = dynamic_cast<SimpleVector*>(xupp.cloneFull());
      SimpleVector* ixlow_old = dynamic_cast<SimpleVector*>(ixlow.cloneFull());
      SimpleVector* ixupp_old = dynamic_cast<SimpleVector*>(ixupp.cloneFull());

      PIPS_MPImaxArrayInPlace(xlow.elements(), xlow.length(), MPI_COMM_WORLD);
      PIPS_MPImaxArrayInPlace(ixlow.elements(), ixlow.length(), MPI_COMM_WORLD);
      PIPS_MPImaxArrayInPlace(xupp.elements(), xupp.length(), MPI_COMM_WORLD);
      PIPS_MPImaxArrayInPlace(ixupp.elements(), ixupp.length(), MPI_COMM_WORLD);

      // this will affect the activities of basically all rows - use with care
      for(int col = 0; col < xlow.length(); ++col)
      {
         double xu_old = std::numeric_limits<double>::infinity();
         double xl_old = -std::numeric_limits<double>::infinity();
         double xu = std::numeric_limits<double>::infinity();
         double xl = -std::numeric_limits<double>::infinity();

         if( PIPSisEQ( (*ixupp_old)[col], 1.0) )
            xu_old = (*xupp_old)[col];
         if( PIPSisEQ( (*ixlow_old)[col], 1.0) )
            xl_old = (*xlow_old)[col];
         if( PIPSisEQ( ixupp[col], 1.0) )
            xu = xupp[col];
         if( PIPSisEQ( ixlow[col], 1.0) )
            xl = xlow[col];

         updateRowActivities(-1, col, xu, xl, xu_old, xl_old);
      }
      
      delete xlow_old;
      delete xupp_old;
      delete ixlow_old;
      delete ixupp_old;
      
      }

   outdated_linking_var_bounds = false;
}

/** allreduce changes in the activities of the linking rows and update the linking row activities */
void PresolveData::allreduceAndApplyLinkingRowActivities()
{
   /* check for new linking rows that have to be computed */
   PIPS_MPIgetLogicOrInPlace(outdated_activities, MPI_COMM_WORLD);
   // todo : better criterion - if there have been this and that many changes etc
   PIPS_MPIgetSumInPlace(linking_rows_need_act_computation, MPI_COMM_WORLD);

   if(!outdated_activities && linking_rows_need_act_computation == 0)
      return;

   /* allreduce unbounded changes, compute rows that now have at most one unbounded entry, strore the local
    * rows in the changes array, allredues MPI_SUM the changes array and update local activities with the global ones
    */
   if(distributed)
   {
      PIPS_MPIsumArrayInPlace(array_act_unbounded_chgs, lenght_array_act_chgs, MPI_COMM_WORLD);
   }

   dynamic_cast<SimpleVectorBase<int>&>(*actmin_eq_ubndd->vecl).axpy( 1.0, *actmin_eq_ubndd_chgs);
   dynamic_cast<SimpleVectorBase<int>&>(*actmax_eq_ubndd->vecl).axpy( 1.0, *actmax_eq_ubndd_chgs);
   dynamic_cast<SimpleVectorBase<int>&>(*actmin_ineq_ubndd->vecl).axpy( 1.0, *actmin_ineq_ubndd_chgs);
   dynamic_cast<SimpleVectorBase<int>&>(*actmax_ineq_ubndd->vecl).axpy( 1.0, *actmax_ineq_ubndd_chgs);

   /* equality system */
   for(int row = 0; row < actmin_eq_ubndd->vecl->n; ++row)
   {
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmin_eq_ubndd->vecl)[row] < 2
            && dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl)[row] == -std::numeric_limits<double>::infinity())
      {
         (*actmin_eq_chgs)[row] = computeLocalLinkingRowMinOrMaxActivity(EQUALITY_SYSTEM, row, false);
         dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl)[row] = 0;
      }

      if( dynamic_cast<SimpleVectorBase<int>&>(*actmax_eq_ubndd->vecl)[row] < 2
            && dynamic_cast<SimpleVector&>(*actmax_eq_part->vecl)[row] == std::numeric_limits<double>::infinity())
      {
         (*actmax_eq_chgs)[row] = computeLocalLinkingRowMinOrMaxActivity(EQUALITY_SYSTEM, row, true);
         dynamic_cast<SimpleVector&>(*actmax_eq_part->vecl)[row] = 0;
      }
   }

   /* inequality system */
   for(int row = 0; row < actmin_ineq_ubndd->vecl->n; ++row)
   {
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmin_ineq_ubndd->vecl)[row] < 2
            && dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl)[row] == -std::numeric_limits<double>::infinity())
      {
         (*actmin_ineq_chgs)[row] = computeLocalLinkingRowMinOrMaxActivity(INEQUALITY_SYSTEM, row, false);
         dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl)[row] = 0;
      }

      if( dynamic_cast<SimpleVectorBase<int>&>(*actmax_ineq_ubndd->vecl)[row] < 2
            && dynamic_cast<SimpleVector&>(*actmax_ineq_part->vecl)[row] == std::numeric_limits<double>::infinity())
      {
         (*actmax_ineq_chgs)[row] = computeLocalLinkingRowMinOrMaxActivity(INEQUALITY_SYSTEM, row, true);
         dynamic_cast<SimpleVector&>(*actmax_ineq_part->vecl)[row] = 0;
      }
   }

   if( distributed )
      PIPS_MPIsumArrayInPlace(array_act_chgs, lenght_array_act_chgs, MPI_COMM_WORLD);

   dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl).axpy( 1.0, *actmin_eq_chgs);
   dynamic_cast<SimpleVector&>(*actmax_eq_part->vecl).axpy( 1.0, *actmax_eq_chgs);
   dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl).axpy( 1.0, *actmin_ineq_chgs);
   dynamic_cast<SimpleVector&>(*actmax_ineq_part->vecl).axpy( 1.0, *actmax_ineq_chgs);

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
         assert( std::fabs(dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl)[i]) != std::numeric_limits<double>::infinity() );
      else
         assert( dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl)[i] == -std::numeric_limits<double>::infinity() );
   }

   for(int i = 0; i < actmin_ineq_ubndd->vecl->n; ++i)
   {
      if( dynamic_cast<SimpleVectorBase<int>&>(*actmin_ineq_ubndd->vecl)[i] < 2 )
         assert( std::fabs(dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl)[i]) != std::numeric_limits<double>::infinity() );
      else
         assert( dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl)[i] == -std::numeric_limits<double>::infinity() );
   }
#endif
}

/** allreduce changes in the nnz counters and apply them locally */
void PresolveData::allreduceAndApplyNnzChanges()
{
   PIPS_MPIgetLogicOrInPlace(outdated_nnzs, MPI_COMM_WORLD);

   if(!outdated_nnzs)
      return;

   if( distributed )
      PIPS_MPIsumArrayInPlace(array_nnz_chgs, length_array_nnz_chgs, MPI_COMM_WORLD);

   /* update local nnzCounters */
   nnzs_col->vec->axpy(1, *nnzs_col_chgs);
   nnzs_row_A->vecl->axpy(1, *nnzs_row_A_chgs);
   nnzs_row_C->vecl->axpy(1, *nnzs_row_C_chgs);

   // todo : this can be done more efficiently, e.g. while substracting
   // todo : this has still flaws - new singleton rows in B0 and D0 are not communicated properly for some reason - might happen else where
   for( int i = 0; i < nnzs_col_chgs->length(); ++i )
   {
      if( (*nnzs_col_chgs)[i] > 0 && dynamic_cast<SimpleVectorBase<int>&>(*nnzs_col->vec)[i] == 1 )
         singleton_cols.push( sCOLINDEX(-1, i) );
   }
   
   for( int i = 0; i < nnzs_row_A_chgs->length(); ++i )
   {
      if( (*nnzs_row_A_chgs)[i] > 0 && dynamic_cast<SimpleVectorBase<int>&>(*nnzs_row_A->vec)[i] == 1 )
         singleton_rows.push( sROWINDEX(EQUALITY_SYSTEM, -2, i) );
   }
   
   for( int i = 0; i < nnzs_row_C_chgs->length(); ++i )
   {
      if( (*nnzs_row_C_chgs)[i] > 0 && dynamic_cast<SimpleVectorBase<int>&>(*nnzs_row_C->vec)[i] == 1 )
         singleton_rows.push( sROWINDEX(INEQUALITY_SYSTEM, -2, i) );
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
   PIPS_MPIgetLogicOrInPlace(outdated_lhsrhs, MPI_COMM_WORLD);

   if(!outdated_lhsrhs)
      return;

   if(distributed)
   {
      PIPS_MPIsumArrayInPlace(array_bound_chgs, lenght_array_bound_chgs, MPI_COMM_WORLD);
   }

   dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*presProb->bA).vecl)->axpy( 1.0, *bound_chgs_A);
   dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*presProb->bl).vecl)->axpy( 1.0, *bound_chgs_C);
   dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*presProb->bu).vecl)->axpy( 1.0, *bound_chgs_C);

   bound_chgs_A->setToZero();
   bound_chgs_C->setToZero();

   outdated_lhsrhs = false;
}

void PresolveData::allreduceAndApplyObjVecChanges()
{
   PIPS_MPIgetLogicOrInPlace(outdated_obj_vector, MPI_COMM_WORLD);

   if(!outdated_obj_vector)
      return;

   if(distributed)
      PIPS_MPIsumArrayInPlace(objective_vec_chgs->elements(), objective_vec_chgs->length(), MPI_COMM_WORLD);

   dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*presProb->g).vec)->axpy(1.0, *objective_vec_chgs);

   objective_vec_chgs->setToZero();
   outdated_obj_vector = false;
}

void PresolveData::allreduceObjOffset()
{
   if(distributed)
      PIPS_MPIgetSumInPlace(obj_offset_chgs, MPI_COMM_WORLD);

   objOffset += obj_offset_chgs;
   obj_offset_chgs = 0;
}

void PresolveData::initNnzCounter(StochVectorBase<int>& nnzs_row_A, StochVectorBase<int>& nnzs_row_C, StochVectorBase<int>& nnzs_col) const
{
   StochGenMatrix& A = dynamic_cast<StochGenMatrix&>(*(presProb->A));
   StochGenMatrix& C = dynamic_cast<StochGenMatrix&>(*(presProb->C));

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
         singleton_rows.push( sROWINDEX( EQUALITY_SYSTEM, -1, i ));
   }

   /* Bl0 */
   if(nnzs_row_A->vecl != NULL)
   {
      for(int i = 0; i < nnzs_row_A->vecl->n; ++i)
      {
         if( dynamic_cast<SimpleVectorBase<int>&>(*nnzs_row_A->vecl)[i] == 1)
            singleton_rows.push( sROWINDEX( EQUALITY_SYSTEM, -2, i ));
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
               singleton_rows.push( sROWINDEX( EQUALITY_SYSTEM, i, j ));
         }
      }
   }

   /* rows of C */
   /* B0 */
   for(int i = 0; i < nnzs_row_C->vec->n; ++i)
      if( dynamic_cast<SimpleVectorBase<int>&>(*nnzs_row_C->vec)[i] == 1)
         singleton_rows.push( sROWINDEX( INEQUALITY_SYSTEM, -1, i ));

   /* Bl0 */
   if( nnzs_row_C->vecl != NULL)
   {
      for( int i = 0; i < nnzs_row_C->vecl->n; ++i )
      {
         if( dynamic_cast<SimpleVectorBase<int>&>(*nnzs_row_C->vecl)[i] == 1 )
            singleton_rows.push(sROWINDEX( INEQUALITY_SYSTEM, -2, i ));
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
               singleton_rows.push(sROWINDEX( INEQUALITY_SYSTEM, i, j ));
         }
      }
   }

   /* columns */
   for(int i = 0; i < nnzs_col->vec->n; ++i)
   {
      if( dynamic_cast<SimpleVectorBase<int>&>(*nnzs_col->vec)[i] == 1)
         singleton_cols.push( sCOLINDEX(-1, i));
   }

   for(int i = 0; i < nChildren; ++i)
   {
      if(nnzs_col->children[i]->vec != NULL)
      {
         for(int j = 0; j < nnzs_col->children[i]->vec->n; ++j)
         {
            if( dynamic_cast<SimpleVectorBase<int>&>(*nnzs_col->children[i]->vec)[j] == 1)
               singleton_cols.push( sCOLINDEX(i, j));
         }
      }
   }

}

bool PresolveData::reductionsEmpty()
{
   double recv = 0;
   if(distributed)
   {
      PIPS_MPIlogicOrArrayInPlace(array_outdated_indicators, length_array_outdated_indicators, MPI_COMM_WORLD);

      recv = PIPS_MPIgetSum(obj_offset_chgs, MPI_COMM_WORLD);
   }
   return !outdated_obj_vector && !outdated_activities && !outdated_lhsrhs && !outdated_linking_var_bounds && !outdated_nnzs && (recv == 0);
}

// todo : if small entry was removed from system no postsolve is necessary - if coefficient was removed because impact of changes in variable are small
// also rhs lhs will be adjusted - this has to be reversed later - there will also be a problem with the reduced costs in than particular row ?
// todo : postsolve if bounds adjusted because of deleted matrix entry simply reverse the adjustment - no changes in multipliers
//- row stays active / inactive ? not true..
void PresolveData::deleteEntry(SystemType system_type, int node, BlockType block_type, int row_index,
      int& index_k, int& row_end)
{
   assert(-1 <= node && node < nChildren);
   assert(node != -1 || block_type != A_MAT);
   const bool linking = (block_type == BL_MAT);

   SparseStorageDynamic* storage = getSparseGenMatrix(system_type, node , block_type)->getStorageDynamic();
   assert(0 <= row_index && row_index <= storage->m);

   SimpleVector& xlower = getSimpleVecFromColStochVec(*presProb->blx, node);

   if(block_type == BL_MAT || block_type == A_MAT || node == -1)
      outdated_nnzs = true;

   /* adjust rhs and lhs */
   adjustMatrixRhsLhsBy(system_type, node, linking, row_index, -storage->M[index_k] * xlower[storage->jcolM[index_k]]);

   /* adjust activity */
   outdated_activities = true;

   /* adjust nnz counters */
   reduceNnzCounterRow(system_type, node, linking, row_index, 1);
   reduceNnzCounterColumn(node, block_type, storage->jcolM[index_k], 1);

   std::swap(storage->M[index_k], storage->M[row_end - 1]);
   std::swap(storage->jcolM[index_k], storage->jcolM[row_end - 1]);
   storage->rowptr[row_index].end--;
   row_end = storage->rowptr[row_index].end;
   index_k--;

   ++elements_deleted; // todo
}

void PresolveData::resetOriginallyFreeVarsBounds(const sData& orig_prob)
{

#ifndef NDEBUG
   if(my_rank == 0)
      std::cout << "Resetting all presolved variable bounds of originally free variables" <<::endl; 
#endif
   long long n = 0;
   for( int node = -1; node < nChildren; ++node )
   {
      n += resetOriginallyFreeVarsBounds( getSimpleVecFromColStochVec(*orig_prob.ixlow, node), getSimpleVecFromColStochVec(*orig_prob.ixupp, node), node);
      if(my_rank != 0 && node == -1)
      {
         n = 0;
      }
   }

#ifndef NDEBUG
   PIPS_MPIgetSumInPlace(n, MPI_COMM_WORLD);
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

void PresolveData::fixEmptyColumn(int node, int col, double val)
{
   assert(-1 <= node && node < nChildren);
   assert(0 <= col);

   postsolver->notifyFixedEmptyColumn(node, col, val);

   removeColumn(node, col, val);

   if( node != -1)
      assert( getSimpleVecFromColStochVec(*nnzs_col, node)[col] == 0 );
}

void PresolveData::fixColumn(int node, int col, double value)
{
   if( TRACK_COLUMN(node, col) )
   {
      std::cout << "TRACKING_COLUMN: column " << col << " node " << node << " got fixed to " << value << std::endl;
   }

   assert( -1 <= node && node < nChildren );
   assert( 0 <= col );
   
   /* current upper and lower bound as well als column - if linking variable then only proc zero stores current root column */
   std::vector<int> idx;
   std::vector<double> val;
   
   // todo
   postsolver->notifyFixedColumn(node, col, value, idx, val);

#ifndef NDEBUG
   double ixlow = getSimpleVecFromColStochVec(*presProb->ixlow, node)[col];
   double ixupp = getSimpleVecFromColStochVec(*presProb->ixupp, node)[col];
   double xlow = PIPSisEQ(ixupp, 1.0) ? getSimpleVecFromColStochVec(*presProb->blx, node)[col] : std::numeric_limits<double>::infinity();
   double xupp = PIPSisEQ(ixlow, 1.0) ? getSimpleVecFromColStochVec(*presProb->bux, node)[col] : std::numeric_limits<double>::infinity();
   assert( PIPSisEQ(ixlow, 1.0) );
   assert( PIPSisEQ(ixupp, 1.0) );
   assert( PIPSisEQ(xlow, xupp, 1e-10) );
   assert( PIPSisEQ(xlow, value, 1e-10) );
#endif

   removeColumn(node, col, value);

   if( node != -1)
      assert( getSimpleVecFromColStochVec(*nnzs_col, node)[col] == 0 );
}

bool PresolveData::rowPropagatedBounds( SystemType system_type, int node_row, BlockType block_type, int row, int col, double ubx, double lbx)
{
   assert( -1 <= node_row && node_row < nChildren );
   
   const int node_var = (block_type == A_MAT) ? -1 : node_row;
   assert( 0 <= col && col < getSimpleVecFromColStochVec( *presProb->ixlow, node_var ).n );

   // todo : choose numerical threshold
   //const double numerical_upper_threshold = std::numeric_limits<double>::max();

   /* check for infeasibility of the newly found bounds */
   const double ixlow = getSimpleVecFromColStochVec( *presProb->ixlow, node_var )[col];
   const double xlow = getSimpleVecFromColStochVec( *presProb->blx, node_var )[col];
   const double ixupp = getSimpleVecFromColStochVec( *presProb->ixupp, node_var )[col];
   const double xupp = getSimpleVecFromColStochVec( *presProb->bux, node_var )[col];
   
   if( TRACK_COLUMN(node_var,col) )
   {
      std::cout << "TRACKING_COLUMN: new bounds [" << lbx << ", " << ubx << "] propagated for column " << 
         col << " from row " << row << " node " << node_row << " in " <<
         ( (system_type == EQUALITY_SYSTEM) ? "EQU_SYS" : "INEQ_SYS") << ":" << std::endl;
      std::cout << "\tbounds were [" << xlow << ", " << xupp << "]" << std::endl;
   }

   if( ( !PIPSisZero(ixlow) && PIPSisLT(ubx, xlow) )
         || ( !PIPSisZero(ixupp) && PIPSisLT(xupp, lbx) )
         || PIPSisLT(ubx, lbx) )
   {
      std::cout << "[" << lbx << ", " << ubx << "] not in [" << ( !PIPSisZero(ixlow) ? xlow : -std::numeric_limits<double>::infinity()) << ", " << 
         ( !PIPSisZero(ixupp) ? xupp : std::numeric_limits<double>::infinity()) << "]" << std::endl;
      PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Row Propagation detected infeasible new bounds!", "PresolveData.C", "rowPropagatedBounds");
   }

   /* adjust bounds */
   bool bounds_changed = false;

   // we do not tighten bounds if impact is too low or bound is bigger than 10e8 // todo : maybe different limit
   // set lower bound
   // if( fabs(lbx) < 1e8 && (PIPSisZERO(ixlow) || feastol * 1e3 <= fabs(xlow - lbx) ) )
   if( ubx < std::numeric_limits<double>::infinity() && ( PIPSisZero(ixupp) || PIPSisLT(ubx, xupp) ) )
   {
      if( updateUpperBoundVariable(node_var, col, ubx) )
      {
         /* store node and row that implied the bound (necessary for resetting bounds later on) */
         getSimpleVecFromColStochVec(*upper_bound_implied_by_system, node_var)[col] = system_type;
         getSimpleVecFromColStochVec(*upper_bound_implied_by_row, node_var)[col] = row;
         getSimpleVecFromColStochVec(*upper_bound_implied_by_node, node_var)[col] = (block_type == BL_MAT) ? -2 : node_row; // -2 for linking rows

         bounds_changed = true;
      }
   }
  // if( fabs(ubx) < 1e8 && (PIPSisZero(ixupp) || feastol * 1e3 <= fabs(xupp- ubx) ) )
   if( -std::numeric_limits<double>::infinity() < lbx && ( PIPSisZero(ixlow) || PIPSisLT(xlow, lbx)) )
   {
      if( updateLowerBoundVariable(node_var, col, lbx) )
      {
         /* store node and row that implied the bound (necessary for resetting bounds later on) */
         getSimpleVecFromColStochVec(*lower_bound_implied_by_system, node_var)[col] = system_type;
         getSimpleVecFromColStochVec(*lower_bound_implied_by_row, node_var)[col] = row;
         getSimpleVecFromColStochVec(*lower_bound_implied_by_node, node_var)[col] = (block_type == BL_MAT) ? -2 : node_row; // -2 for linking rows

         bounds_changed = true;
      }
   }

   /// every process should have the same root node data thus all of them should propagate it's rows similarly
   if( bounds_changed && (node_row == -1 || block_type == A_MAT) )
      assert(outdated_linking_var_bounds == true);

// todo : how to undo propagations from linking constraint rows..
// todo : in case a linking row propagated we'll have to store the whole linking row
//   SparseGenMatrix* mat = getSparseGenMatrix(system_type, node_row, block_type);
//   assert(row < mat->getStorageDynamic()->m );
//
//   int row_start = mat->getStorageDynamic()->rowptr[row].start;
//   int row_end = mat->getStorageDynamic()->rowptr[row].end;
//
//   assert(row_start < row_end);
// if bounds_changed
//   postsolver->notifyRowPropagated(system_type, node_row, row, (block_type == BL_MAT), col, lbx, ubx, mat->getStorageDynamic()->M + row_start,
//         mat->getStorageDynamic()->jcolM + row_start, row_end - row_start );

   return bounds_changed;
}

void PresolveData::tightenRowBoundsParallelRow(SystemType system_type, int node, int row, double lhs, double rhs, bool linking)
{
   if( TRACK_ROW(node, row, system_type, linking) )
   {
      std::cout << "TRACKING_ROW: before RHS LHS adjustment" << std::endl;
      writeRowLocalToStreamDense(std::cout, system_type, node, linking, row);
   }

   assert(-1 <= node && node < nChildren);
   assert(!linking);
   assert(system_type == INEQUALITY_SYSTEM);


   if( linking )
   {
      assert(false);
      //outdated_lhsrhs = true;
      //(system_type == EQUALITY_SYSTEM) ? (*bound_chgs_A)[row] += value : (*bound_chgs_C)[row] += value;
      return;
   }

   if(system_type == EQUALITY_SYSTEM)
   {
      assert(false);
      //getSimpleVecFromRowStochVec(*presProb->bA, node, linking)[row] += value;
   }
   else
   {
      if( lhs != -std::numeric_limits<double>::infinity() )
      {
         if( !PIPSisEQ(getSimpleVecFromRowStochVec(*presProb->iclow, node, linking)[row], 1.0) )
         {
            getSimpleVecFromRowStochVec(*presProb->iclow, node, linking)[row] = 1.0;
            getSimpleVecFromRowStochVec(*presProb->bl, node, linking)[row] = lhs;

         }
         else
            getSimpleVecFromRowStochVec(*presProb->bl, node, linking)[row] = 
               std::max(getSimpleVecFromRowStochVec(*presProb->bl, node, linking)[row], lhs);
      }
      if( rhs != std::numeric_limits<double>::infinity() )
      {
         if( !PIPSisEQ(getSimpleVecFromRowStochVec(*presProb->icupp, node, linking)[row], 1.0) )
         {
            getSimpleVecFromRowStochVec(*presProb->icupp, node, linking)[row] = 1.0;
         }
         else
         {
            getSimpleVecFromRowStochVec(*presProb->bu, node, linking)[row] = 
               std::min(getSimpleVecFromRowStochVec(*presProb->bu, node, linking)[row], rhs);
         }
      }
      // todo!
   }

   if( TRACK_ROW(node, row, system_type, linking) )
   {
      std::cout << "TRACKING_ROW: after RHS LHS adjustment " << std::endl;
      writeRowLocalToStreamDense(std::cout, system_type, node, linking, row);
   }
}

/** this methods does not call any postsolve procedures but simply changes the bounds (lhs, rhs) of either A or B by value */
void PresolveData::adjustMatrixRhsLhsBy(SystemType system_type, int node, bool linking, int row, double value)
{
   if( TRACK_ROW(node, row, system_type, linking) )
   {
      std::cout << "TRACKING_ROW: RHS LHS of row " << row << " being adjusted by " << value << std::endl;
      writeRowLocalToStreamDense(std::cout, system_type, node, linking, row);
   }

   assert( -1 <= node && node < nChildren);
   if( PIPSisEQ(value, 0.0) )
      return;


   if(linking)
   {
      outdated_lhsrhs = true;
      (system_type == EQUALITY_SYSTEM) ? (*bound_chgs_A)[row] += value : (*bound_chgs_C)[row] += value;
      return;
   }

   if(system_type == EQUALITY_SYSTEM)
   {
      getSimpleVecFromRowStochVec(*presProb->bA, node, linking )[row] += value;
   }
   else
   {
      if( PIPSisEQ(getSimpleVecFromRowStochVec(*presProb->icupp, node, linking )[row], 1.0) )
         getSimpleVecFromRowStochVec(*presProb->bu, node, linking )[row] += value;

      if( PIPSisEQ(getSimpleVecFromRowStochVec(*presProb->iclow, node, linking )[row], 1.0 ) )
         getSimpleVecFromRowStochVec(*presProb->bl, node, linking )[row] += value;
   }

   if( TRACK_ROW(node, row, system_type, linking) )
   {
      std::cout << "TRACKING_ROW: after RHS LHS adjustment " << std::endl;
      writeRowLocalToStreamDense(std::cout, system_type, node, linking, row);
   }
}


// todo : make a finish block_deletion ? // todo - remove this and always update both matrices..
void PresolveData::updateTransposedSubmatrix( SystemType system_type, int node, BlockType block_type, std::vector<std::pair<int, int> >& elements)
{
   assert(-1 <= node && node < nChildren);
   assert(node != -1 || block_type != A_MAT);
   SparseStorageDynamic* transposed = getSparseGenMatrix(system_type, node, block_type)->getStorageDynamicTransposed();

   for( size_t i = 0; i < elements.size(); ++i )
   {
      std::pair<int, int> entry = elements.at(i);
      const int row_A = entry.first;
      const int row_At = entry.second;

      const int start = transposed->rowptr[row_At].start;
      const int end = transposed->rowptr[row_At].end;
      int col_At;

      for( col_At = start; col_At < end; col_At++ )
      {
         if( transposed->jcolM[col_At] == row_A )
            break;
      }

      std::swap(transposed->M[col_At], transposed->M[end - 1]);
      std::swap(transposed->jcolM[col_At], transposed->jcolM[end - 1]);
      transposed->rowptr[row_At].end--;

      ++elements_deleted_transposed;
   }

   assert(elements_deleted == elements_deleted_transposed);
   elements_deleted = 0;
   elements_deleted_transposed = 0;
}

void PresolveData::reduceNnzCounterRow(SystemType system_type, int node, bool linking, int row_index, int amount)
{
   assert(-1 <= node && node < nChildren);
   assert(0 <= amount);

   if(amount == 0)
      return;

   /* linking constraints get stored */
   if( linking )
   {
      if(my_rank == 0 || node != -1)
      {
         (system_type == EQUALITY_SYSTEM) ? ((*nnzs_row_A_chgs)[row_index] -= amount) : ((*nnzs_row_C_chgs)[row_index] -= amount);
         outdated_nnzs = true;
      }
   }
   else
   {
      if(system_type == EQUALITY_SYSTEM)
      {
         getSimpleVecFromRowStochVec(*nnzs_row_A, node, linking )[row_index] -= amount;
         if( getSimpleVecFromRowStochVec(*nnzs_row_A, node, linking )[row_index] == 1)
            singleton_rows.push( sROWINDEX( system_type, node, row_index ) );
         assert( 0 <= getSimpleVecFromRowStochVec(*nnzs_row_A, node, linking )[row_index] );
      }
      else
      {
         getSimpleVecFromRowStochVec(*nnzs_row_C, node, linking )[row_index] -= amount;
         if( getSimpleVecFromRowStochVec(*nnzs_row_C, node, linking )[row_index] == 1)
            singleton_rows.push( sROWINDEX( system_type, node, row_index ) );
         assert( 0 <= getSimpleVecFromRowStochVec(*nnzs_row_C, node, linking )[row_index] );
      }
   }
}

void PresolveData::reduceNnzCounterColumn(int node, BlockType block_type, int col_index, int amount)
{
   assert(-1 <= node && node < nChildren);
   if(amount == 0)
      return;

   /* linking constraints get stored */
   if(node == -1 || block_type == A_MAT)
   {
      if(my_rank == 0 || node != -1)
      {
         (*nnzs_col_chgs)[col_index] -= amount;
         outdated_nnzs = true;
      }
   }
   else
   {
      getSimpleVecFromColStochVec( *nnzs_col, node )[col_index] -= amount;
      if( getSimpleVecFromColStochVec( *nnzs_col, node )[col_index] == 1)
         singleton_cols.push( sCOLINDEX( node, col_index ) );
      assert(0 <= getSimpleVecFromColStochVec( *nnzs_col, node )[col_index] );
   }
}

/** removes a column from the whole system A, C by fixing x to fixation
 * updates non-zero counters, rhs, lhs, objective offset and activities
 * does not call postsolve routines
 */
void PresolveData::removeColumn(int node, int col, double fixation)
{
   assert( -1 <= node && node < nChildren );

   if(node == -1)
   {
      removeColumnFromMatrix(EQUALITY_SYSTEM, -1, B_MAT, col, fixation);
      removeColumnFromMatrix(INEQUALITY_SYSTEM, -1, B_MAT, col, fixation);

      if(hasLinking(EQUALITY_SYSTEM))
         removeColumnFromMatrix(EQUALITY_SYSTEM, -1, BL_MAT, col, fixation);
      if( hasLinking(INEQUALITY_SYSTEM) )
         removeColumnFromMatrix(INEQUALITY_SYSTEM, -1, BL_MAT, col, fixation);

      for(int i = 0; i < nChildren; ++i)
      {
         if(!nodeIsDummy(i))
         {
            removeColumnFromMatrix(EQUALITY_SYSTEM, i, A_MAT, col, fixation);
            removeColumnFromMatrix(INEQUALITY_SYSTEM, i, A_MAT, col, fixation);
         }
      }
   }
   else
   {
      removeColumnFromMatrix(EQUALITY_SYSTEM, node, B_MAT, col, fixation);
      removeColumnFromMatrix(INEQUALITY_SYSTEM, node, B_MAT, col, fixation);

      if(hasLinking(EQUALITY_SYSTEM))
         removeColumnFromMatrix(EQUALITY_SYSTEM, node, BL_MAT, col, fixation);
      if(hasLinking(INEQUALITY_SYSTEM))
         removeColumnFromMatrix(INEQUALITY_SYSTEM, node, BL_MAT, col, fixation);
   }

   /* adjust objective function */
   if(node != -1 || my_rank == 0)
   {
      double objective_factor = getSimpleVecFromColStochVec(*presProb->g, node)[col];
      obj_offset_chgs += objective_factor * fixation;

   }

   /* mark column as removed */
   getSimpleVecFromColStochVec(*presProb->g, node)[col] = 0.0;
   getSimpleVecFromColStochVec(*presProb->ixlow, node)[col] = 0.0;
   getSimpleVecFromColStochVec(*presProb->ixupp, node)[col] = 0.0;
   getSimpleVecFromColStochVec(*presProb->blx, node)[col] = 0.0;
   getSimpleVecFromColStochVec(*presProb->bux, node)[col] = 0.0;
}

/** remove column - adjust lhs, rhs and activity as well as nnz_counters */
void PresolveData::removeColumnFromMatrix(SystemType system_type, int node, BlockType block_type, int col, double fixation)
{
   assert(-1 <= node && node < nChildren);
   assert(node != -1 || block_type != A_MAT);
   const bool linking = (block_type == BL_MAT);
   SparseGenMatrix* mat = getSparseGenMatrix(system_type, node, block_type);

   SparseStorageDynamic& matrix = mat->getStorageDynamicRef();
   SparseStorageDynamic& matrix_transp = mat->getStorageDynamicTransposedRef();

   assert(0 <= col && col < matrix_transp.m);

   /* remove all entries in column from the sparse storage dynamic */
   for( int j = matrix_transp.rowptr[col].start; j < matrix_transp.rowptr[col].end; j++ )
   {
      const int row = matrix_transp.jcolM[j];
      const double coeff = matrix_transp.M[j];

      assert( !PIPSisEQ(0.0, coeff) );

      if( TRACK_ROW(node, row, system_type, linking) )
      {
         std::cout << "TRACKING_ROW: fixation of column " << col << " in tracked row"<< std::endl;
         writeRowLocalToStreamDense(std::cout, system_type, node, linking, row);
      }

      /* remove the entry, adjust activity and row counters and rhs/lhs */
      removeEntryInDynamicStorage(matrix, row, col);

      reduceNnzCounterRow(system_type, node, linking, row, 1);

      adjustMatrixRhsLhsBy(system_type, node, linking, row, - coeff * fixation);

      adjustRowActivityFromDeletion(system_type, node, block_type, row, col, coeff);

      if( TRACK_ROW(node, row, system_type, linking) )
      {
         std::cout << "TRACKING_ROW: after removal of column" << std::endl;
         writeRowLocalToStreamDense(std::cout, system_type, node, block_type, row);

         double act_min = (system_type == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmin_eq_part, node, linking)[row] :
               getSimpleVecFromRowStochVec(*actmin_ineq_part, node, linking)[row];
         double act_max = (system_type == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmax_eq_part, node, linking)[row] :
               getSimpleVecFromRowStochVec(*actmax_ineq_part, node, linking)[row];

         int act_min_ubndd = (system_type == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmin_eq_ubndd, node, linking)[row] :
               getSimpleVecFromRowStochVec(*actmin_ineq_ubndd, node, linking)[row];
         int act_max_ubndd = (system_type == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmax_eq_ubndd, node, linking)[row] :
               getSimpleVecFromRowStochVec(*actmax_ineq_ubndd, node, linking)[row];

         std::cout << "TRACKING_ROW: New activity of row " << row << std::endl;
         std::cout << "\tnew min/max activity is: " << act_min << "/" << act_max << ", min/max unbounded counters are " << act_min_ubndd << "/" << act_max_ubndd << std::endl;
      }
   }

   /* adjust column counters */
   reduceNnzCounterColumn(node, block_type, col, matrix_transp.rowptr[col].end - matrix_transp.rowptr[col].start);

   /* update the transposed */
   matrix_transp.rowptr[col].end = matrix_transp.rowptr[col].start;

   // todo assert(transposed and normal matrix are in sync)
}

// todo
void PresolveData::removeParallelRow(SystemType system_type, int node, int row, bool linking)
{
   if( TRACK_ROW(node, row, system_type, linking) )
   {
      std::cout << "TRACKING_ROW: removal of tracked row as parallel row" << std::endl;
   }

   throw std::runtime_error("Not yet implemented");
//   if(postsolver)
//      postsolver->notifyParallelRow()

   removeRow(system_type, node, row, linking);
}

/* a singleton variable is substituted out of the problem and then it's original row can be removed from the problem */
void PresolveData::substituteVariableParallelRows(SystemType system_type, int node, int var1, int row1, int node_var1, int var2, int row2, int node_var2,
   double scalar, double translation)
{
   // todo : track row
   
   postsolver->notifyParallelRowSubstitution(system_type, node, var1, row1, node_var1, var2, row2, node_var2, scalar, translation);

   // delete the equality constraint which contained var2 (the substituted variable)
   removeRedundantRow( system_type, node, row2, false);
   assert( PIPSisZero(getSimpleVecFromColStochVec(*nnzs_col, node_var2)[var2]) );
   
   const double obj_var2 = getSimpleVecFromColStochVec(*presProb->g, node_var2)[var2];
   const double val_offset = translation * obj_var2;
   const double change_obj_var1 = scalar * obj_var2;

   removeColumn( node_var2, var2, 0.0 );

   if( node_var1 != -1 )
   {
      getSimpleVecFromColStochVec(*presProb->g, node)[var1] += change_obj_var1;
      obj_offset_chgs += val_offset;
   }
   else if( node == -1 )  
   {
      /* parallel rows in parent block - all processes should have detected this */
      getSimpleVecFromColStochVec(*presProb->g, -1)[var1] += change_obj_var1;

      // only add the objective offset for root as process ZERO:
      if( my_rank == 0 ) 
         obj_offset_chgs += val_offset;
   }
   else
   {
      /* var1 is a linking variable - store objective adaptions and allreduce them */
      (*objective_vec_chgs)[var1] += change_obj_var1;
      outdated_obj_vector = true; 
      obj_offset_chgs += val_offset;
   }
}

void PresolveData::removeRedundantRow(SystemType system_type, int node, int row, bool linking)
{
   if(postsolver)
   {
      // todo : adjust interface here
      std::vector<int> idx_row;
      std::vector<double> val_row;

      postsolver->notifyRedundantRow(system_type, node, row, linking, idx_row, val_row);
   }
 
#ifdef TRACK_ROW
   if( TRACK_ROW(node, row, system_type, linking) )
   {
      std::cout << "TRACKING_ROW: removal of tracked row as redundant row" << std::endl;
   }
#endif

   removeRow(system_type, node, row, linking);
}

/* removes row from local system - sets rhs lhs and activities to zero */
void PresolveData::removeRow(SystemType system_type, int node, int row, bool linking)
{
   assert(-1 <= node && node < nChildren);
   assert(!nodeIsDummy(node));

   if(linking)
   {
      assert(node == -1);

      /* Bl0 */
      removeRowFromMatrix(system_type, -1, BL_MAT, row);

      /* linking rows Bli */
      for(int child = 0; child < nChildren; ++child)
      {
         if(!nodeIsDummy(child))
            removeRowFromMatrix(system_type, child, BL_MAT, row);
      }
   }
   else
   {
      /* Bmat */
      removeRowFromMatrix(system_type, node, B_MAT, row);

      /* Amat */
      if(node != -1)
         removeRowFromMatrix(system_type, node, A_MAT, row);
   }


   /* set lhs rhs to zero */
   if(system_type == EQUALITY_SYSTEM)
      getSimpleVecFromRowStochVec(*presProb->bA, node, linking)[row] = 0.0;
   else
   {
      getSimpleVecFromRowStochVec(*presProb->bl, node, linking)[row] = 0.0;
      getSimpleVecFromRowStochVec(*presProb->bu, node, linking)[row] = 0.0;
   }

   /* set activities and unbounded counters to zero */
   if(system_type == EQUALITY_SYSTEM)
   {
      getSimpleVecFromRowStochVec(*actmax_eq_part, node, linking)[row] = 0.0;
      getSimpleVecFromRowStochVec(*actmin_eq_part, node, linking)[row] = 0.0;
      getSimpleVecFromRowStochVec(*actmax_eq_ubndd, node, linking)[row] = 0;
      getSimpleVecFromRowStochVec(*actmin_eq_ubndd, node, linking)[row] = 0;

      if(linking)
      {
         (*actmax_eq_chgs)[row] = 0.0;
         (*actmin_eq_chgs)[row] = 0.0;
         (*actmax_eq_ubndd_chgs)[row] = 0;
         (*actmin_eq_ubndd_chgs)[row] = 0;
      }
   }
   else
   {
      getSimpleVecFromRowStochVec(*actmax_ineq_part, node, linking)[row] = 0.0;
      getSimpleVecFromRowStochVec(*actmin_ineq_part, node, linking)[row] = 0.0;
      getSimpleVecFromRowStochVec(*actmax_ineq_ubndd, node, linking)[row] = 0;
      getSimpleVecFromRowStochVec(*actmin_ineq_ubndd, node, linking)[row] = 0;

      if(linking)
      {
         (*actmax_ineq_chgs)[row] = 0.0;
         (*actmin_ineq_chgs)[row] = 0.0;
         (*actmax_ineq_ubndd_chgs)[row] = 0;
         (*actmin_ineq_ubndd_chgs)[row] = 0;
      }
   }


#ifndef NDEBUG
   /* assert non-zero counters of row are zero - only works for non-linking rows */
   if(system_type == EQUALITY_SYSTEM)
   {
      if(!linking)
         assert( getSimpleVecFromRowStochVec(*nnzs_row_A, node, linking)[row] == 0 );   
   }
   else
   {
      if(!linking)
         assert( getSimpleVecFromRowStochVec(*nnzs_row_C, node, linking)[row] == 0 );
   }  
#endif
}

void PresolveData::removeRowFromMatrix(SystemType system_type, int node, BlockType block_type, int row)
{
   assert(!nodeIsDummy(node));
   assert(-1 <= node && node < nChildren);
   assert(node != -1 || block_type != A_MAT);
   SparseGenMatrix* mat = getSparseGenMatrix(system_type, node, block_type);

   assert(mat);
   assert(mat->hasDynamicStorage());

   SparseStorageDynamic* mat_storage = mat->getStorageDynamic();
   SparseStorageDynamic* mat_transp_storage = mat->getStorageDynamicTransposed();
   assert(mat_storage);
   assert( 0 <= row && row < mat_storage->m);
   assert(mat_transp_storage);

   const int row_start = mat_storage->rowptr[row].start;
   const int row_end = mat_storage->rowptr[row].end;

   for(int k = row_start; k < row_end; k++)
   {
      const int col_idx = mat_storage->jcolM[k];

      removeEntryInDynamicStorage(*mat_transp_storage, col_idx, row);
      reduceNnzCounterColumn(node, block_type, col_idx, 1);
   }
   
   const bool linking = (block_type == BL_MAT);

   reduceNnzCounterRow(system_type, node, linking, row, mat_storage->rowptr[row].end - mat_storage->rowptr[row].start);
   mat_storage->rowptr[row].end = mat_storage->rowptr[row].start;
}

/** removes row col from dynamic storage */
void PresolveData::removeEntryInDynamicStorage(SparseStorageDynamic& storage, int row, int col) const
{
   int i = -1;
   int end = storage.rowptr[row].end;
   int start = storage.rowptr[row].start;

   for( i = start; i < end; i++)
   {
      if( storage.jcolM[i] == col )
         break;
   }
   assert( storage.jcolM[i] == col);

   std::swap(storage.M[i], storage.M[end-1]);
   std::swap(storage.jcolM[i], storage.jcolM[end-1]);
   storage.rowptr[row].end--;
}

bool PresolveData::verifyActivities()
{
   assert(!outdated_activities && linking_rows_need_act_computation == 0);

   bool activities_correct = true;

   StochVectorHandle actmax_eq_part_old(dynamic_cast<StochVector*>(actmax_eq_part->cloneFull()));
   StochVectorHandle actmin_eq_part_old(dynamic_cast<StochVector*>(actmin_eq_part->cloneFull()));

   StochVectorBaseHandle<int> actmax_eq_ubndd_old(dynamic_cast<StochVectorBase<int>*>(actmax_eq_ubndd->cloneFull()));
   StochVectorBaseHandle<int> actmin_eq_ubndd_old(dynamic_cast<StochVectorBase<int>*>(actmin_eq_ubndd->cloneFull()));

   StochVectorHandle actmax_ineq_part_old(dynamic_cast<StochVector*>(actmax_ineq_part->cloneFull()));
   StochVectorHandle actmin_ineq_part_old(dynamic_cast<StochVector*>(actmin_ineq_part->cloneFull()));

   StochVectorBaseHandle<int> actmax_ineq_ubndd_old(dynamic_cast<StochVectorBase<int>*>(actmax_ineq_ubndd->cloneFull()));
   StochVectorBaseHandle<int> actmin_ineq_ubndd_old(dynamic_cast<StochVectorBase<int>*>(actmin_ineq_ubndd->cloneFull()));

   actmax_eq_part->setToZero();
   actmin_eq_part->setToZero();

   actmax_eq_ubndd->setToZero();
   actmin_eq_ubndd->setToZero();

   actmax_ineq_part->setToZero();
   actmin_ineq_part->setToZero();

   actmax_ineq_ubndd->setToZero();
   actmin_ineq_ubndd->setToZero();

   recomputeActivities();

   if( !actmax_eq_part_old->componentEqual(*actmax_eq_part, feastol))
   {
      if(my_rank == 0)
         std::cout << "actmax_eq_part not correct" << std::endl;
      activities_correct = false;
   }

   if( !actmin_eq_part_old->componentEqual(*actmin_eq_part, feastol))
   {
      if(my_rank == 0)
         std::cout << "actmin_eq_part not correct" << std::endl;
      activities_correct = false;
   }

   if( !actmax_eq_ubndd_old->componentEqual(*actmax_eq_ubndd, feastol))
   {
      if(my_rank == 0)
         std::cout << "actmax_eq_ubndd not correct" << std::endl;
      activities_correct = false;
   }

   if( !actmin_eq_ubndd_old->componentEqual(*actmin_eq_ubndd, feastol))
   {
      if(my_rank == 0)
         std::cout << "actmin_eq_ubndd not correct" << std::endl;
      activities_correct = false;
   }

   if( !actmax_ineq_part_old->componentEqual(*actmax_ineq_part, feastol))
   {
      if(my_rank == 0)
         std::cout << "actmax_ineq_part not correct" << std::endl;
      activities_correct = false;
   }

   if( !actmin_ineq_part_old->componentEqual(*actmin_ineq_part, feastol))
   {
      if(my_rank == 0)
         std::cout << "actmin_ineq_part not correct" << std::endl;
      activities_correct = false;
   }

   if( !actmax_ineq_ubndd_old->componentEqual(*actmax_ineq_ubndd, feastol))
   {
      if(my_rank == 0)
         std::cout << "actmax_ineq_ubndd not correct" << std::endl;
      activities_correct = false;
   }

   if( !actmin_ineq_ubndd_old->componentEqual(*actmin_ineq_ubndd, feastol))
   {
      if(my_rank == 0)
         std::cout << "actmin_ineq_ubndd not correct" << std::endl;
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

   assert(dynamic_cast<StochGenMatrix&>(*presProb->A).children[node]->isKindOf(kStochGenDummyMatrix)
		   == dynamic_cast<StochGenMatrix&>(*presProb->C).children[node]->isKindOf(kStochGenDummyMatrix));

   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*presProb->A);
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
   if( system_type == EQUALITY_SYSTEM )
   {
      dynamic_cast<StochGenMatrix&>(*(presProb->A)).Blmat->getSize(mlink, nlink);
      if( mlink > 0 )
      {
         // todo: assert that all vectors and matrices have linking part
         return true;
      }
   }
   else
   {
      dynamic_cast<StochGenMatrix&>(*(presProb->C)).Blmat->getSize(mlink, nlink);
      if( mlink > 0 )
      {
         // todo: assert that all vectors and matrices have linking part
         return true;
      }
   }
   return false;
}

/** adjusts unbounded counters of row as well as activity (if applicable) */ // todo refactor
void PresolveData::adjustRowActivityFromDeletion(SystemType system_type, int node, BlockType block_type, int row, int col, double coeff)
{
   assert(-1 <= node && node < nChildren);
   assert(0 <= col);
   assert(0 <= row);
   assert( !PIPSisZero(coeff) );
   assert( !nodeIsDummy(node) );

   /* get upper and lower bound on variable */
   const SimpleVector& ixlow = getSimpleVecFromColStochVec(*(presProb->ixlow), (block_type == A_MAT) ? -1 : node);
   const SimpleVector& xlow = getSimpleVecFromColStochVec(*(presProb->blx), (block_type == A_MAT) ? -1 : node);
   const SimpleVector& ixupp = getSimpleVecFromColStochVec(*(presProb->ixupp), (block_type == A_MAT) ? -1 : node);
   const SimpleVector& xupp = getSimpleVecFromColStochVec(*(presProb->bux), (block_type == A_MAT) ? -1 : node);

   const bool linking = ( block_type == BL_MAT );

   assert(PIPSisEQ(ixlow[col], 1.0));
   assert(PIPSisEQ(ixupp[col], 1.0));
   assert(PIPSisEQ(xlow[col], xupp[col], 1e-10));

   /* get unbounded counters */
   int* actmax_ubndd = (system_type == EQUALITY_SYSTEM) ? &getSimpleVecFromRowStochVec(*actmax_eq_ubndd, node, linking)[row]
         : &getSimpleVecFromRowStochVec(*actmax_ineq_ubndd, node, linking)[row];
   int* actmin_ubndd = (system_type == EQUALITY_SYSTEM) ? &getSimpleVecFromRowStochVec(*actmin_eq_ubndd, node, linking)[row]
         : &getSimpleVecFromRowStochVec(*actmin_ineq_ubndd, node, linking)[row];

   /* sum of reduction and counters */
   int actmax_ubndd_val = *actmax_ubndd;
   int actmin_ubndd_val = *actmin_ubndd;

   if( linking )
   {
      actmax_ubndd = (system_type == EQUALITY_SYSTEM) ? &(*actmax_eq_ubndd_chgs)[row]
            : &(*actmax_ineq_ubndd_chgs)[row];
      actmin_ubndd = (system_type == EQUALITY_SYSTEM) ? &(*actmin_eq_ubndd_chgs)[row]
            : &(*actmin_ineq_ubndd_chgs)[row];

      actmax_ubndd_val += *actmax_ubndd;
      actmin_ubndd_val += *actmin_ubndd;
   }

   assert(0 <= actmax_ubndd_val);
   assert(0 <= actmin_ubndd_val);

   /* get partial activities */
   double* actmax_part = (system_type == EQUALITY_SYSTEM) ? &getSimpleVecFromRowStochVec(*actmax_eq_part, node, linking)[row]
         : &getSimpleVecFromRowStochVec(*actmax_ineq_part, node, linking)[row];
   double* actmin_part = (system_type == EQUALITY_SYSTEM) ? &getSimpleVecFromRowStochVec(*actmin_eq_part, node, linking)[row]
         : &getSimpleVecFromRowStochVec(*actmin_ineq_part, node, linking)[row];

   if( linking )
   {
      actmax_part = (system_type == EQUALITY_SYSTEM) ? &(*actmax_eq_chgs)[row] : &(*actmax_ineq_chgs)[row];
      actmin_part = (system_type == EQUALITY_SYSTEM) ? &(*actmin_eq_chgs)[row] : &(*actmin_ineq_chgs)[row];
   }

   assert(actmin_part);
   assert(actmax_part);
   assert(actmin_ubndd);
   assert(actmax_ubndd);
   assert(col < ixlow.n);

   /* depending on the sign of coeff add / substract lbx * coeff/ ubx * coeff from actmin and actmax */
   if( PIPSisLT(0.0, coeff) )
   {
      if( PIPSisEQ(ixupp[col], 1.0) )
         (*actmax_part) -= coeff * xupp[col];
      else
      {
         assert(1 <= actmax_ubndd_val);
         --(*actmax_ubndd);
         if( actmax_ubndd_val == 2)
            computeRowMinOrMaxActivity(system_type, node, linking, row, true);
      }

      if( PIPSisEQ(ixlow[col], 1.0) )
         (*actmin_part) -= coeff * xlow[col];
      else
      {
         assert(1 <= actmin_ubndd_val);
         --(*actmin_ubndd);
         if( actmin_ubndd_val == 2)
            computeRowMinOrMaxActivity(system_type, node, linking, row, false);
      }
   }
   else
   {
      if( PIPSisEQ(ixlow[col], 1.0) )
         (*actmax_part) -= coeff*xlow[col];
      else
      {
         assert(1 <= actmax_ubndd_val);
         --(*actmax_ubndd);
         if( actmax_ubndd_val == 2)
            computeRowMinOrMaxActivity(system_type, node, linking, row, true);
      }

      if( PIPSisEQ(ixupp[col], 1.0) )
         (*actmin_part) -= coeff*xupp[col];
      else
      {
         assert(1 <= actmin_ubndd_val);
         --(*actmin_ubndd);
         if( actmin_ubndd_val == 2)
            computeRowMinOrMaxActivity(system_type, node, linking, row, false);
      }
   }

   if( linking )
      outdated_activities = true;
}

/// computes min or max activity of linking row regardless of unbounded counters
double PresolveData::computeLocalLinkingRowMinOrMaxActivity(SystemType system_type, int row, bool upper) const
{
   double act_part = 0.0;

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
      SparseStorageDynamic& mat = getSparseGenMatrix(system_type, node, BL_MAT)->getStorageDynamicRef();

      for( int j = mat.rowptr[row].start; j < mat.rowptr[row].end; j++ )
      {
         const int col = mat.jcolM[j];
         const double entry = mat.M[j];

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
void PresolveData::computeRowMinOrMaxActivity(SystemType system_type, int node, bool linking, int row, bool upper)
{
   assert( -1 <= node && node < nChildren);
   
   if( TRACK_ROW(node, row, system_type, linking) )
   {
      if( linking)
         std::cout << "TRACKING_ROW: tracked row is linking constraint and needs " << (upper ? "upper" : "lower") << " activity recomputation " << std::endl;
      else
         std::cout << "TRACKING_ROW: computing " << (upper ? "upper" : "lower") << " activity of tracked row " << std::endl;
   }

   /* single linking rows that have to be computed for the first time */
   if( linking )
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
   const SimpleVector& ixlow = getSimpleVecFromColStochVec(*(presProb->ixlow), node);
   const SimpleVector& xlow = getSimpleVecFromColStochVec(*(presProb->blx), node);
   const SimpleVector& ixupp = getSimpleVecFromColStochVec(*(presProb->ixupp), node);
   const SimpleVector& xupp = getSimpleVecFromColStochVec(*(presProb->bux), node);

   /* get matrix */
   SparseStorageDynamic& Bmat = getSparseGenMatrix(system_type, node, B_MAT)->getStorageDynamicRef();

   /* get activity vector */
   SimpleVector* act_vec;
   if(system_type == EQUALITY_SYSTEM)
   {
      act_vec = (!upper) ? &getSimpleVecFromRowStochVec(*actmin_eq_part, node, linking) :
            &getSimpleVecFromRowStochVec(*actmax_eq_part, node, linking);
   }
   else
   {
      act_vec = (!upper) ? &getSimpleVecFromRowStochVec(*actmin_ineq_part, node, linking) :
            &getSimpleVecFromRowStochVec(*actmax_ineq_part, node, linking);
   }

   double& act_part = (*act_vec)[row];
   act_part = 0;

   /* Bmat */
   for( int j = Bmat.rowptr[row].start; j < Bmat.rowptr[row].end; j++)
   {
      const int col = Bmat.jcolM[j];
      const double entry = Bmat.M[j];

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
   if( node != -1 )
   {
      SparseStorageDynamic& Amat = getSparseGenMatrix(system_type, node, A_MAT)->getStorageDynamicRef();
      for( int j = Amat.rowptr[row].start; j < Amat.rowptr[row].end; j++ )
      {
         const int col = Amat.jcolM[j];
         const double entry = Amat.M[j];

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

   if( TRACK_ROW(node, row, system_type, linking) )
         std::cout << "TRACKING_ROW: new " << (upper ? "upper" : "lower") << " activity is " << act_part << std::endl;
}

/** updates the bounds on a variable as well as activities */
bool PresolveData::updateBoundsVariable(int node, int col, double xu, double xl)
{
   SimpleVector& ixlow = getSimpleVecFromColStochVec(*(presProb->ixlow), node);
   SimpleVector& xlow = getSimpleVecFromColStochVec(*(presProb->blx), node);
   SimpleVector& ixupp = getSimpleVecFromColStochVec(*(presProb->ixupp), node);
   SimpleVector& xupp = getSimpleVecFromColStochVec(*(presProb->bux), node);

   if( TRACK_COLUMN(node, col) )
      std::cout << "TRACKING_COLUMN: updating column bounds from [" << xlow[col] << ", " << 
         xupp[col] << "] for column " << col << " node " << node << " with [" <<
         xl << ", " << xu << "]" << std::endl;

   if( ( !PIPSisZero(ixlow[col]) && PIPSisLT(xu, xlow[col]) )
         || ( !PIPSisZero(ixupp[col]) && PIPSisLT(xupp[col], xl) )
         || PIPSisLT(xu, xl) )
      PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Varbounds update detected infeasible new bounds!", "PresolveData.C", "updateBoundsVariable");

   bool updated = false;
   double xu_old = std::numeric_limits<double>::infinity();
   double xl_old = -std::numeric_limits<double>::infinity();

   if( xu < std::numeric_limits<double>::max() && ( PIPSisZero(ixupp[col]) || PIPSisLT(xu, xupp[col])) )
   {
      updated = true;
      if( PIPSisEQ(ixupp[col], 1.0) )
         xu_old = xupp[col];

      xupp[col] = xu;
      ixupp[col] = 1.0;
   }

   if( -std::numeric_limits<double>::max() < xl && ( PIPSisZero(ixlow[col]) || PIPSisLE( xlow[col], xl)) )
   {
      updated = true;
      if( PIPSisEQ(ixlow[col], 1.0) )
         xl_old = xlow[col];

      xlow[col] = xl;
      ixlow[col] = 1.0;
   }

   if( updated )
   {
      if(  TRACK_COLUMN(node, col) )
      {
         std::cout << "TRACKING_COLUMN: bounds are now [" << xlow[col] << ", " << xupp[col] << "]" << std::endl;
         std::cout << "TRACKING_COLUMN: moving on to update activities" << std::endl;
      }
      updateRowActivities(node, col, xu, xl, xu_old, xl_old);
   }
   else
   {
      if( TRACK_COLUMN(node, col) )
         std::cout << "TRACKING_COLUMN: col " << col << " was not updated" << std::endl;
   }


   if( updated && node == -1 )
      outdated_linking_var_bounds = true;

   return updated;
}

/** goes through given column and adjusts activities by a_ij * (old_ubx - ubx)
 *  If the variable previously was unbounded in upper/lower direction unbounded_counters
 *  get adjusted and the row activity is computed for the first time if necessary
 */
void PresolveData::updateRowActivities(int node, int col, double ubx, double lbx, double old_ubx, double old_lbx)
{
   assert(-1 <= node && node < nChildren);

   if( TRACK_COLUMN(node,col) )
   {
      std::cout << "TRACKING_COLUMN: col " << col << " node " << node << " gets it's activities updated" << std::endl;
      std::cout << "\t bounds changed from [" << old_ubx << ", " << old_lbx << "] to [" << lbx << ", " << ubx << "]" << std::endl;
   }

   /* if node == -1 go through all linking var blocks of both systems */
   if( node == -1 )
   {
      if( TRACK_COLUMN(node,col) )
         std::cout << "TRACKING_COLUMN: the column is linking (root node)" << std::endl;
      /* A0/B0 and C0/D0block */
      updateRowActivitiesBlock(EQUALITY_SYSTEM, -1, B_MAT, col, ubx, lbx, old_ubx, old_lbx);
      updateRowActivitiesBlock(INEQUALITY_SYSTEM, -1, B_MAT, col, ubx, lbx, old_ubx, old_lbx);

      /* Bl0 and Dl0 */
      updateRowActivitiesBlock(EQUALITY_SYSTEM, -1, BL_MAT, col, ubx, lbx, old_ubx, old_lbx);
      updateRowActivitiesBlock(INEQUALITY_SYSTEM, -1, BL_MAT, col, ubx, lbx, old_ubx, old_lbx);

      for(int child = 0; child < nChildren; ++child)
      {
         /* Ai and Ci */
         updateRowActivitiesBlock(EQUALITY_SYSTEM, child, A_MAT, col, ubx, lbx, old_ubx, old_lbx);
         updateRowActivitiesBlock(INEQUALITY_SYSTEM, child, A_MAT, col, ubx, lbx, old_ubx, old_lbx);
      }
   }
   else
   {
      if( TRACK_COLUMN(node,col) )
         std::cout << "TRACKING_COLUMN: the column is non-linking (non-root)" << std::endl;
      /* Bmat, Blmat */
      /* Bmat */
      updateRowActivitiesBlock(EQUALITY_SYSTEM, node, B_MAT, col, ubx, lbx, old_ubx, old_lbx);

      /* Blmat */
      updateRowActivitiesBlock(EQUALITY_SYSTEM, node, BL_MAT, col, ubx, lbx, old_ubx, old_lbx);

      /* Dmat Dlmat */

      /* Dmat */
      updateRowActivitiesBlock(INEQUALITY_SYSTEM, node, B_MAT, col, ubx, lbx, old_ubx, old_lbx);

      /* Dlmat */
      updateRowActivitiesBlock(INEQUALITY_SYSTEM, node, BL_MAT, col, ubx, lbx, old_ubx, old_lbx);
   }
}


void PresolveData::updateRowActivitiesBlock(SystemType system_type, int node, BlockType block_type, int col,
		 double ubx, double lbx, double old_ubx, double old_lbx)
{
	updateRowActivitiesBlock(system_type, node, block_type, col, lbx, old_lbx, false);
	updateRowActivitiesBlock(system_type, node, block_type, col, ubx, old_ubx, true);
}


void PresolveData::updateRowActivitiesBlock(SystemType system_type, int node, BlockType block_type, int col, double bound, double old_bound, bool upper)
{
   assert(-1 <= node && node < nChildren);
   assert(0 <= col);

   /* dummies do not adjust anything */
   if( nodeIsDummy(node) )
      return;

   /* we do not have to adjust activities if no new bound was found */
   if( bound == std::numeric_limits<double>::infinity() || bound == -std::numeric_limits<double>::infinity() )
      return;

   /* no changes if a worse bound has been found */
   if( (upper && !PIPSisLT( bound, old_bound)) || (!upper && PIPSisLT(bound, old_bound)) )
      return;

   SparseStorageDynamic& mat_transp = getSparseGenMatrix(system_type, node, block_type)->getStorageDynamicTransposedRef();
   const bool linking = ( block_type == BL_MAT );

   assert(col < mat_transp.m);

   SimpleVectorBase<int>& actmax_ubndd = (system_type == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmax_eq_ubndd, node, linking) :
         getSimpleVecFromRowStochVec(*actmax_ineq_ubndd, node, linking);
   SimpleVectorBase<int>& actmin_ubndd = (system_type == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmin_eq_ubndd, node, linking) :
         getSimpleVecFromRowStochVec(*actmin_ineq_ubndd, node, linking);
   SimpleVector& actmax_part = (system_type == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmax_eq_part, node, linking) :
         getSimpleVecFromRowStochVec(*actmax_ineq_part, node, linking);
   SimpleVector& actmin_part = (system_type == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmin_eq_part, node, linking) :
         getSimpleVecFromRowStochVec(*actmin_ineq_part, node, linking);

   /* we always set these variables but only use them in case an actual linking row gets it's activities updated */
   SimpleVectorBase<int>& actmax_ubndd_chgs = (system_type == EQUALITY_SYSTEM) ? *actmax_eq_ubndd_chgs : *actmax_ineq_ubndd_chgs;
   SimpleVectorBase<int>& actmin_ubndd_chgs = (system_type == EQUALITY_SYSTEM) ? *actmin_eq_ubndd_chgs : *actmin_ineq_ubndd_chgs;
   SimpleVector& actmax_part_chgs = (system_type == EQUALITY_SYSTEM) ? *actmax_eq_chgs : *actmax_ineq_chgs;
   SimpleVector& actmin_part_chgs = (system_type == EQUALITY_SYSTEM) ? *actmin_eq_chgs : *actmin_ineq_chgs;

   if( TRACK_COLUMN(node, col) )
      std::cout << "TRACKING_COLUMN: updating activities column " << col << " node " << node << " system " << ( (system_type == EQUALITY_SYSTEM) ? "EQ_SYS" : "INEQ_SYS" ) <<
         " with new " << ( upper ? "upper" : "lower" ) << " bound " << bound << " from " << old_bound << std::endl;

   for( int j = mat_transp.rowptr[col].start; j < mat_transp.rowptr[col].end; ++j )
   {
      const int row = mat_transp.jcolM[j];
      const double entry = mat_transp.M[j];

      assert( !PIPSisZero(entry) );

      /* get affected partial activity and act_ubndd */
      bool switch_upperlower = upper;
      if( PIPSisLE(entry, 0.0) )
         switch_upperlower = !upper;

      SimpleVector& act_part = (switch_upperlower) ? actmax_part : actmin_part;
      SimpleVector& act_part_chgs = (switch_upperlower) ? actmax_part_chgs : actmin_part_chgs;
      SimpleVectorBase<int>& act_ubndd = (switch_upperlower) ? actmax_ubndd : actmin_ubndd;
      SimpleVectorBase<int>& act_ubndd_chgs = (switch_upperlower) ? actmax_ubndd_chgs : actmin_ubndd_chgs;


      if( TRACK_ROW(node, row, system_type, linking) )
      {
         std::cout << "TRACKING_ROW: activities of tracked row are getting updated. Current state:" << std::endl;
         if( linking)
         {
            std::cout << "\twas linking row with " << ( switch_upperlower ? "upper" : "lower" ) << " activity " << act_part[row] << " + " << act_part_chgs[row] << "(changes) and " <<
                  (switch_upperlower ? "upper" : "lower" ) << " unbounded counters " << act_ubndd[row] << " + " << act_ubndd_chgs[row] << "(changes)" << std::endl;
         }
         else
         {
            std::cout << "\twas non-linking row with " << ( switch_upperlower ? "upper" : "lower" ) << " activity " << act_part[row] << " and " << (switch_upperlower ? "upper" : "lower" ) <<
                  " unbounded counters " << act_ubndd[row] << std::endl;
         }
      }

      /* if the old bound was not set we have to modify the unbounded counters */
      if( old_bound == std::numeric_limits<double>::infinity() || old_bound == -std::numeric_limits<double>::infinity() )
      {
         if( TRACK_ROW(node, row, system_type, block_type == linking ) )
               std::cout << "TRACKING_ROW: ubndd counters are being changed" << std::endl;

         /* every process works the root node - even in the linking row case */
         if( node != -1 && linking  )
         {
            --act_ubndd_chgs[row];
            outdated_activities = true;
         }
         else
            --act_ubndd[row];

         if(linking )
            assert( 0 <= act_ubndd[row] + act_ubndd_chgs[row]);
         else
            assert( 0 <= act_ubndd[row]);

         /* from now on activity of row has to be computed */
         if( (linking  && (act_ubndd[row] + act_ubndd_chgs[row] == 1))
               || (act_ubndd[row] == 1 && !linking) )
         {

            if( TRACK_ROW(node, row, system_type, linking ) )
            {
               if( linking  )
                  std::cout << "TRACKING_ROW: " << ( switch_upperlower ? "upper" : "lower" ) << " now needs act computation " << std::endl;
               else
                  std::cout << "TRACKING_ROW: first time computation of " << ( switch_upperlower ? "upper" : "lower" ) << " activity " << std::endl;
            }

            if(linking )
               ++linking_rows_need_act_computation;
            else
               computeRowMinOrMaxActivity(system_type, node, linking, row, switch_upperlower);
         }
         else if( (linking  && act_ubndd[row] + act_ubndd_chgs[row] == 0) ||
               (!linking && act_ubndd[row] == 0) )
         {
            if( TRACK_ROW(node, row, system_type, linking ) )
                  std::cout << "TRACKING_ROW: adjusting activities " << std::endl;

            if( node != -1 && linking )
               act_part_chgs[row] += bound * entry;
            else
               act_part[row] += bound * entry;
         }
      }
      else
      {
            if( TRACK_ROW(node, row, system_type, linking ) )
                  std::cout << "TRACKING_ROW: " << " no changes in ubndd counters - only adjust activities" << std::endl;

         /* better bound was found */
         if( (linking  && act_ubndd[row] + act_ubndd_chgs[row] <= 1) ||
               (!linking && act_ubndd[row] <= 1) )
         {

            if( node != -1 && linking )
            {
               outdated_activities = true;
               act_part_chgs[row] += (bound - old_bound) * entry;
            }
            else
               act_part[row] += (bound - old_bound) * entry;
         }
      }

      if( TRACK_ROW(node, row, system_type, linking ) )
      {
         std::cout << "TRACKING_ROW: activities of tracked row have been updated: " << std::endl;
         if( linking )
         {
            std::cout << "\tnow linking row with " << ( switch_upperlower ? "upper" : "lower" ) <<
               " activity " << act_part[row] << " + " << act_part_chgs[row] << "(changes) and " <<
               (switch_upperlower ? "upper" : "lower" ) << " unbounded counters " << 
               act_ubndd[row] << " + " << act_ubndd_chgs[row] << "(changes)" << std::endl;
         }
         else
         {
            std::cout << "\tnow non-linking row with " << ( switch_upperlower ? "upper" : "lower" ) <<
            " activity " << act_part[row] << " and " << (switch_upperlower ? "upper" : "lower" ) <<
            " unbounded counters " << act_ubndd[row] << std::endl;
         }
      }
   }
}

/// returns activity and unbounded counters for specified row
/// +/- infinity if there are two or more unbouded entries in a row 
void PresolveData::getRowActivities(SystemType system_type, int node, bool linking, int row, double& max_act,
      double& min_act, int& max_ubndd, int& min_ubndd) const
{
   assert(!nodeIsDummy(node));

   max_ubndd = (system_type == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmax_eq_ubndd, node, linking)[row]
         : getSimpleVecFromRowStochVec(*actmax_ineq_ubndd, node, linking)[row];
   min_ubndd = (system_type == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmin_eq_ubndd, node, linking)[row]
         : getSimpleVecFromRowStochVec(*actmin_ineq_ubndd, node, linking)[row];

   if( linking )
   {
      max_ubndd += (system_type == EQUALITY_SYSTEM) ? (*actmax_eq_ubndd_chgs)[row]
            : (*actmax_ineq_ubndd_chgs)[row];
      min_ubndd += (system_type == EQUALITY_SYSTEM) ? (*actmin_eq_ubndd_chgs)[row]
            : (*actmin_ineq_ubndd_chgs)[row];
   }

   max_act = (system_type == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmax_eq_part, node, linking)[row]
         : getSimpleVecFromRowStochVec(*actmax_ineq_part, node, linking)[row];
   min_act = (system_type == EQUALITY_SYSTEM) ? getSimpleVecFromRowStochVec(*actmin_eq_part, node, linking)[row]
         : getSimpleVecFromRowStochVec(*actmin_ineq_part, node, linking)[row];

   if( linking )
   {
      max_act += (system_type == EQUALITY_SYSTEM) ? (*actmax_eq_chgs)[row] : (*actmax_ineq_chgs)[row];
      min_act += (system_type == EQUALITY_SYSTEM) ? (*actmin_eq_chgs)[row] : (*actmin_ineq_chgs)[row];
   }

   if( max_ubndd >= 2)
      assert( max_act == std::numeric_limits<double>::infinity() );
   else
   {
      if( !linking )
      {
         assert(max_act < std::numeric_limits<double>::infinity());
      }
   }

   if( min_ubndd >= 2)
      assert( min_act == -std::numeric_limits<double>::infinity() );
   else
   {
      if( !linking )
      {
         assert(-std::numeric_limits<double>::infinity() < min_act);
      }
   }
}

SparseGenMatrix* PresolveData::getSparseGenMatrix(SystemType system_type, int node, BlockType block_type) const
{
   assert( -1 <= node && node < nChildren );
   assert(!nodeIsDummy(node));

   StochGenMatrix& sMat = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<StochGenMatrix&>(*presProb->A) : dynamic_cast<StochGenMatrix&>(*presProb->C);

   if(node == -1)
   {
      if(block_type == BL_MAT)
         return sMat.Blmat;
      else 
      {
         assert(block_type == B_MAT);
         return sMat.Bmat;
      }
   }
   else
   {
      if(block_type == A_MAT)
         return sMat.children[node]->Amat;
      else if(block_type == B_MAT)
         return sMat.children[node]->Bmat;
      else if(block_type == BL_MAT)
         return sMat.children[node]->Blmat;
   }
   return NULL;
}

/* only prints the part of a linking constraint the current process knows about */
// todo does not yet print linking constraints
void PresolveData::writeRowLocalToStreamDense(std::ostream& out, SystemType system_type, int node, bool linking, int row) const
{
   if(nodeIsDummy(node))
      return;

   if(node == -1 && !linking && my_rank != 0)
      return;

   out << "SystemType: " << system_type << "\tnode: " << node << "\tLinkingCons: " << linking << "\trow: " << row << std::endl;

   if(system_type == EQUALITY_SYSTEM)
   {
      out << getSimpleVecFromRowStochVec(*presProb->bA, node, linking)[row] << " = ";
   }
   else
   {
      double iclow = getSimpleVecFromRowStochVec(*presProb->iclow, node, linking)[row];
      double clow = PIPSisEQ(iclow, 1.0) ? getSimpleVecFromRowStochVec(*presProb->bl, node, linking)[row] : -std::numeric_limits<double>::infinity();

      out << clow << " <= ";
   }

   if(node != -1 && !linking)
   {
      writeMatrixRowToStreamDense(out, *getSparseGenMatrix(system_type, node, A_MAT), node, row, getSimpleVecFromColStochVec(*presProb->ixupp, -1),
            getSimpleVecFromColStochVec(*presProb->bux, -1), getSimpleVecFromColStochVec(*presProb->ixlow, -1), getSimpleVecFromColStochVec(*presProb->blx, -1));
   }

   writeMatrixRowToStreamDense(out, *getSparseGenMatrix(system_type, node, B_MAT), node, row, getSimpleVecFromColStochVec(*presProb->ixupp, node),
         getSimpleVecFromColStochVec(*presProb->bux, node),getSimpleVecFromColStochVec(*presProb->ixlow, node),getSimpleVecFromColStochVec(*presProb->blx, node));


   if(system_type == INEQUALITY_SYSTEM)
   {
      double icupp = getSimpleVecFromRowStochVec(*presProb->icupp, node, linking)[row];
      double cupp = PIPSisEQ(icupp, 1.0) ? getSimpleVecFromRowStochVec(*presProb->bu, node, linking)[row] : std::numeric_limits<double>::infinity();

      out << " <= " << cupp;
   }
   out << std::endl;
}

void PresolveData::writeMatrixRowToStreamDense(std::ostream& out, const SparseGenMatrix& mat, int node, int row, const SimpleVector& ixupp, const SimpleVector& xupp,
      const SimpleVector& ixlow, const SimpleVector& xlow) const
{
   const SparseStorageDynamic& storage = mat.getStorageDynamicRef();

   const int start = storage.rowptr[row].start;
   const int end = storage.rowptr[row].end;

   for(int k = start; k < end; ++k)
   {
      const double val = storage.M[k];
      const int col = storage.jcolM[k];

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
