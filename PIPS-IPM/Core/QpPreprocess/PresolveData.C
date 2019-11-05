/*
 * PresolveData.C
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

// todo : find some useful lenth to initialize singleton rows and singleton cols to - prealloc some storage

#include "PresolveData.h"
#include "StochGenMatrix.h"
#include "DoubleMatrixTypes.h"
#include "pipsdef.h"

#include <stdexcept>
#include <limits>
#include <algorithm>
#include <string>
#include <cmath>

/* can specify a cloumn here which bound changes we wanna track maybe */

// #ifndef NDEBUG
//   #define TRACK_COLUMN
//   #define COLUMN 103327
//   #define NODE 102
// #endif

// #ifndef NDEBUG
//    #define TRACK_ROW
//    #define ROW 457762
//    #define ROW_NODE 28
//    #define ROW_BLOCK LINKING_VARS_BLOCK
//    #define ROW_IS_LINK false
//    #define ROW_SYS EQUALITY_SYSTEM
// #endif

PresolveData::PresolveData(const sData* sorigprob, StochPostsolver* postsolver) :
      postsolver(postsolver),
      outdated_lhsrhs(false),
      outdated_nnzs(false),
      outdated_linking_var_bounds(false),
      outdated_activities(true),
      outdated_obj_vector(false),
      linking_rows_need_act_computation(0),
      nnzs_row_A(dynamic_cast<StochVector*>(sorigprob->bA->clone())),
      nnzs_row_C(dynamic_cast<StochVector*>(sorigprob->icupp->clone())),
      nnzs_col(dynamic_cast<StochVector*>(sorigprob->g->clone())),
      actmax_eq_part(nnzs_row_A->clone()),
      actmin_eq_part(nnzs_row_A->clone()),
      actmax_eq_ubndd(nnzs_row_A->clone()),
      actmin_eq_ubndd(nnzs_row_A->clone()),
      actmax_ineq_part(nnzs_row_C->clone()),
      actmin_ineq_part(nnzs_row_C->clone()),
      actmax_ineq_ubndd(nnzs_row_C->clone()),
      actmin_ineq_ubndd(nnzs_row_C->clone()),
      obj_vector_chgs(dynamic_cast<SimpleVector*>(dynamic_cast<SimpleVector*>(nnzs_col->vec)->cloneFull())),
      nChildren(nnzs_col->children.size()),
      objOffset(0.0), obj_offset_chgs(0.0),
      objective_vec_chgs(dynamic_cast<SimpleVector*>(nnzs_col->vec->cloneFull())),
      lower_bound_implied_by_system(nnzs_col->clone()),
      lower_bound_implied_by_row(nnzs_col->clone()),
      lower_bound_implied_by_node(nnzs_col->clone()),
      upper_bound_implied_by_system(nnzs_col->clone()),
      upper_bound_implied_by_row(nnzs_col->clone()),
      upper_bound_implied_by_node(nnzs_col->clone()),
      elements_deleted(0), elements_deleted_transposed(0)
{
   lower_bound_implied_by_system->setToConstant(-10);
   lower_bound_implied_by_row->setToConstant(-10);
   lower_bound_implied_by_node->setToConstant(-10);
   upper_bound_implied_by_system->setToConstant(-10);
   upper_bound_implied_by_row->setToConstant(-10);
   upper_bound_implied_by_node->setToConstant(-10);
   obj_vector_chgs->setToZero();

   presProb = sorigprob->cloneFull(true);

   int n_linking_vars = (nnzs_col->vec) ? nnzs_col->vec->n : 0;

   int n_linking_A = (nnzs_row_A->vecl) ? nnzs_row_A->vecl->n : 0;
   int n_linking_C = (nnzs_row_C->vecl) ? nnzs_row_C->vecl->n : 0;

   assert( n_linking_vars + 2 * n_linking_A + 2 * n_linking_C < std::numeric_limits<int>::max() );

   length_array_nnz_chgs = n_linking_vars + n_linking_A + n_linking_C;
   array_nnz_chgs = new double[length_array_nnz_chgs]();

   nnzs_col_chgs = SimpleVectorHandle(new SimpleVector(array_nnz_chgs, n_linking_vars));
   nnzs_row_A_chgs = SimpleVectorHandle(new SimpleVector(array_nnz_chgs + n_linking_vars, n_linking_A));
   nnzs_row_C_chgs = SimpleVectorHandle(new SimpleVector(array_nnz_chgs + n_linking_vars + n_linking_A, n_linking_C));

   lenght_array_act_chgs = n_linking_A * 2 + n_linking_C * 2;
   array_act_chgs = new double[lenght_array_act_chgs]();

   actmax_eq_chgs = SimpleVectorHandle( new SimpleVector(array_act_chgs, n_linking_A));
   actmin_eq_chgs = SimpleVectorHandle( new SimpleVector(array_act_chgs + n_linking_A, n_linking_A));
   actmax_ineq_chgs = SimpleVectorHandle( new SimpleVector(array_act_chgs + 2 * n_linking_A, n_linking_C));
   actmin_ineq_chgs = SimpleVectorHandle( new SimpleVector(array_act_chgs + 2 * n_linking_A + n_linking_C, n_linking_C));

   array_act_unbounded_chgs = new double[lenght_array_act_chgs]();

   actmax_eq_ubndd_chgs = SimpleVectorHandle( new SimpleVector(array_act_unbounded_chgs, n_linking_A));
   actmin_eq_ubndd_chgs = SimpleVectorHandle( new SimpleVector(array_act_unbounded_chgs + n_linking_A, n_linking_A));
   actmax_ineq_ubndd_chgs = SimpleVectorHandle( new SimpleVector(array_act_unbounded_chgs + 2 * n_linking_A, n_linking_C));
   actmin_ineq_ubndd_chgs = SimpleVectorHandle( new SimpleVector(array_act_unbounded_chgs + 2 * n_linking_A + n_linking_C, n_linking_C));

   lenght_array_bound_chgs = n_linking_A + n_linking_C;
   array_bound_chgs = new double[lenght_array_bound_chgs]();

   bound_chgs_A = SimpleVectorHandle( new SimpleVector(array_bound_chgs, n_linking_A));
   bound_chgs_C = SimpleVectorHandle( new SimpleVector(array_bound_chgs + n_linking_A, n_linking_C));

   getRankDistributed(MPI_COMM_WORLD, my_rank, distributed);

   // initialize all dynamic transposed sub matrices
   dynamic_cast<StochGenMatrix&>(*presProb->A).initTransposed(true);
   dynamic_cast<StochGenMatrix&>(*presProb->C).initTransposed(true);

   recomputeActivities();
   initNnzCounter( *nnzs_row_A, *nnzs_row_C, *nnzs_col);

#ifdef TRACK_COLUMN
   if( (my_rank == 0 || NODE != -1) && !nodeIsDummy(NODE, EQUALITY_SYSTEM) )
   {
      std::cout << "inital column" << std::endl;
      std::cout << "ixlow: " << getSimpleVecColFromStochVec( *presProb->ixlow, NODE)[COLUMN] << std::endl;
      std::cout << "xlow: " << getSimpleVecColFromStochVec( *presProb->blx, NODE)[COLUMN] << std::endl;
      std::cout << "ixupp: " << getSimpleVecColFromStochVec( *presProb->ixupp, NODE)[COLUMN] << std::endl;
      std::cout << "xupp: " << getSimpleVecColFromStochVec( *presProb->bux, NODE)[COLUMN] << std::endl;
   }
#endif

#ifdef TRACK_ROW
   if(!nodeIsDummy(ROW_NODE, ROW_SYS))
   {
      std::cout << "TRACKING: row " << ROW << " node " << ROW_NODE << " in BlockType " << ROW_BLOCK << " SystemType " << ROW_SYS << std::endl;
      writeRowLocalToStreamDense(std::cout, ROW_SYS, ROW_NODE, ROW_BLOCK, ROW);
   }
#endif

   initSingletons();
   setUndefinedVarboundsTo(std::numeric_limits<double>::max());
}

PresolveData::~PresolveData()
{
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
      if( vec_ixlow[i] == 0.0)
         vec_xlow[i] = -value;

      if( vec_ixupp[i] == 0.0)
         vec_xupp[i] = value;
   }
}

sData* PresolveData::finalize()
{

#ifndef NDEBUG
   if(distributed)
   {
      // todo make one array?
      MPI_Allreduce(MPI_IN_PLACE, &outdated_activities, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_lhsrhs, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_nnzs, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_linking_var_bounds, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_obj_vector, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
   }
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

      SimpleVector& actmin_eq_root_ubndd = dynamic_cast<SimpleVector&>(*actmin_eq_ubndd->vec);
      SimpleVector& actmax_eq_root_ubndd = dynamic_cast<SimpleVector&>(*actmax_eq_ubndd->vec);
      SimpleVector& actmin_ineq_root_ubndd = dynamic_cast<SimpleVector&>(*actmin_ineq_ubndd->vec);
      SimpleVector& actmax_ineq_root_ubndd = dynamic_cast<SimpleVector&>(*actmax_ineq_ubndd->vec);

      addActivityOfBlock(mat_A.Bmat->getStorageDynamicRef(), actmin_eq_root_part, actmin_eq_root_ubndd, actmax_eq_root_part, actmax_eq_root_ubndd,
            xlow_root, ixlow_root, xupp_root, ixupp_root);

      addActivityOfBlock(mat_C.Bmat->getStorageDynamicRef(), actmin_ineq_root_part, actmin_ineq_root_ubndd, actmax_ineq_root_part, actmax_ineq_root_ubndd,
            xlow_root, ixlow_root, xupp_root, ixupp_root);
   }

   SimpleVector& actmin_eq_link_part = dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl);
   SimpleVector& actmax_eq_link_part = dynamic_cast<SimpleVector&>(*actmax_eq_part->vecl);
   SimpleVector& actmin_ineq_link_part = dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl);
   SimpleVector& actmax_ineq_link_part = dynamic_cast<SimpleVector&>(*actmax_ineq_part->vecl);

   SimpleVector& actmin_eq_link_ubndd = dynamic_cast<SimpleVector&>(*actmin_eq_ubndd->vecl);
   SimpleVector& actmax_eq_link_ubndd = dynamic_cast<SimpleVector&>(*actmax_eq_ubndd->vecl);
   SimpleVector& actmin_ineq_link_ubndd = dynamic_cast<SimpleVector&>(*actmin_ineq_ubndd->vecl);
   SimpleVector& actmax_ineq_link_ubndd = dynamic_cast<SimpleVector&>(*actmax_ineq_ubndd->vecl);

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

            SimpleVector& actmin_eq_child_ubndd = dynamic_cast<SimpleVector&>(*actmin_eq_ubndd->children[node]->vec);
            SimpleVector& actmax_eq_child_ubndd = dynamic_cast<SimpleVector&>(*actmax_eq_ubndd->children[node]->vec);

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

            SimpleVector& actmin_ineq_child_ubndd = dynamic_cast<SimpleVector&>(*actmin_ineq_ubndd->children[node]->vec);
            SimpleVector& actmax_ineq_child_ubndd = dynamic_cast<SimpleVector&>(*actmax_ineq_ubndd->children[node]->vec);

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
      MPI_Allreduce(MPI_IN_PLACE, dynamic_cast<SimpleVector*>(actmin_eq_part->vecl)->elements(), actmin_eq_part->vecl->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, dynamic_cast<SimpleVector*>(actmax_eq_part->vecl)->elements(), actmax_eq_part->vecl->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, dynamic_cast<SimpleVector*>(actmin_ineq_part->vecl)->elements(), actmin_ineq_part->vecl->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, dynamic_cast<SimpleVector*>(actmax_ineq_part->vecl)->elements(), actmax_ineq_part->vecl->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      MPI_Allreduce(MPI_IN_PLACE, dynamic_cast<SimpleVector*>(actmin_eq_ubndd->vecl)->elements(), actmin_eq_ubndd->vecl->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, dynamic_cast<SimpleVector*>(actmax_eq_ubndd->vecl)->elements(), actmax_eq_ubndd->vecl->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, dynamic_cast<SimpleVector*>(actmin_ineq_ubndd->vecl)->elements(), actmin_ineq_ubndd->vecl->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, dynamic_cast<SimpleVector*>(actmax_ineq_ubndd->vecl)->elements(), actmax_ineq_ubndd->vecl->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }

   /* set activities to infinity // theoretically not necessary but for debugging */
   for(int row = 0; row < actmin_eq_part->vecl->n; ++row)
   {
      if( dynamic_cast<SimpleVector&>(*actmin_eq_ubndd->vecl)[row] >= 2 )
         dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl)[row] = -std::numeric_limits<double>::infinity();
      if( dynamic_cast<SimpleVector&>(*actmax_eq_ubndd->vecl)[row] >= 2 )
         dynamic_cast<SimpleVector&>(*actmax_eq_part->vecl)[row] = std::numeric_limits<double>::infinity();
   }
   for(int row = 0; row < actmin_ineq_part->vecl->n; ++row)
   {
      if( dynamic_cast<SimpleVector&>(*actmin_ineq_ubndd->vecl)[row] >= 2 )
         dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl)[row] = -std::numeric_limits<double>::infinity();
      if( dynamic_cast<SimpleVector&>(*actmax_ineq_ubndd->vecl)[row] >= 2 )
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

#ifdef TRACK_ROW
   if(!nodeIsDummy(ROW_NODE, ROW_SYS))
   {
   double act_min = (ROW_SYS == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmin_eq_part, ROW_NODE, ROW_BLOCK)[ROW] :
         getSimpleVecRowFromStochVec(*actmin_ineq_part, ROW_NODE, ROW_BLOCK)[ROW];
   double act_max = (ROW_SYS == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmax_eq_part, ROW_NODE, ROW_BLOCK)[ROW] :
         getSimpleVecRowFromStochVec(*actmax_ineq_part, ROW_NODE, ROW_BLOCK)[ROW];

   double act_min_ubndd = (ROW_SYS == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmin_eq_ubndd, ROW_NODE, ROW_BLOCK)[ROW] :
         getSimpleVecRowFromStochVec(*actmin_ineq_ubndd, ROW_NODE, ROW_BLOCK)[ROW];
   double act_max_ubndd = (ROW_SYS == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmax_eq_ubndd, ROW_NODE, ROW_BLOCK)[ROW] :
         getSimpleVecRowFromStochVec(*actmax_ineq_ubndd, ROW_NODE, ROW_BLOCK)[ROW];


   std::cout << "TRACKING: Recomputed activity of row " << ROW << std::endl;
   std::cout << "\tnew min/max activity is: " << act_min << "/" << act_max << ", min/max unbounded counters are " << act_min_ubndd << "/" << act_max_ubndd << std::endl;
   }
#endif
}


/** Computes minimal and maximal activity of all rows in given matrix. Adds activities to min/max_activities accordingly. */
void PresolveData::addActivityOfBlock( const SparseStorageDynamic& matrix, SimpleVector& min_partact, SimpleVector& unbounded_min, SimpleVector& max_partact,
      SimpleVector& unbounded_max, const SimpleVector& xlow, const SimpleVector& ixlow, const SimpleVector& xupp, const SimpleVector& ixupp) const
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
            if( ixlow[col] != 0.0)
               min_partact[row] += entry * xlow[col];
            else
               ++unbounded_min[row];
            if( ixupp[col] != 0.0 )
               max_partact[row] += entry * xupp[col];
            else
               ++unbounded_max[row];
         }
         else
         {
            if( ixupp[col] != 0.0 )
               min_partact[row] += entry * xupp[col];
            else
               ++unbounded_min[row];
            if( ixlow[col] != 0.0 )
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
   MPI_Allreduce(MPI_IN_PLACE, &outdated_linking_var_bounds, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);

   if(!outdated_linking_var_bounds)
      return;

   if(distributed)
   {

      SimpleVector& xlow = getSimpleVecColFromStochVec(*presProb->blx, -1);
      SimpleVector& xupp = getSimpleVecColFromStochVec(*presProb->bux, -1);
      SimpleVector& ixlow = getSimpleVecColFromStochVec(*presProb->ixlow, -1);
      SimpleVector& ixupp = getSimpleVecColFromStochVec(*presProb->ixupp, -1);

      /* copy old values for later compairson */
      SimpleVector* xlow_old = xlow.cloneFull();
      SimpleVector* xupp_old = xupp.cloneFull();
      SimpleVector* ixlow_old =  xlow.cloneFull();
      SimpleVector* ixupp_old =  xupp.cloneFull();

      MPI_Allreduce(MPI_IN_PLACE, xlow.elements(), xlow.length(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
      MPI_Allreduce(MPI_IN_PLACE, ixlow.elements(), ixlow.length(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
      MPI_Allreduce(MPI_IN_PLACE, xupp.elements(), xupp.length(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
      MPI_Allreduce(MPI_IN_PLACE, ixupp.elements(), ixupp.length(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );


      // this will affect the activities of basically all rows - use with care
      for(int col = 0; col < xlow.length(); ++col)
      {
         double xu_old = std::numeric_limits<double>::infinity();
         double xl_old = -std::numeric_limits<double>::infinity();
         double xu = std::numeric_limits<double>::infinity();
         double xl = -std::numeric_limits<double>::infinity();

         if( (*ixupp_old)[col] == 1.0 )
            xu_old = (*xupp_old)[col];
         if( (*ixlow_old)[col] == 1.0 )
            xl_old = (*xlow_old)[col];
         if( ixupp[col] == 1.0 )
            xu = xupp[col];
         if( ixlow[col] == 1.0 )
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
   MPI_Allreduce(MPI_IN_PLACE, &outdated_activities, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);

   // todo : better criterion - if this and that many changes etc
   MPI_Allreduce(MPI_IN_PLACE, &linking_rows_need_act_computation, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   if(!outdated_activities && linking_rows_need_act_computation == 0)
      return;

   /* allreduce unbounded changes, compute rows that now have at most one unbounded entry, strore the local
    * rows in the changes array, allredues MPI_SUM the changes array and update local activities with the global ones
    */
   if(distributed)
   {
      MPI_Allreduce(MPI_IN_PLACE, array_act_unbounded_chgs, lenght_array_act_chgs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }

   dynamic_cast<SimpleVector&>(*actmin_eq_ubndd->vecl).axpy( 1.0, *actmin_eq_ubndd_chgs);
   dynamic_cast<SimpleVector&>(*actmax_eq_ubndd->vecl).axpy( 1.0, *actmax_eq_ubndd_chgs);
   dynamic_cast<SimpleVector&>(*actmin_ineq_ubndd->vecl).axpy( 1.0, *actmin_ineq_ubndd_chgs);
   dynamic_cast<SimpleVector&>(*actmax_ineq_ubndd->vecl).axpy( 1.0, *actmax_ineq_ubndd_chgs);

   /* equality system */
   for(int row = 0; row < actmin_eq_ubndd->vecl->n; ++row)
   {
      if( dynamic_cast<SimpleVector&>(*actmin_eq_ubndd->vecl)[row] < 2
            && dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl)[row] == -std::numeric_limits<double>::infinity())
      {
         (*actmin_eq_chgs)[row] = computeLocalLinkingRowMinOrMaxActivity(EQUALITY_SYSTEM, row, false);
         dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl)[row] = 0;
      }

      if( dynamic_cast<SimpleVector&>(*actmax_eq_ubndd->vecl)[row] < 2
            && dynamic_cast<SimpleVector&>(*actmax_eq_part->vecl)[row] == std::numeric_limits<double>::infinity())
      {
         (*actmax_eq_chgs)[row] = computeLocalLinkingRowMinOrMaxActivity(EQUALITY_SYSTEM, row, true);
         dynamic_cast<SimpleVector&>(*actmax_eq_part->vecl)[row] = 0;
      }
   }

   /* inequality system */
   for(int row = 0; row < actmin_ineq_ubndd->vecl->n; ++row)
   {
      if( dynamic_cast<SimpleVector&>(*actmin_ineq_ubndd->vecl)[row] < 2
            && dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl)[row] == -std::numeric_limits<double>::infinity())
      {
         (*actmin_ineq_chgs)[row] = computeLocalLinkingRowMinOrMaxActivity(INEQUALITY_SYSTEM, row, false);
         dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl)[row] = 0;
      }

      if( dynamic_cast<SimpleVector&>(*actmax_ineq_ubndd->vecl)[row] < 2
            && dynamic_cast<SimpleVector&>(*actmax_ineq_part->vecl)[row] == std::numeric_limits<double>::infinity())
      {
         (*actmax_ineq_chgs)[row] = computeLocalLinkingRowMinOrMaxActivity(INEQUALITY_SYSTEM, row, true);
         dynamic_cast<SimpleVector&>(*actmax_ineq_part->vecl)[row] = 0;
      }
   }

   if( distributed )
      MPI_Allreduce(MPI_IN_PLACE, array_act_chgs, lenght_array_act_chgs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

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
      if( dynamic_cast<SimpleVector&>(*actmin_eq_ubndd->vecl)[i] < 2 )
         assert( std::fabs(dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl)[i]) != std::numeric_limits<double>::infinity() );
      else
         assert( dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl)[i] == -std::numeric_limits<double>::infinity() );
   }

   for(int i = 0; i < actmin_ineq_ubndd->vecl->n; ++i)
   {
      if( dynamic_cast<SimpleVector&>(*actmin_ineq_ubndd->vecl)[i] < 2 )
         assert( std::fabs(dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl)[i]) != std::numeric_limits<double>::infinity() );
      else
         assert( dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl)[i] == -std::numeric_limits<double>::infinity() );
   }
#endif
}

/** allreduce changes in the nnz counters and apply them locally */
void PresolveData::allreduceAndApplyNnzChanges()
{
   MPI_Allreduce(MPI_IN_PLACE, &outdated_nnzs, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);

   MPI_Barrier(MPI_COMM_WORLD);
   if(!outdated_nnzs)
      return;
   MPI_Barrier(MPI_COMM_WORLD);

   if( distributed )
      MPI_Allreduce(MPI_IN_PLACE, array_nnz_chgs, length_array_nnz_chgs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   /* update local nnzCounters */
   nnzs_col->vec->axpy(1.0, *nnzs_col_chgs);
   nnzs_row_A->vecl->axpy(1.0, *nnzs_row_A_chgs);
   nnzs_row_C->vecl->axpy(1.0, *nnzs_row_C_chgs);

   // todo : this can be done more efficiently, e.g. while substracting
   // todo : this has still flaws - new singleton rows in B0 and D0 are not communicated properly for some reason - might happen else where
   for( int i = 0; i < nnzs_col_chgs->length(); ++i )
   {
      if( (*nnzs_col_chgs)[i] > 0.0 && dynamic_cast<SimpleVector&>(*nnzs_col->vec)[i] == 1 )
         singleton_cols.push_back( sCOLINDEX(-1, i) );
   }
   
   for( int i = 0; i < nnzs_row_A_chgs->length(); ++i )
   {
      if( (*nnzs_row_A_chgs)[i] > 0.0 && dynamic_cast<SimpleVector&>(*nnzs_row_A->vec)[i] == 1 )
         singleton_rows.push_back( sROWINDEX(EQUALITY_SYSTEM, -2, i) );
   }
   
   for( int i = 0; i < nnzs_row_C_chgs->length(); ++i )
   {
      if( (*nnzs_row_C_chgs)[i] > 0.0 && dynamic_cast<SimpleVector&>(*nnzs_row_C->vec)[i] == 1 )
         singleton_rows.push_back( sROWINDEX(INEQUALITY_SYSTEM, -2, i) );
   }

#ifndef NDEBUG
   double minval = -1.0;
   int index = -1;
   nnzs_col->min(minval, index);
   assert( minval >= 0.0 );
   nnzs_row_A->vecl->min(minval, index);
   assert(minval >= 0.0);
   nnzs_row_C->vecl->min(minval, index);
   assert(minval >= 0.0);
#endif

   nnzs_col_chgs->setToZero();
   nnzs_row_A_chgs->setToZero();
   nnzs_row_C_chgs->setToZero();

   outdated_nnzs = false;
}

void PresolveData::allreduceAndApplyBoundChanges()
{
   MPI_Allreduce(MPI_IN_PLACE, &outdated_lhsrhs, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);

   if(!outdated_lhsrhs)
      return;

   if(distributed)
   {
      MPI_Allreduce(MPI_IN_PLACE, array_bound_chgs, lenght_array_bound_chgs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
   MPI_Allreduce(MPI_IN_PLACE, &outdated_obj_vector, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);

   if(!outdated_obj_vector)
      return;

   if(distributed)
      MPI_Allreduce(MPI_IN_PLACE, objective_vec_chgs->elements(), objective_vec_chgs->length(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*presProb->g).vec)->axpy(1.0, *objective_vec_chgs);

   objective_vec_chgs->setToZero();
   outdated_obj_vector = false;
}

void PresolveData::allreduceObjOffset()
{
   if(distributed)
      MPI_Allreduce(MPI_IN_PLACE, &obj_offset_chgs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   objOffset += obj_offset_chgs;
   obj_offset_chgs = 0;
}

void PresolveData::initNnzCounter(StochVector& nnzs_row_A, StochVector& nnzs_row_C, StochVector& nnzs_col) const
{
   StochGenMatrix& A = dynamic_cast<StochGenMatrix&>(*(presProb->A));
   StochGenMatrix& C = dynamic_cast<StochGenMatrix&>(*(presProb->C));

   StochVectorHandle colClone(dynamic_cast<StochVector*>(nnzs_col.clone()));

   A.getNnzPerRow(nnzs_row_A);
   C.getNnzPerRow(nnzs_row_C);
   A.getNnzPerCol(nnzs_col);
   C.getNnzPerCol(*colClone);

   nnzs_col.axpy(1.0, *colClone);
}

void PresolveData::initSingletons()
{
   /* rows of A */
   /* B0 */
   for(int i = 0; i < nnzs_row_A->vec->n; ++i)
   if( dynamic_cast<SimpleVector&>(*nnzs_row_A->vec)[i] == 1.0)
   {
      singleton_rows.push_back( sROWINDEX( EQUALITY_SYSTEM, -1, i ));
   }

   /* Bl0 */
   if(nnzs_row_A->vecl != NULL)
   {
      for(int i = 0; i < nnzs_row_A->vecl->n; ++i)
      {
         if( dynamic_cast<SimpleVector&>(*nnzs_row_A->vecl)[i] == 1.0)
            singleton_rows.push_back( sROWINDEX( EQUALITY_SYSTEM, -2, i ));
      }
   }

   /* children An + Bn */
   for(int i = 0; i < nChildren; ++i)
   {
      if( !nodeIsDummy(i, EQUALITY_SYSTEM) )
      {
         for(int j = 0; j < nnzs_row_A->children[i]->vec->n; ++j)
         {
            if( dynamic_cast<SimpleVector&>(*nnzs_row_A->children[i]->vec)[j] == 1.0)
               singleton_rows.push_back( sROWINDEX( EQUALITY_SYSTEM, i, j ));
         }
      }
   }

   /* rows of C */
   /* B0 */
   for(int i = 0; i < nnzs_row_C->vec->n; ++i)
      if( dynamic_cast<SimpleVector&>(*nnzs_row_C->vec)[i] == 1.0)
         singleton_rows.push_back( sROWINDEX( INEQUALITY_SYSTEM, -1, i ));

   /* Bl0 */
   if( nnzs_row_C->vecl != NULL)
   {
      for( int i = 0; i < nnzs_row_C->vecl->n; ++i )
      {
         if( dynamic_cast<SimpleVector&>(*nnzs_row_C->vecl)[i] == 1.0 )
            singleton_rows.push_back(sROWINDEX( INEQUALITY_SYSTEM, -2, i ));
      }
   }

   /* children An + Bn */
   for( int i = 0; i < nChildren; ++i )
   {
      if( !nodeIsDummy(i, INEQUALITY_SYSTEM) )
      {
         for( int j = 0; j < nnzs_row_C->children[i]->vec->n; ++j )
         {
            if( dynamic_cast<SimpleVector&>(*nnzs_row_C->children[i]->vec)[j] == 1.0 )
               singleton_rows.push_back(sROWINDEX( INEQUALITY_SYSTEM, i, j ));
         }
      }
   }

   /* columns */
   for(int i = 0; i < nnzs_col->vec->n; ++i)
   {
      if( dynamic_cast<SimpleVector&>(*nnzs_col->vec)[i] == 1.0)
         singleton_cols.push_back( sCOLINDEX(-1, i));
   }

   for(int i = 0; i < nChildren; ++i)
   {
      if(nnzs_col->children[i]->vec != NULL)
      {
         for(int j = 0; j < nnzs_col->children[i]->vec->n; ++j)
         {
            if( dynamic_cast<SimpleVector&>(*nnzs_col->children[i]->vec)[j] == 1.0)
               singleton_cols.push_back( sCOLINDEX(i, j));
         }
      }
   }

}

bool PresolveData::reductionsEmpty()
{
   double recv = 0;
   if(distributed)
   {
      MPI_Allreduce(MPI_IN_PLACE, &outdated_activities, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_lhsrhs, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_nnzs, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_linking_var_bounds, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_obj_vector, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);

      MPI_Allreduce(&obj_offset_chgs, &recv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_objvector, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD );
   }
   return !outdated_obj_vector && !outdated_activities && !outdated_lhsrhs && !outdated_linking_var_bounds && !outdated_nnzs && (recv == 0);
}

// todo : if small entry was removed from system no postsolve is necessary - if coefficient was removed because impact of changes in variable are small
// also rhs lhs will be adjusted - this has to be reversed later - there will also be a problem with the reduced costs in than particular row ?
// todo : postsolve if bounds adjusted because of deleted matrix entry simply reverse the adjustment - no changes in multipliers
//- row stays active / inactive ? not true..
void PresolveData::deleteEntry(SystemType system_type, int node, BlockType block_type, int row,
      int& col_idx, int& row_end)
{
   assert(-1 <= node && node < nChildren);

   SparseStorageDynamic& storage = getSparseGenMatrix(system_type, node , block_type)->getStorageDynamicRef();
   
   double val = storage.getMat(col_idx);
   int col = storage.getJcolM(col_idx);

   storage.removeEntryAtIndex(row, col_idx);

   --col_idx;
   --row_end;
   ++elements_deleted;

   SimpleVector& xlower = getSimpleVecColFromStochVec(*presProb->blx, node);

   if(block_type == LINKING_CONS_BLOCK || block_type == LINKING_VARS_BLOCK || node == -1)
      outdated_nnzs = true;

   /* adjust rhs and lhs */
   adjustMatrixRhsLhsBy(system_type, node, block_type, row, -val * xlower[col]);

   /* adjust activity */
   // todo : necessary ? impact should be low

   /* adjust nnz counters */
   removeIndexRow(system_type, node, block_type, row, 1);
   removeIndexColumn(node, block_type, col, 1);

}

void PresolveData::resetOriginallyFreeVarsBounds(const sData& orig_prob)
{

#ifndef NDEBUG
   if(my_rank == 0)
      std::cout << "Resetting all presolved variable bounds of originally free variables" <<::endl; 
#endif
   unsigned long long n = 0;
   for( int node = -1; node < nChildren; ++node )
   {
      n += resetOriginallyFreeVarsBounds( getSimpleVecColFromStochVec(*orig_prob.ixlow, node), getSimpleVecColFromStochVec(*orig_prob.ixupp, node), node);
      if(my_rank != 0 && node == -1)
      {
         n = 0;
      }
   }

#ifndef NDEBUG
   MPI_Allreduce(MPI_IN_PLACE, &n, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

   if(my_rank == 0)
      std::cout << "Reset " << n << " bounds" << std::endl;
#endif
}

long PresolveData::resetOriginallyFreeVarsBounds(const SimpleVector& ixlow_orig, const SimpleVector& ixupp_orig, int node)
{
   long reset_bounds = 0;

   if( nodeIsDummy( node, EQUALITY_SYSTEM ) && nodeIsDummy( node, INEQUALITY_SYSTEM ) )
      return reset_bounds;

   SimpleVector& ixlow = getSimpleVecColFromStochVec(*presProb->ixlow, node);
   SimpleVector& ixupp = getSimpleVecColFromStochVec(*presProb->ixupp, node);

   SimpleVector& xlow = getSimpleVecColFromStochVec(*presProb->blx, node);
   SimpleVector& xupp = getSimpleVecColFromStochVec(*presProb->bux, node);

   SimpleVector& nnzs_col_vec = getSimpleVecColFromStochVec(*nnzs_col, node);

   /* check whether row that implied bound is still there - if so, and if the variable is still in that row we can remove the bound again
    * since it should still be implied 
    */

   /* todo: actually need to check whether a bound is still an implied one - if so we can reset it - this needs more mechanisms */
   /* store row that implied bound - if row still there - check if bound still implied (or even better bound implied) - if so - reset bound */
   /* print stats on how many bounds were reset */
   for(int col = 0; col < ixlow.n; ++col)
   {
      /* do not reset fixed columns */
      if( nnzs_col_vec[col] != 0 && (ixupp_orig[col] == 0.0 || ixlow_orig[col] == 0.0) )
      {
         int sys_row_lower = getSimpleVecColFromStochVec(*lower_bound_implied_by_system, node)[col];
         int row_lower = getSimpleVecColFromStochVec(*lower_bound_implied_by_row, node)[col];
         int node_row_lower = getSimpleVecColFromStochVec(*lower_bound_implied_by_node, node)[col];
         int sys_row_upper = getSimpleVecColFromStochVec(*upper_bound_implied_by_system, node)[col];
         int row_upper = getSimpleVecColFromStochVec(*upper_bound_implied_by_row, node)[col];
         int node_row_upper = getSimpleVecColFromStochVec(*upper_bound_implied_by_node, node)[col];

         /* do not reset bounds implied by singleton rows since these rows are already removed from the problem */
         if( ixupp_orig[col] == 0.0 && ixupp[col] == 1.0 )
         {
            if( !(0 <= sys_row_upper && sys_row_upper <= 1) )
               std::cout << sys_row_upper << std::endl;
            assert( 0 <= sys_row_upper && sys_row_upper <= 1);
            assert(row_upper != -10);
            assert(node_row_upper != -10);

            const StochVector& nnzs = (sys_row_upper == 0) ? *nnzs_row_A : *nnzs_row_C;

            if( getSimpleVecRowFromStochVec(nnzs, (node_row_upper == -2) ? -1 : node_row_upper, 
               (node_row_upper == -2) ? LINKING_CONS_BLOCK : LINKING_VARS_BLOCK)[row_upper] != 0.0 )
            {
               ixupp[col] = 0;
               xupp[col] = 0.0;
               ++reset_bounds;
            }
         }

         if( ixlow_orig[col] == 0.0 && ixlow[col] == 1.0 )
         {
            assert( 0 <= sys_row_lower && sys_row_lower <= 1);
            assert(row_lower != -10);
            assert(node_row_lower != -10);

            const StochVector& nnzs = (sys_row_lower == 0) ? *nnzs_row_A : *nnzs_row_C;

            // todo : problably should also check whether variable is still in row - not necessary for so far implemented presolvers
            if( getSimpleVecRowFromStochVec(nnzs, (node_row_lower == -2) ? -1 : node_row_lower, 
               (row_lower == -2) ? LINKING_CONS_BLOCK : LINKING_VARS_BLOCK)[row_lower] != 0.0 )
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

bool PresolveData::varBoundImpliedFreeBy( bool upper, int node_col, int col, SystemType system_type, int node_row, int row, bool linking_row )
{
   node_row = ( !linking_row) ? node_row : -2;

   if( upper )
   {
      return (getSimpleVecColFromStochVec( *upper_bound_implied_by_system, node_col)[col] == system_type && 
         getSimpleVecColFromStochVec( *upper_bound_implied_by_node, node_col)[col] == node_row &&
         getSimpleVecColFromStochVec( *upper_bound_implied_by_row, node_col)[col] == row);
   }
   else
   {
      return (getSimpleVecColFromStochVec( *lower_bound_implied_by_system, node_col)[col] == system_type && 
         getSimpleVecColFromStochVec( *lower_bound_implied_by_node, node_col)[col] == node_row &&
         getSimpleVecColFromStochVec( *lower_bound_implied_by_row, node_col)[col] == row);
   }

}

void PresolveData::fixEmptyColumn(int node, int col, double val)
{
   assert(-1 <= node && node < nChildren);
   assert(0 <= col);
   //todo
   postsolver->notifyFixedEmptyColumn(node, col, val);

   removeColumn(node, col, val);

   if( node != -1)
      assert( getSimpleVecColFromStochVec(*nnzs_col, node)[col] == 0.0 );
}

void PresolveData::fixColumn(int node, int col, double value)
{
   assert(-1 <= node && node < nChildren);
   assert(0 <= col);

   /* current upper and lower bound as well als column - if linking variable then only proc zero stores current root column */
   std::vector<int> idx;
   std::vector<double> val;
   
   // todo
   //buildColForPostsolve( node, col, idx, val);
   postsolver->notifyFixedColumn(node, col, value, idx, val);

#ifndef NDEBUG
   double ixlow = getSimpleVecColFromStochVec(*presProb->ixlow, node)[col];
   double ixupp = getSimpleVecColFromStochVec(*presProb->ixupp, node)[col];
   double xlow = (ixupp == 1.0) ? getSimpleVecColFromStochVec(*presProb->blx, node)[col] : std::numeric_limits<double>::infinity();
   double xupp = (ixlow == 1.0) ? getSimpleVecColFromStochVec(*presProb->bux, node)[col] : std::numeric_limits<double>::infinity();
   assert(ixlow == 1.0);
   assert(ixupp == 1.0);
   assert(PIPSisEQ(xlow, xupp, 1e-10));
   assert(PIPSisEQ(xlow, value, 1e-10));
#endif

   removeColumn(node, col, value);

   if( node != -1)
      assert( getSimpleVecColFromStochVec(*nnzs_col, node)[col] == 0.0 );
}

bool PresolveData::rowPropagatedBounds( SystemType system_type, int node_row, BlockType block_type, int row, int col, double ubx, double lbx)
{
   assert( -1 <= node_row && node_row < nChildren );

   //const double numerical_upper_threshold = std::numeric_limits<double>::max();

   const int node_var = (block_type == LINKING_VARS_BLOCK) ? -1 : node_row;

   assert( 0 <= col && col < getSimpleVecColFromStochVec( *presProb->ixlow, node_var ).n );

   /* check for infeasibility of the newly found bounds */
   const double ixlow = getSimpleVecColFromStochVec( *presProb->ixlow, node_var )[col];
   const double xlow = getSimpleVecColFromStochVec( *presProb->blx, node_var )[col];
   const double ixupp = getSimpleVecColFromStochVec( *presProb->ixupp, node_var )[col];
   const double xupp = getSimpleVecColFromStochVec( *presProb->bux, node_var )[col];

#ifdef TRACK_COLUMN
   if( NODE == node_var && COLUMN == col && (my_rank == 0 || node_var != -1) && !nodeIsDummy(NODE, EQUALITY_SYSTEM) )
   {
      std::cout << "TRACKING COLUMN: new bounds [" << lbx << ", " << ubx << "] propagated for column " << COLUMN << " from row " << row << " node " << node_row << " in " <<
            ( (system_type == EQUALITY_SYSTEM) ? "EQU_SYS" : "INEQ_SYS") << ":" << std::endl;
      writeRowLocalToStreamDense(std::cout, system_type, node_row, block_type, row);
      std::cout << "\tbounds were [" << xlow<< ", " << xupp<< "]" << std::endl;
   }
#endif

   if( ( ixlow != 0.0 && PIPSisLT(ubx, xlow) )
         || ( ixupp != 0.0 && PIPSisLT(xupp, lbx) )
         || PIPSisLT(ubx, lbx) )
   {
      std::cout << "[" << lbx << ", " << ubx << "] not in [" << ((ixlow != 0.0) ? xlow : -std::numeric_limits<double>::infinity()) << ", " << 
         ((ixupp != 0.0) ? xupp : std::numeric_limits<double>::infinity()) << "]" << std::endl;
      abortInfeasible(MPI_COMM_WORLD, "Row Propagation detected infeasible new bounds!", "PresolveData.C", "rowPropagatedBounds");
   }

   /* adjust bounds */
   bool bounds_changed = false;

   // we do not tighten bounds if impact is too low or bound is bigger than 10e8 // todo : maybe different limit
   // set lower bound
   // if( fabs(lbx) < 1e8 && (ixlow== 0.0  || feastol * 1e3 <= fabs(xlow - lbx) ) )
   if( ubx < std::numeric_limits<double>::infinity() && ( ixupp == 0.0 || PIPSisLT(ubx, xupp) ) )
   {
#ifdef TRACK_COLUMN
      if( NODE == node_var && COLUMN == col && (my_rank == 0 || node_var != -1) && !nodeIsDummy(NODE, EQUALITY_SYSTEM) )
         std::cout << "TRACKING COLUMN: new upper bound through propagation" << std::endl;
#endif
      if( updateUpperBoundVariable(node_var, col, ubx) )
      {
         /* store node and row that implied the bound (necessary for resetting bounds later on) */
         getSimpleVecColFromStochVec(*upper_bound_implied_by_system, node_var)[col] = system_type;
         getSimpleVecColFromStochVec(*upper_bound_implied_by_row, node_var)[col] = row;
         getSimpleVecColFromStochVec(*upper_bound_implied_by_node, node_var)[col] = (block_type == LINKING_CONS_BLOCK) ? -2 : node_row; // -2 for linking rows

         bounds_changed = true;
      }
   }
  // if( fabs(ubx) < 1e8 && (ixupp== 0.0  || feastol * 1e3 <= fabs(xupp- ubx) ) )
   if( - std::numeric_limits<double>::infinity() < lbx && ( ixlow == 0.0 || PIPSisLT(xlow, lbx)) )
   {
#ifdef TRACK_COLUMN
      if( NODE == node_var && COLUMN == col && (my_rank == 0 || node_var != -1) && !nodeIsDummy(NODE, EQUALITY_SYSTEM) )
         std::cout << "TRACKING COLUMN: new lower bound through propagation" << std::endl;
#endif
      if( updateLowerBoundVariable(node_var, col, lbx) )
      {
         /* store node and row that implied the bound (necessary for resetting bounds later on) */
         getSimpleVecColFromStochVec(*lower_bound_implied_by_system, node_var)[col] = system_type;
         getSimpleVecColFromStochVec(*lower_bound_implied_by_row, node_var)[col] = row;
         getSimpleVecColFromStochVec(*lower_bound_implied_by_node, node_var)[col] = (block_type == LINKING_CONS_BLOCK) ? -2 : node_row; // -2 for linking rows

         bounds_changed = true;
      }
   }

   /// every process should have the same root node data thus all of them should propagate it's rows similarly
   if( bounds_changed && (node_row == -1 || block_type == LINKING_VARS_BLOCK) )
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
//   postsolver->notifyRowPropagated(system_type, node_row, row, (block_type == LINKING_CONS_BLOCK), col, lbx, ubx, mat->getStorageDynamic()->M + row_start,
//         mat->getStorageDynamic()->jcolM + row_start, row_end - row_start );

   return bounds_changed;
}

void PresolveData::tightenRowBoundsParallelRow(SystemType system_type, int node, int row, double lhs, double rhs, bool linking)
{
   assert( -1 <= node && node <= nChildren);
   assert( !linking );
   assert(system_type == INEQUALITY_SYSTEM);

#ifdef TRACK_ROW
/*   if(row == ROW && node == ROW_NODE && system_type == ROW_SYS && !nodeIsDummy(ROW_NODE, ROW_SYS))
   {
      std::cout << "TRACKING: RHS LHS of row " << ROW << " being adjusted by " << value << std::endl;
      writeRowLocalToStreamDense(std::cout, ROW_SYS, ROW_NODE, ROW_BLOCK, ROW);
   }*/
#endif

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
      //getSimpleVecRowFromStochVec(*presProb->bA, node, block_type)[row] += value;
   }
   else
   {
      if( lhs != -std::numeric_limits<double>::infinity() )
      {
         if( !PIPSisEQ(getSimpleVecRowFromStochVec(*presProb->iclow, node, CHILD_BLOCK)[row], 1.0) )
         {
            getSimpleVecRowFromStochVec(*presProb->iclow, node, CHILD_BLOCK)[row] = 1.0;
            getSimpleVecRowFromStochVec(*presProb->bl, node, CHILD_BLOCK)[row] = lhs;

         }
         else
            getSimpleVecRowFromStochVec(*presProb->bl, node, CHILD_BLOCK)[row] = 
               std::max(getSimpleVecRowFromStochVec(*presProb->bl, node, CHILD_BLOCK)[row], lhs);
      }
      if( rhs != std::numeric_limits<double>::infinity() )
      {
         if( !PIPSisEQ(getSimpleVecRowFromStochVec(*presProb->icupp, node, CHILD_BLOCK)[row], 1.0) )
         {
            getSimpleVecRowFromStochVec(*presProb->icupp, node, CHILD_BLOCK)[row] = 1.0;
         }
         else
         {
            getSimpleVecRowFromStochVec(*presProb->bu, node, CHILD_BLOCK)[row] = 
               std::min(getSimpleVecRowFromStochVec(*presProb->bu, node, CHILD_BLOCK)[row], rhs);
         }
      }
      // todo!
   }

#ifdef TRACK_ROW
/*   if(row == ROW && node == ROW_NODE && system_type == ROW_SYS && !nodeIsDummy(ROW_NODE, ROW_SYS))
   {
      std::cout << "TRACKING: after RHS LHS adjustment " << std::endl;
      writeRowLocalToStreamDense(std::cout, ROW_SYS, ROW_NODE, ROW_BLOCK, ROW);
   }*/
#endif
}

/** this methods does not call any postsolve procedures but simply changes the bounds (lhs, rhs) of either A or B by value */
void PresolveData::adjustMatrixRhsLhsBy(SystemType system_type, int node, BlockType block_type, int row, double value)
{
   assert( -1 <= node && node <= nChildren);

   if( PIPSisEQ(value, 0.0) )
      return;

#ifdef TRACK_ROW
   if(row == ROW && node == ROW_NODE && system_type == ROW_SYS && !nodeIsDummy(ROW_NODE, ROW_SYS))
   {
      std::cout << "TRACKING: RHS LHS of row " << ROW << " being adjusted by " << value << std::endl;
      writeRowLocalToStreamDense(std::cout, ROW_SYS, ROW_NODE, ROW_BLOCK, ROW);
   }
#endif

   if(block_type == LINKING_CONS_BLOCK)
   {
      outdated_lhsrhs = true;
      (system_type == EQUALITY_SYSTEM) ? (*bound_chgs_A)[row] += value : (*bound_chgs_C)[row] += value;
      return;
   }

   if(system_type == EQUALITY_SYSTEM)
   {
      getSimpleVecRowFromStochVec(*presProb->bA, node, block_type)[row] += value;
   }
   else
   {
      if( getSimpleVecRowFromStochVec(*presProb->icupp, node, block_type)[row] == 1.0)
         getSimpleVecRowFromStochVec(*presProb->bu, node, block_type)[row] += value;

      if( getSimpleVecRowFromStochVec(*presProb->iclow, node, block_type)[row] == 1.0 )
         getSimpleVecRowFromStochVec(*presProb->bl, node, block_type)[row] += value;
   }

#ifdef TRACK_ROW
   if(row == ROW && node == ROW_NODE && system_type == ROW_SYS && !nodeIsDummy(ROW_NODE, ROW_SYS))
   {
      std::cout << "TRACKING: after RHS LHS adjustment " << std::endl;
      writeRowLocalToStreamDense(std::cout, ROW_SYS, ROW_NODE, ROW_BLOCK, ROW);
   }
#endif
}


// todo : make a finish block_deletion ?
void PresolveData::updateTransposedSubmatrix( SystemType system_type, int node, BlockType block_type, std::vector<std::pair<int, int> >& elements)
{
   SparseStorageDynamic& transposed = getSparseGenMatrix(system_type, node, block_type)->getStorageDynamicTransposedRef();

   for( size_t i = 0; i < elements.size(); ++i )
   {
      std::pair<int, int> entry = elements.at(i);
      const int row_A = entry.first;
      const int row_At = entry.second;

      const int start = transposed.getRowPtr(row_At).start;
      const int end = transposed.getRowPtr(row_At).end;
      int col_At;

      for( col_At = start; col_At < end; col_At++ )
      {
         if( transposed.getJcolM(col_At) == row_A )
            break;
      }

      transposed.removeEntryAtIndex(row_At, col_At);

      ++elements_deleted_transposed;
   }

   assert(elements_deleted == elements_deleted_transposed);
   elements_deleted = 0;
   elements_deleted_transposed = 0;
}

void PresolveData::removeIndexRow(SystemType system_type, int node, BlockType block_type, int row_index, int amount)
{
   assert(-1 <= node && node <= nChildren);
   assert(0 <= amount);

   if(amount == 0)
      return;

   /* linking constraints get stored */
   if(block_type == LINKING_CONS_BLOCK)
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
         getSimpleVecRowFromStochVec(*nnzs_row_A, node, block_type)[row_index] -= amount;
         if( getSimpleVecRowFromStochVec(*nnzs_row_A, node, block_type)[row_index]  == 1)
            singleton_rows.push_back( sROWINDEX( system_type, node, row_index ) );
         assert( 0 <= getSimpleVecRowFromStochVec(*nnzs_row_A, node, block_type)[row_index] );
      }
      else
      {
         getSimpleVecRowFromStochVec(*nnzs_row_C, node, block_type)[row_index] -= amount;
         if( getSimpleVecRowFromStochVec(*nnzs_row_C, node, block_type)[row_index]  == 1)
            singleton_rows.push_back( sROWINDEX( system_type, node, row_index ) );
         assert( 0 <= getSimpleVecRowFromStochVec(*nnzs_row_C, node, block_type)[row_index] );
      }
   }
}

void PresolveData::removeIndexColumn(int node, BlockType block_type, int col_index, int amount)
{
   assert(-1 <= node && node <= nChildren);
   if(amount == 0)
      return;

   /* linking constraints get stored */
   if(node == -1 || block_type == LINKING_VARS_BLOCK)
   {
      if(my_rank == 0 || node != -1)
      {
         (*nnzs_col_chgs)[col_index] -= amount;
         outdated_nnzs = true;
      }
   }
   else
   {
      getSimpleVecColFromStochVec( *nnzs_col, node )[col_index] -= amount;
      if( getSimpleVecColFromStochVec( *nnzs_col, node )[col_index]  == 1)
         singleton_cols.push_back( sCOLINDEX( node, col_index ) );
      assert(0 <= getSimpleVecColFromStochVec( *nnzs_col, node )[col_index] );
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
      removeColumnFromMatrix(EQUALITY_SYSTEM, -1, LINKING_VARS_BLOCK, col, fixation);
      removeColumnFromMatrix(INEQUALITY_SYSTEM, -1, LINKING_VARS_BLOCK, col, fixation);

      if(hasLinking(EQUALITY_SYSTEM))
         removeColumnFromMatrix(EQUALITY_SYSTEM, -1, LINKING_CONS_BLOCK, col, fixation);
      if( hasLinking(INEQUALITY_SYSTEM) )
         removeColumnFromMatrix(INEQUALITY_SYSTEM, -1, LINKING_CONS_BLOCK, col, fixation);

      for(int i = 0; i < nChildren; ++i)
      {
         if(!nodeIsDummy(i, EQUALITY_SYSTEM))
            removeColumnFromMatrix(EQUALITY_SYSTEM, i, LINKING_VARS_BLOCK, col, fixation);
         if( !nodeIsDummy(i, INEQUALITY_SYSTEM) )
            removeColumnFromMatrix(INEQUALITY_SYSTEM, i, LINKING_VARS_BLOCK, col, fixation);
      }
   }
   else
   {
      removeColumnFromMatrix(EQUALITY_SYSTEM, node, CHILD_BLOCK, col, fixation);
      removeColumnFromMatrix(INEQUALITY_SYSTEM, node, CHILD_BLOCK, col, fixation);

      if(hasLinking(EQUALITY_SYSTEM))
         removeColumnFromMatrix(EQUALITY_SYSTEM, node, LINKING_CONS_BLOCK, col, fixation);
      if(hasLinking(INEQUALITY_SYSTEM))
         removeColumnFromMatrix(INEQUALITY_SYSTEM, node, LINKING_CONS_BLOCK, col, fixation);
   }

   /* adjust objective function */
   if(node != -1 || my_rank == 0)
   {
      double objective_factor = getSimpleVecColFromStochVec(*presProb->g, node)[col];
      obj_offset_chgs += objective_factor * fixation;

   }

   /* mark column as removed */
   getSimpleVecColFromStochVec(*presProb->g, node)[col] = 0.0;
   getSimpleVecColFromStochVec(*presProb->ixlow, node)[col] = 0.0;
   getSimpleVecColFromStochVec(*presProb->ixupp, node)[col] = 0.0;
   getSimpleVecColFromStochVec(*presProb->blx, node)[col] = 0.0;
   getSimpleVecColFromStochVec(*presProb->bux, node)[col] = 0.0;
}

/** remove column - adjust lhs, rhs and activity as well as nnz_counters */
void PresolveData::removeColumnFromMatrix(SystemType system_type, int node, BlockType block_type, int col, double fixation)
{
   SparseGenMatrix* mat = getSparseGenMatrix(system_type, node, block_type);

   SparseStorageDynamic& matrix = mat->getStorageDynamicRef();
   SparseStorageDynamic& matrix_transp = mat->getStorageDynamicTransposedRef();

   assert(0 <= col && col < matrix_transp.getM());

   /* remove all entries in column from the sparse storage dynamic */
   for( int j = matrix_transp.getRowPtr(col).start; j < matrix_transp.getRowPtr(col).end; j++ )
   {
      const int row = matrix_transp.getJcolM(j);
      const double coeff = matrix_transp.getMat(j);

      assert( !PIPSisEQ(0.0, coeff) );

#ifdef TRACK_ROW
   if(row == ROW && node == ROW_NODE && system_type == ROW_SYS && !nodeIsDummy(ROW_NODE, ROW_SYS))
   {
      std::cout << "TRACKING: fixation of column " << col << " in tracked row"<< std::endl;
      writeRowLocalToStreamDense(std::cout, ROW_SYS, ROW_NODE, ROW_BLOCK, ROW);
   }
#endif

      /* remove the entry, adjust activity and row counters and rhs/lhs */
      matrix.removeEntryAtRowCol(row, col);

      removeIndexRow(system_type, node, block_type, row, 1);

      adjustMatrixRhsLhsBy(system_type, node, block_type, row, - coeff * fixation);

      adjustRowActivityFromDeletion(system_type, node, block_type, row, col, coeff);

#ifdef TRACK_ROW
   if( row == ROW && node == ROW_NODE && system_type == ROW_SYS && !nodeIsDummy(ROW_NODE, ROW_SYS) )
   {
      std::cout << "TRACKING ROW: after removal" << std::endl;
      writeRowLocalToStreamDense(std::cout, ROW_SYS, ROW_NODE, ROW_BLOCK, ROW);

      double act_min = (ROW_SYS == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmin_eq_part, ROW_NODE, ROW_BLOCK)[ROW] :
            getSimpleVecRowFromStochVec(*actmin_ineq_part, ROW_NODE, ROW_BLOCK)[ROW];
      double act_max = (ROW_SYS == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmax_eq_part, ROW_NODE, ROW_BLOCK)[ROW] :
            getSimpleVecRowFromStochVec(*actmax_ineq_part, ROW_NODE, ROW_BLOCK)[ROW];

      double act_min_ubndd = (ROW_SYS == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmin_eq_ubndd, ROW_NODE, ROW_BLOCK)[ROW] :
            getSimpleVecRowFromStochVec(*actmin_ineq_ubndd, ROW_NODE, ROW_BLOCK)[ROW];
      double act_max_ubndd = (ROW_SYS == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmax_eq_ubndd, ROW_NODE, ROW_BLOCK)[ROW] :
            getSimpleVecRowFromStochVec(*actmax_ineq_ubndd, ROW_NODE, ROW_BLOCK)[ROW];

      std::cout << "TRACKING: New activity of row " << ROW << std::endl;
      std::cout << "\tnew min/max activity is: " << act_min << "/" << act_max << ", min/max unbounded counters are " << act_min_ubndd << "/" << act_max_ubndd << std::endl;
   }
#endif

   }

   /* adjust column counters */
   removeIndexColumn(node, block_type, col, matrix_transp.getRowPtr(col).end - matrix_transp.getRowPtr(col).start);

   /* update the transposed */
   matrix_transp.clearRow( col );

   // todo assert(transposed and normal matrix are in sync)
}

// todo
void PresolveData::removeParallelRow(SystemType system_type, int node, int row, bool linking)
{
#ifdef TRACK_ROW
   if(row == ROW && node == ROW_NODE && system_type == ROW_SYS && !nodeIsDummy(ROW_NODE, ROW_SYS))
   {
      std::cout << "TRACKING: removal of tracked row as parallel row" << std::endl;
   }
#endif

   throw std::runtime_error("Not yet implemented");
//   if(postsolver)
//      postsolver->notifyParallelRow()

   removeRow(system_type, node, row, linking);
}

/* a singleton variable is substituted out of the problem and then it's original row can be removed from the problem */
void PresolveData::substituteVariableParallelRows(SystemType system_type, int node, int var1, int row1, int node_var1, int var2, int row2, int node_var2,
   double scalar, double translation)
{
#ifdef TRACK_ROW
// todo
#endif
   
   postsolver->notifyParallelRowSubstitution(system_type, node, var1, row1, node_var1, var2, row2, node_var2, scalar, translation);

   // delete the equality constraint which contained var2 (the substituted variable)
   removeRedundantRow( system_type, node, row2, false);
   assert( PIPSisZero(getSimpleVecColFromStochVec(*nnzs_col, node_var2)[var2]) );
   
   const double obj_var2 = getSimpleVecColFromStochVec(*presProb->g, node_var2)[var2];
   const double val_offset = translation * obj_var2;
   const double change_obj_var1 = scalar * obj_var2;

   removeColumn( node_var2, var2, 0.0 );

   if( node_var1 != -1 )
   {
      getSimpleVecColFromStochVec(*presProb->g, node)[var1] += change_obj_var1;
      obj_offset_chgs += val_offset;
   }
   else if( node == -1 )  
   {
      /* parallel rows in parent block - all processes should have detected this */
      getSimpleVecColFromStochVec(*presProb->g, -1)[var1] += change_obj_var1;

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
      std::vector<int> idx_row;
      std::vector<double> val_row;
      //BlockType block_type = (linking) ? LINKING_CONS_BLOCK : CHILD_BLOCK;
      // todo
      //buildRowForPostsolve( system_type, node, block_type, row, idx_row, val_row);

      postsolver->notifyRedundantRow(system_type, node, row, linking, idx_row, val_row);
   }
 
#ifdef TRACK_ROW
   if(row == ROW && node == ROW_NODE && system_type == ROW_SYS && !nodeIsDummy(ROW_NODE, ROW_SYS))
   {
      std::cout << "TRACKING: removal of tracked row as redundant row" << std::endl;
   }
#endif

   removeRow(system_type, node, row, linking);
}

void PresolveData::removeImpliedFreeColumnSingleton( SystemType system_type, int node_row, int row, bool linking_row, int node_col, int col )
{
   assert( system_type == EQUALITY_SYSTEM );

   const double rhs = getSimpleVecRowFromStochVec( *presProb->bA, node_row, (linking_row) ? LINKING_CONS_BLOCK : CHILD_BLOCK )[row];

   postsolver->notifyFreeColumnSingleton( system_type, node_row, row, linking_row, rhs, node_col, col, 
      dynamic_cast<const StochGenMatrix&>( (system_type == EQUALITY_SYSTEM) ? *presProb->A : *presProb->C ) );

#ifdef TRACK_COLUMN
  if( NODE == node_col && COLUMN == col && (my_rank == 0 || node_col != -1) && !nodeIsDummy(NODE, EQUALITY_SYSTEM) )
  {
     std::cout << "TRACKING: tracked column removed as (implied) free column singelton" << std::endl;
  }
#endif
#ifdef TRACK_ROW
   if(row == ROW && node_row == ROW_NODE && system_type == ROW_SYS && !nodeIsDummy(ROW_NODE, ROW_SYS))
   {
      std::cout << "TRACKING: removal of tracked row since it contained an (implied) free column singleton" << std::endl;
   }
#endif
   assert( !linking_row );
   assert( system_type == EQUALITY_SYSTEM );

   adaptObjectiveSubstitutedRow( system_type, node_row, row, linking_row, node_col, col );

   removeRow( system_type, node_row, row, linking_row );
   removeColumn( node_col, col, 0.0);
}

/* column col getting substituted with row row */ 
void PresolveData::adaptObjectiveSubstitutedRow( SystemType system_type, int node_row, int row, bool linking_row, int node_col, int col )
{
   assert( system_type == EQUALITY_SYSTEM );
   assert( !linking_row );
   assert( -1 <= node_row && node_row < nChildren );
   assert( -1 <= node_col && node_col < nChildren );

   BlockType block_col = ( node_col == -1) ? LINKING_VARS_BLOCK : CHILD_BLOCK; 
   const SparseStorageDynamic& col_mat_tp = getSparseGenMatrix(system_type, node_row, block_col)->getStorageDynamicTransposedRef();
   const double col_coef = col_mat_tp.getMat(col_mat_tp.getRowPtr(col).start);
   const double obj_coef = getSimpleVecColFromStochVec( *presProb->g, node_col)[col];

   assert( (col_mat_tp.getRowPtr(col).end - col_mat_tp.getRowPtr(col).start) == 1 );
   assert( row == col_mat_tp.getJcolM(col_mat_tp.getRowPtr(col).start) );
   assert( ! PIPSisZero(col_coef) );
   
   if( PIPSisZero(obj_coef) )
      return;

   const double rhs = getSimpleVecRowFromStochVec( *presProb->bA, node_row, CHILD_BLOCK)[row];
   
   const SparseStorageDynamic& a_mat = getSparseGenMatrix(system_type, node_row, LINKING_VARS_BLOCK)->getStorageDynamicRef();

   for(int i = a_mat.getRowPtr(row).start ; i < a_mat.getRowPtr(row).end; ++i)
   {
      const int col_idx = a_mat.getJcolM(i);

      if(col_idx != col || node_col != -1)
      {
         (*obj_vector_chgs)[col_idx] -= obj_coef * a_mat.getMat(i) / col_coef;
         outdated_objvector = true;
      }
   }

   if( node_row != -1 )
   {
      const SparseStorageDynamic& b_mat = getSparseGenMatrix(system_type, node_row, CHILD_BLOCK)->getStorageDynamicRef();
      
      for(int i = b_mat.getRowPtr(row).start; i < b_mat.getRowPtr(row).end; ++i)
      {
         const int col_idx = b_mat.getJcolM(i);

         if(col_idx != col || node_col != node_row)
            getSimpleVecColFromStochVec( *presProb->g, node_row)[col_idx] -= obj_coef * b_mat.getMat(i) / col_coef;
      }
   }

   obj_offset_chgs += obj_coef * rhs / col_coef;

   getSimpleVecColFromStochVec( *presProb->g, node_col)[col] = 0.0;
}

/* removes row from local system - sets rhs lhs and activities to zero */
void PresolveData::removeRow(SystemType system_type, int node, int row, bool linking)
{
   assert(-1 <= node && node <= nChildren);
   assert(!nodeIsDummy(node, system_type));

   if(linking)
   {
      assert(node == -1);

      /* Bl0 */
      removeRowFromMatrix(system_type, -1, LINKING_CONS_BLOCK, row);

      /* linking rows Bli */
      for(int child = 0; child < nChildren; ++child)
      {
         if(!nodeIsDummy(child, system_type))
            removeRowFromMatrix(system_type, child, LINKING_CONS_BLOCK, row);
      }
   }
   else
   {
      /* Amat */
      removeRowFromMatrix(system_type, node, LINKING_VARS_BLOCK, row);

      /* Bmat */
      if(node != -1)
         removeRowFromMatrix(system_type, node, CHILD_BLOCK, row);
   }


   BlockType block_type = ( linking ) ? LINKING_CONS_BLOCK : CHILD_BLOCK;

   /* set lhs rhs to zero */
   if(system_type == EQUALITY_SYSTEM)
      getSimpleVecRowFromStochVec(*presProb->bA, node, block_type)[row] = 0.0;
   else
   {
      getSimpleVecRowFromStochVec(*presProb->bl, node, block_type)[row] = 0.0;
      getSimpleVecRowFromStochVec(*presProb->bu, node, block_type)[row] = 0.0;
   }

   /* set activities and unbounded counters to zero */
   if(system_type == EQUALITY_SYSTEM)
   {
      getSimpleVecRowFromStochVec(*actmax_eq_part, node, block_type)[row] = 0.0;
      getSimpleVecRowFromStochVec(*actmin_eq_part, node, block_type)[row] = 0.0;
      getSimpleVecRowFromStochVec(*actmax_eq_ubndd, node, block_type)[row] = 0.0;
      getSimpleVecRowFromStochVec(*actmin_eq_ubndd, node, block_type)[row] = 0.0;

      if(linking)
      {
         (*actmax_eq_chgs)[row] = 0.0;
         (*actmin_eq_chgs)[row] = 0.0;
         (*actmax_eq_ubndd_chgs)[row] = 0.0;
         (*actmin_eq_ubndd_chgs)[row] = 0.0;
      }
   }
   else
   {
      getSimpleVecRowFromStochVec(*actmax_ineq_part, node, block_type)[row] = 0.0;
      getSimpleVecRowFromStochVec(*actmin_ineq_part, node, block_type)[row] = 0.0;
      getSimpleVecRowFromStochVec(*actmax_ineq_ubndd, node, block_type)[row] = 0.0;
      getSimpleVecRowFromStochVec(*actmin_ineq_ubndd, node, block_type)[row] = 0.0;

      if(linking)
      {
         (*actmax_ineq_chgs)[row] = 0.0;
         (*actmin_ineq_chgs)[row] = 0.0;
         (*actmax_ineq_ubndd_chgs)[row] = 0.0;
         (*actmin_ineq_ubndd_chgs)[row] = 0.0;
      }
   }


#ifndef NDEBUG
   /* assert non-zero counters of row are zero - only works for non-linking rows */
   if(system_type == EQUALITY_SYSTEM)
   {
      if(!linking)
         assert( getSimpleVecRowFromStochVec(*nnzs_row_A, node, block_type)[row] == 0.0 );   
   }
   else
   {
      if(!linking)
         assert( getSimpleVecRowFromStochVec(*nnzs_row_C, node, block_type)[row] == 0.0 );
   }  
#endif
}

void PresolveData::removeRowFromMatrix(SystemType system_type, int node, BlockType block_type, int row)
{
   assert(!nodeIsDummy(node, system_type));
   SparseGenMatrix* mat = getSparseGenMatrix(system_type, node, block_type);

   assert(mat);
   assert(mat->hasDynamicStorage());

   SparseStorageDynamic& mat_storage = mat->getStorageDynamicRef();
   SparseStorageDynamic& mat_transp_storage = mat->getStorageDynamicTransposedRef();
   assert( 0 <= row && row < mat_storage.getM());

   const int row_start = mat_storage.getRowPtr(row).start;
   const int row_end = mat_storage.getRowPtr(row).end;

   for(int k = row_start; k < row_end; k++)
   {
      const int col = mat_storage.getJcolM(k);

      mat_transp_storage.removeEntryAtRowCol(col, row);
      removeIndexColumn(node, block_type, col, 1);
   }

   removeIndexRow(system_type, node, block_type, row, mat_storage.getRowPtr(row).end - mat_storage.getRowPtr(row).start);
   mat_storage.clearRow(row);
}

bool PresolveData::verifyActivities()
{
   assert(!outdated_activities && linking_rows_need_act_computation == 0);

   bool activities_correct = true;

   StochVectorHandle actmax_eq_part_old(actmax_eq_part->cloneFull());
   StochVectorHandle actmin_eq_part_old(actmin_eq_part->cloneFull());

   StochVectorHandle actmax_eq_ubndd_old(actmax_eq_ubndd->cloneFull());
   StochVectorHandle actmin_eq_ubndd_old(actmin_eq_ubndd->cloneFull());

   StochVectorHandle actmax_ineq_part_old(actmax_ineq_part->cloneFull());
   StochVectorHandle actmin_ineq_part_old(actmin_ineq_part->cloneFull());

   StochVectorHandle actmax_ineq_ubndd_old(actmax_ineq_ubndd->cloneFull());
   StochVectorHandle actmin_ineq_ubndd_old(actmin_ineq_ubndd->cloneFull());

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
      if(my_rank == 14)
         std::cout << "actmax_eq_part not correct" << std::endl;
      activities_correct = false;
   }
   MPI_Barrier(MPI_COMM_WORLD);

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
   StochVectorHandle nnzs_col_new(nnzs_col->cloneFull());
   StochVectorHandle nnzs_row_A_new(nnzs_row_A->cloneFull());
   StochVectorHandle nnzs_row_C_new(nnzs_row_C->cloneFull());

   nnzs_col_new->setToZero();
   nnzs_row_A_new->setToZero();
   nnzs_row_C_new->setToZero();

   initNnzCounter(*nnzs_row_A_new, *nnzs_row_C_new, *nnzs_col_new);

   // linking variables:
   SimpleVector* nColOrigSimple = dynamic_cast<SimpleVector*>(nnzs_col_new->vec);
   SimpleVector* nColUpdatedSimple = dynamic_cast<SimpleVector*>(nnzs_col->vec);
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
      nColOrigSimple = dynamic_cast<SimpleVector*>(nnzs_col_new->children[it]->vec);
      nColUpdatedSimple = dynamic_cast<SimpleVector*>(nnzs_col->children[it]->vec);
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
   SimpleVector* nRowAOrigSimple = dynamic_cast<SimpleVector*>(nnzs_row_A_new->vec);
   SimpleVector* nRowAUpdatedSimple = dynamic_cast<SimpleVector*>(nnzs_row_A->vec);
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
      nRowAOrigSimple = dynamic_cast<SimpleVector*>(nnzs_row_A_new->children[it]->vec);
      nRowAUpdatedSimple = dynamic_cast<SimpleVector*>(nnzs_row_A->children[it]->vec);
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
      nRowAOrigSimple = dynamic_cast<SimpleVector*>(nnzs_row_A_new->vecl);
      nRowAUpdatedSimple = dynamic_cast<SimpleVector*>(nnzs_row_A->vecl);
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
   SimpleVector* nRowCOrigSimple = dynamic_cast<SimpleVector*>(nnzs_row_C_new->vec);
   SimpleVector* nRowCUpdatedSimple = dynamic_cast<SimpleVector*>(nnzs_row_C->vec);
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
      nRowCOrigSimple = dynamic_cast<SimpleVector*>(nnzs_row_C_new->children[it]->vec);
      nRowCUpdatedSimple = dynamic_cast<SimpleVector*>(nnzs_row_C->children[it]->vec);
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
      nRowCOrigSimple = dynamic_cast<SimpleVector*>(nnzs_row_C_new->vecl);
      nRowCUpdatedSimple = dynamic_cast<SimpleVector*>(nnzs_row_C->vecl);
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

bool PresolveData::nodeIsDummy(int node, SystemType system_type) const // todo change order
{
   assert( -1 <= node && node < nChildren );
   if( node == -1 )
      return false;
   StochGenMatrix& matrix = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<StochGenMatrix&>(*presProb->A) : dynamic_cast<StochGenMatrix&>(*presProb->C);
   // todo : asserts
   if( matrix.children[node]->isKindOf(kStochGenDummyMatrix))
   {
      assert( dynamic_cast<StochVector&>(*(presProb->bux)).children[node]->isKindOf(kStochDummy) );
      assert( dynamic_cast<StochVector&>(*(presProb->blx)).children[node]->isKindOf(kStochDummy) );

      if( system_type == EQUALITY_SYSTEM)
      {
         assert( dynamic_cast<StochVector&>(*(presProb->bA)).children[node]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->bux)).children[node]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->blx)).children[node]->isKindOf(kStochDummy) );
         assert( nnzs_row_A->children[node]->isKindOf(kStochDummy) );
      }
      else
      {
         assert( dynamic_cast<StochVector&>(*(presProb->bu)).children[node]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->bl)).children[node]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->icupp)).children[node]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->iclow)).children[node]->isKindOf(kStochDummy) );
         assert( nnzs_row_C->children[node]->isKindOf(kStochDummy) );
      }
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
   assert( !nodeIsDummy(node, system_type) );

   /* get upper and lower bound on variable */
   const SimpleVector& ixlow = getSimpleVecColFromStochVec(*(presProb->ixlow), (block_type == LINKING_VARS_BLOCK) ? -1 : node);
   const SimpleVector& xlow = getSimpleVecColFromStochVec(*(presProb->blx), (block_type == LINKING_VARS_BLOCK) ? -1 : node);
   const SimpleVector& ixupp = getSimpleVecColFromStochVec(*(presProb->ixupp), (block_type == LINKING_VARS_BLOCK) ? -1 : node);
   const SimpleVector& xupp = getSimpleVecColFromStochVec(*(presProb->bux), (block_type == LINKING_VARS_BLOCK) ? -1 : node);

   assert(ixlow[col] == 1.0);
   assert(ixupp[col] == 1.0);
   assert(PIPSisEQ(xlow[col], xupp[col], 1e-10));

   /* get unbounded counters */
   double* actmax_ubndd = (system_type == EQUALITY_SYSTEM) ? &getSimpleVecRowFromStochVec(*actmax_eq_ubndd, node, block_type)[row]
         : &getSimpleVecRowFromStochVec(*actmax_ineq_ubndd, node, block_type)[row];
   double* actmin_ubndd = (system_type == EQUALITY_SYSTEM) ? &getSimpleVecRowFromStochVec(*actmin_eq_ubndd, node, block_type)[row]
         : &getSimpleVecRowFromStochVec(*actmin_ineq_ubndd, node, block_type)[row];

   /* sum of reduction and counters */
   double actmax_ubndd_val = *actmax_ubndd;
   double actmin_ubndd_val = *actmin_ubndd;

   if(block_type == LINKING_CONS_BLOCK)
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
   double* actmax_part = (system_type == EQUALITY_SYSTEM) ? &getSimpleVecRowFromStochVec(*actmax_eq_part, node, block_type)[row]
         : &getSimpleVecRowFromStochVec(*actmax_ineq_part, node, block_type)[row];
   double* actmin_part = (system_type == EQUALITY_SYSTEM) ? &getSimpleVecRowFromStochVec(*actmin_eq_part, node, block_type)[row]
         : &getSimpleVecRowFromStochVec(*actmin_ineq_part, node, block_type)[row];

   if(block_type == LINKING_CONS_BLOCK)
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
      if(ixupp[col] == 1.0)
         (*actmax_part) -= coeff * xupp[col];
      else
      {
         assert(1 <= actmax_ubndd_val);
         --(*actmax_ubndd);
         if( actmax_ubndd_val == 2)
            computeRowMinOrMaxActivity(system_type, node, block_type, row, true);
      }

      if(ixlow[col] == 1.0)
         (*actmin_part) -= coeff * xlow[col];
      else
      {
         assert(1 <= actmin_ubndd_val);
         --(*actmin_ubndd);
         if( actmin_ubndd_val == 2)
            computeRowMinOrMaxActivity(system_type, node, block_type, row, false);
      }
   }
   else
   {
      if(ixlow[col] == 1.0)
         (*actmax_part) -= coeff*xlow[col];
      else
      {
         assert(1 <= actmax_ubndd_val);
         --(*actmax_ubndd);
         if( actmax_ubndd_val == 2)
            computeRowMinOrMaxActivity(system_type, node, block_type, row, true);
      }

      if(ixupp[col] == 1.0)
         (*actmin_part) -= coeff*xupp[col];
      else
      {
         assert(1 <= actmin_ubndd_val);
         --(*actmin_ubndd);
         if( actmin_ubndd_val == 2)
            computeRowMinOrMaxActivity(system_type, node, block_type, row, false);
      }
   }

   if( block_type == LINKING_CONS_BLOCK)
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

      if( nodeIsDummy(node, system_type) )
         continue;

      /* get upper, lower bounds */
      const SimpleVector& ixlow = getSimpleVecColFromStochVec(*(presProb->ixlow), node);
      const SimpleVector& xlow = getSimpleVecColFromStochVec(*(presProb->blx), node);
      const SimpleVector& ixupp = getSimpleVecColFromStochVec(*(presProb->ixupp), node);
      const SimpleVector& xupp = getSimpleVecColFromStochVec(*(presProb->bux), node);

      /* get matrix */
      SparseStorageDynamic& mat = getSparseGenMatrix(system_type, node, LINKING_CONS_BLOCK)->getStorageDynamicRef();

      for( int j = mat.getRowPtr(row).start; j < mat.getRowPtr(row).end; j++ )
      {
         const int col = mat.getJcolM(j);
         const double entry = mat.getMat(j);

         assert( 0 <= col && col < ixlow.n);
         assert( !PIPSisZero(entry));

         if( PIPSisLT(0.0, entry) )
         {
            if( !upper && ixlow[col] != 0.0 )
               act_part += entry * xlow[col];
            if( upper && ixupp[col] != 0.0 )
               act_part += entry * xupp[col];
         }
         else
         {
            if( !upper && ixupp[col] != 0.0 )
               act_part += entry * xupp[col];
            if( upper && ixlow[col] != 0.0 )
               act_part += entry * xlow[col];
         }
      }
   }

   return act_part;
}

/// recomputes and updates the activity of a locally available row or increases the recompute counter for linking constraints
/// does not compute activities for linking rows
void PresolveData::computeRowMinOrMaxActivity(SystemType system_type, int node, BlockType block_type, int row, bool upper)
{
#ifdef TRACK_ROW
   if( ROW_NODE == node && ROW == row && ROW_SYS == system_type &&
         ( (block_type == LINKING_CONS_BLOCK && ROW_IS_LINK) || (!ROW_IS_LINK)) && (my_rank == 0 || node != -1) )
   {
      if( block_type == LINKING_CONS_BLOCK)
         std::cout << "TRACKING ROW: tracked row is linking constraint and needs " << (upper ? "upper" : "lower") << " activity recomputation " << std::endl;
      else
         std::cout << "TRACKING ROW: computing " << (upper ? "upper" : "lower") << " activity of tracked row " << std::endl;
   }
#endif

   /* single linking rows that have to be computed for the first time */
   if( block_type == LINKING_CONS_BLOCK )
   {
      ++linking_rows_need_act_computation;
      return;
   }

   /* upper lower bounds linking vars */
   const SimpleVector& ixlow_root = getSimpleVecColFromStochVec(*(presProb->ixlow), -1);
   const SimpleVector& xlow_root = getSimpleVecColFromStochVec(*(presProb->blx), -1);
   const SimpleVector& ixupp_root = getSimpleVecColFromStochVec(*(presProb->ixupp), -1);
   const SimpleVector& xupp_root = getSimpleVecColFromStochVec(*(presProb->bux), -1);

   /* get upper, lower bounds */
   const SimpleVector& ixlow = getSimpleVecColFromStochVec(*(presProb->ixlow), node);
   const SimpleVector& xlow = getSimpleVecColFromStochVec(*(presProb->blx), node);
   const SimpleVector& ixupp = getSimpleVecColFromStochVec(*(presProb->ixupp), node);
   const SimpleVector& xupp = getSimpleVecColFromStochVec(*(presProb->bux), node);

   /* get matrix */
   SparseStorageDynamic& Amat = getSparseGenMatrix(system_type, node, LINKING_VARS_BLOCK)->getStorageDynamicRef();
   SparseStorageDynamic& Bmat = getSparseGenMatrix(system_type, node, CHILD_BLOCK)->getStorageDynamicRef();

   /* get activity vector */
   SimpleVector* act_vec;
   if(system_type == EQUALITY_SYSTEM)
   {
      act_vec = (!upper) ? &getSimpleVecRowFromStochVec(*actmin_eq_part, node, block_type) :
            &getSimpleVecRowFromStochVec(*actmax_eq_part, node, block_type);
   }
   else
   {
      act_vec = (!upper) ? &getSimpleVecRowFromStochVec(*actmin_ineq_part, node, block_type) :
            &getSimpleVecRowFromStochVec(*actmax_ineq_part, node, block_type);
   }

   double& act_part = (*act_vec)[row];
   act_part = 0;

   for( int j = Amat.getRowPtr(row).start; j < Amat.getRowPtr(row).end; j++)
   {
      const int col = Amat.getJcolM(j);
      const double entry = Amat.getMat(j);

      assert( 0 <= col && col < ixlow_root.n );
      assert( !PIPSisZero(entry) );

      if( PIPSisLT(0.0, entry) )
      {
         if( !upper && ixlow_root[col] != 0.0)
            act_part += entry * xlow_root[col];
         if( upper && ixupp_root[col] != 0.0)
            act_part += entry * xupp_root[col];
      }
      else
      {
         if( !upper && ixupp_root[col] != 0.0 )
            act_part += entry * xupp_root[col];
         if( upper && ixlow_root[col] != 0.0 )
            act_part += entry * xlow_root[col];
      }
   }

   if( node != -1 )
   {
      for( int j = Bmat.getRowPtr(row).start; j < Bmat.getRowPtr(row).end; j++ )
      {
         const int col = Bmat.getJcolM(j);
         const double entry = Bmat.getMat(j);

         assert( 0 <= col && col < ixlow.n);
         assert( !PIPSisZero(entry));

         if( PIPSisLT(0.0, entry) )
         {
            if( !upper && ixlow[col] != 0.0 )
               act_part += entry * xlow[col];
            if( upper && ixupp[col] != 0.0 )
               act_part += entry * xupp[col];
         }
         else
         {
            if( !upper && ixupp[col] != 0.0 )
               act_part += entry * xupp[col];
            if( upper && ixlow[col] != 0.0 )
               act_part += entry * xlow[col];
         }
      }
   }

#ifdef TRACK_ROW
   if( ROW_NODE == node && ROW == row && ROW_SYS == system_type &&
         ( (block_type == LINKING_CONS_BLOCK && ROW_IS_LINK) || (!ROW_IS_LINK)) && (my_rank == 0 || node != -1) )
         std::cout << "TRACKING ROW: new " << (upper ? "upper" : "lower") << " activity is " << act_part << std::endl;
#endif
}

/** updates the bounds on a variable as well as activities */
bool PresolveData::updateBoundsVariable(int node, int col, double xu, double xl)
{
   SimpleVector& ixlow = getSimpleVecColFromStochVec(*(presProb->ixlow), node);
   SimpleVector& xlow = getSimpleVecColFromStochVec(*(presProb->blx), node);
   SimpleVector& ixupp = getSimpleVecColFromStochVec(*(presProb->ixupp), node);
   SimpleVector& xupp = getSimpleVecColFromStochVec(*(presProb->bux), node);

#ifdef TRACK_COLUMN
   if( NODE == node && COLUMN == col && (my_rank == 0 || node != -1) && !nodeIsDummy(NODE, EQUALITY_SYSTEM) )
      std::cout << "TRACKING COLUMN: updating column bounds from [" << xlow[col] << ", " << xupp[col] << "] for column " << COLUMN << " node " << node << " with [" <<
            xl << ", " << xu << "]" << std::endl;
#endif

   if( ( ixlow[col] != 0.0 && PIPSisLT(xu, xlow[col]) )
         || (ixupp[col] != 0.0 && PIPSisLT(xupp[col], xl) )
         || PIPSisLT(xu, xl) )
      abortInfeasible(MPI_COMM_WORLD, "Varbounds update detected infeasible new bounds!", "PresolveData.C", "updateBoundsVariable");

   bool updated = false;
   double xu_old = std::numeric_limits<double>::infinity();
   double xl_old = -std::numeric_limits<double>::infinity();

   if( xu < std::numeric_limits<double>::max() && (ixupp[col] == 0.0 || PIPSisLT(xu, xupp[col])) )
   {
      updated = true;
      if(ixupp[col] == 1.0)
         xu_old = xupp[col];

      xupp[col] = xu;
      ixupp[col] = 1.0;
   }

   if( -std::numeric_limits<double>::max() < xl && (ixlow[col] == 0.0 || PIPSisLE( xlow[col], xl)) )
   {
      updated = true;
      if(ixlow[col] == 1.0)
         xl_old = xlow[col];

      xlow[col] = xl;
      ixlow[col] = 1.0;
   }

   if( updated )
   {
#ifdef TRACK_COLUMN
      if( NODE == node && COLUMN == col && (my_rank == 0 || node != -1) && !nodeIsDummy(NODE, EQUALITY_SYSTEM) )
      {
         std::cout << "TRACKING COLUMN: bounds are now [" << xlow[col] << ", " << xupp[col] << "]" << std::endl;
         std::cout << "TRACKING COLUMN: moving on to update activities" << std::endl;
      }
#endif
      updateRowActivities(node, col, xu, xl, xu_old, xl_old);
   }
   else
   {
#ifdef TRACK_COLUMN
      if( NODE == node && COLUMN == col && (my_rank == 0 || node != -1) && !nodeIsDummy(NODE, EQUALITY_SYSTEM) )
         std::cout << "TRACKING COLUMN: col " << COLUMN << " was not updated" << std::endl;
#endif
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
   assert(-1 <= node && node <= nChildren);

#ifdef TRACK_COLUMN
   if( NODE == node && COLUMN == col && (my_rank == 0 || node != -1) && !nodeIsDummy(NODE, EQUALITY_SYSTEM) )
   {
      std::cout << "TRACKING COLUMN: col " << COLUMN << " node " << NODE << " gets it's activities updated" << std::endl;
      std::cout << "\t bounds changed from [" << old_ubx << ", " << old_lbx << "] to [" << lbx << ", " << ubx << "]" << std::endl;
   }
#endif

   /* if node == -1 go through all linking var blocks of both systems */
   if( node == -1 )
   {
#ifdef TRACK_COLUMN
      if( NODE == node && COLUMN == col && (my_rank == 0 || node != -1) && !nodeIsDummy(NODE, EQUALITY_SYSTEM) )
         std::cout << "TRACKING COLUMN: the column is linking (root node)" << std::endl;
#endif
      /* A0/B0 and C0/D0block */
      updateRowActivitiesBlock(EQUALITY_SYSTEM, -1, LINKING_VARS_BLOCK, col, ubx, old_ubx, true);
      updateRowActivitiesBlock(EQUALITY_SYSTEM, -1, LINKING_VARS_BLOCK, col, lbx, old_lbx, false);

      updateRowActivitiesBlock(INEQUALITY_SYSTEM, -1, LINKING_VARS_BLOCK, col, ubx, old_ubx, true);
      updateRowActivitiesBlock(INEQUALITY_SYSTEM, -1, LINKING_VARS_BLOCK, col, lbx, old_lbx, false);

      /* Bl0 and Dl0 */
      updateRowActivitiesBlock(EQUALITY_SYSTEM, -1, LINKING_CONS_BLOCK, col, ubx, old_ubx, true);
      updateRowActivitiesBlock(EQUALITY_SYSTEM, -1, LINKING_CONS_BLOCK, col, lbx, old_lbx, false);

      updateRowActivitiesBlock(INEQUALITY_SYSTEM, -1, LINKING_CONS_BLOCK, col, ubx, old_ubx, true);
      updateRowActivitiesBlock(INEQUALITY_SYSTEM, -1, LINKING_CONS_BLOCK, col, lbx, old_lbx, false);

      for(int child = 0; child < nChildren; ++child)
      {
         /* Ai and Ci */
         updateRowActivitiesBlock(EQUALITY_SYSTEM, child, LINKING_VARS_BLOCK, col, ubx, old_ubx, true);
         updateRowActivitiesBlock(EQUALITY_SYSTEM, child, LINKING_VARS_BLOCK, col, lbx, old_lbx, false);

         updateRowActivitiesBlock(INEQUALITY_SYSTEM, child, LINKING_VARS_BLOCK, col, ubx, old_ubx, true);
         updateRowActivitiesBlock(INEQUALITY_SYSTEM, child, LINKING_VARS_BLOCK, col, lbx, old_lbx, false);
      }
   }
   else
   {
#ifdef TRACK_COLUMN
      if( NODE == node && COLUMN == col && (my_rank == 0 || node != -1) && !nodeIsDummy(NODE, EQUALITY_SYSTEM) )
         std::cout << "TRACKING COLUMN: the column is non-linking (non-root)" << std::endl;
#endif
      /* Bmat, Blmat */
      /* Bmat */
      updateRowActivitiesBlock(EQUALITY_SYSTEM, node, CHILD_BLOCK, col, ubx, old_ubx, true);
      updateRowActivitiesBlock(EQUALITY_SYSTEM, node, CHILD_BLOCK, col, lbx, old_lbx, false);

      /* Blmat */
      updateRowActivitiesBlock(EQUALITY_SYSTEM, node, LINKING_CONS_BLOCK, col, ubx, old_ubx, true);
      updateRowActivitiesBlock(EQUALITY_SYSTEM, node, LINKING_CONS_BLOCK, col, lbx, old_lbx, false);

      /* Dmat Dlmat */

      /* Dmat */
      updateRowActivitiesBlock(INEQUALITY_SYSTEM, node, CHILD_BLOCK, col, ubx, old_ubx, true);
      updateRowActivitiesBlock(INEQUALITY_SYSTEM, node, CHILD_BLOCK, col, lbx, old_lbx, false);

      /* Dlmat */
      updateRowActivitiesBlock(INEQUALITY_SYSTEM, node, LINKING_CONS_BLOCK, col, ubx, old_ubx, true);
      updateRowActivitiesBlock(INEQUALITY_SYSTEM, node, LINKING_CONS_BLOCK, col, lbx, old_lbx, false);
   }
}

void PresolveData::updateRowActivitiesBlock(SystemType system_type, int node, BlockType block_type, int col, double bound, double old_bound, bool upper)
{
   assert(-1 <= node && node < nChildren);
   assert(0 <= col);

   /* dummies do not adjust anything */
   if( nodeIsDummy(node, system_type) )
      return;

   /* we do not have to adjust activities if no new bound was found */
   if( bound == std::numeric_limits<double>::infinity() || bound == -std::numeric_limits<double>::infinity() )
      return;

   /* no changes if a worse bound has been found */
   if( (upper && !PIPSisLT( bound, old_bound)) || (!upper && PIPSisLT(bound, old_bound)) )
      return;

   SparseStorageDynamic& mat_transp = getSparseGenMatrix(system_type, node, block_type)->getStorageDynamicTransposedRef();

   assert(col < mat_transp.getM() );

   SimpleVector& actmax_ubndd = (system_type == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmax_eq_ubndd, node, block_type) :
         getSimpleVecRowFromStochVec(*actmax_ineq_ubndd, node, block_type);
   SimpleVector& actmin_ubndd = (system_type == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmin_eq_ubndd, node, block_type) :
         getSimpleVecRowFromStochVec(*actmin_ineq_ubndd, node, block_type);
   SimpleVector& actmax_part = (system_type == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmax_eq_part, node, block_type) :
         getSimpleVecRowFromStochVec(*actmax_ineq_part, node, block_type);
   SimpleVector& actmin_part = (system_type == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmin_eq_part, node, block_type) :
         getSimpleVecRowFromStochVec(*actmin_ineq_part, node, block_type);

   /* we always set these variables but only use them in case an actual linking row gets it's activities updated */
   SimpleVector& actmax_ubndd_chgs = (system_type == EQUALITY_SYSTEM) ? *actmax_eq_ubndd_chgs : *actmax_ineq_ubndd_chgs;
   SimpleVector& actmin_ubndd_chgs = (system_type == EQUALITY_SYSTEM) ? *actmin_eq_ubndd_chgs : *actmin_ineq_ubndd_chgs;
   SimpleVector& actmax_part_chgs = (system_type == EQUALITY_SYSTEM) ? *actmax_eq_chgs : *actmax_ineq_chgs;
   SimpleVector& actmin_part_chgs = (system_type == EQUALITY_SYSTEM) ? *actmin_eq_chgs : *actmin_ineq_chgs;

#ifdef TRACK_COLUMN
   if( NODE == node && COLUMN == col && (my_rank == 0 || node != -1) && !nodeIsDummy(NODE, EQUALITY_SYSTEM) )
      std::cout << "TRACKING COLUMN: updating activities column " << col << " node " << node << " system " << ( (system_type == EQUALITY_SYSTEM) ? "EQ_SYS" : "INEQ_SYS" ) <<
         " with new " << ( upper ? "upper" : "lower" ) << " bound " << bound << " from " << old_bound << std::endl;
#endif

   for( int j = mat_transp.getRowPtr(col).start; j < mat_transp.getRowPtr(col).end; ++j )
   {
      const int row = mat_transp.getJcolM(j);
      const double entry = mat_transp.getMat(j);

      assert( !PIPSisZero(entry) );

      /* get affected partial activity and act_ubndd */
      bool switch_upperlower = upper;
      if( PIPSisLE(entry, 0.0) )
         switch_upperlower = !upper;

      SimpleVector& act_part = (switch_upperlower) ? actmax_part : actmin_part;
      SimpleVector& act_part_chgs = (switch_upperlower) ? actmax_part_chgs : actmin_part_chgs;
      SimpleVector& act_ubndd = (switch_upperlower) ? actmax_ubndd : actmin_ubndd;
      SimpleVector& act_ubndd_chgs = (switch_upperlower) ? actmax_ubndd_chgs : actmin_ubndd_chgs;


#ifdef TRACK_ROW
   if( ROW_NODE == node && ROW == row && ROW_SYS == system_type &&
         ( (block_type == LINKING_CONS_BLOCK && ROW_IS_LINK) || (!ROW_IS_LINK)) && (my_rank == 0 || node != -1) )
   {
      std::cout << "TRACKING ROW: activities of tracked row are getting updated. Current state:" << std::endl;
      if( ROW_BLOCK == LINKING_CONS_BLOCK)
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
#endif

      /* if the old bound was not set we have to modify the unbounded counters */
      if( old_bound == std::numeric_limits<double>::infinity() || old_bound == -std::numeric_limits<double>::infinity() )
      {
#ifdef TRACK_ROW
            if( ROW_NODE == node && ROW == row && ROW_SYS == system_type &&
                  ( (block_type == LINKING_CONS_BLOCK && ROW_IS_LINK) || (!ROW_IS_LINK)) && (my_rank == 0 || node != -1) )
                  std::cout << "TRACKING ROW: ubndd counters are being changed" << std::endl;
#endif

         /* every process works the root node - even in the linking row case */
         if( node != -1 && block_type == LINKING_CONS_BLOCK )
         {
            --act_ubndd_chgs[row];
            outdated_activities = true;
         }
         else
            --act_ubndd[row];

         if(block_type == LINKING_CONS_BLOCK)
            assert( 0 <= act_ubndd[row] + act_ubndd_chgs[row]);
         else
            assert( 0 <= act_ubndd[row]);

         /* from now on activity of row has to be computed */
         if( (block_type == LINKING_CONS_BLOCK && (act_ubndd[row] + act_ubndd_chgs[row] == 1))
               || (act_ubndd[row] == 1 && block_type != LINKING_CONS_BLOCK) )
         {
#ifdef TRACK_ROW
            if( ROW_NODE == node && ROW == row && ROW_SYS == system_type &&
                  ( (block_type == LINKING_CONS_BLOCK && ROW_IS_LINK) || (!ROW_IS_LINK)) && (my_rank == 0 || node != -1) )
            {
               if( ROW_IS_LINK )
                  std::cout << "TRACKING ROW: " << ( switch_upperlower ? "upper" : "lower" ) << " now needs act computation " << std::endl;
               else
                  std::cout << "TRACKING ROW: first time computation of " << ( switch_upperlower ? "upper" : "lower" ) << " activity " << std::endl;
            }
#endif
            if(block_type == LINKING_CONS_BLOCK)
               ++linking_rows_need_act_computation;
            else
               computeRowMinOrMaxActivity(system_type, node, block_type, row, switch_upperlower);
         }
         else if( (block_type == LINKING_CONS_BLOCK && act_ubndd[row] + act_ubndd_chgs[row] == 0) ||
               (block_type != LINKING_CONS_BLOCK && act_ubndd[row] == 0) )
         {
#ifdef TRACK_ROW
            if( ROW_NODE == node && ROW == row && ROW_SYS == system_type &&
                  ( (block_type == LINKING_CONS_BLOCK && ROW_IS_LINK) || (!ROW_IS_LINK)) && (my_rank == 0 || node != -1) )
                  std::cout << "TRACKING ROW: adjusting activities " << std::endl;
#endif
            if( node != -1 && block_type == LINKING_CONS_BLOCK)
               act_part_chgs[row] += bound * entry;
            else
               act_part[row] += bound * entry;
         }
      }
      else
      {
#ifdef TRACK_ROW
            if( ROW_NODE == node && ROW == row && ROW_SYS == system_type &&
                  ( (block_type == LINKING_CONS_BLOCK && ROW_IS_LINK) || (!ROW_IS_LINK)) && (my_rank == 0 || node != -1) )
                  std::cout << "TRACKING ROW: " << " no changes in ubndd counters - only adjust activities" << std::endl;
#endif
         /* better bound was found */
         if( (block_type == LINKING_CONS_BLOCK && act_ubndd[row] + act_ubndd_chgs[row] <= 1) ||
               (block_type != LINKING_CONS_BLOCK && act_ubndd[row] <= 1) )
         {

            if( node != -1 && block_type == LINKING_CONS_BLOCK)
            {
               outdated_activities = true;
               act_part_chgs[row] += (bound - old_bound) * entry;
            }
            else
               act_part[row] += (bound - old_bound) * entry;
         }
      }

#ifdef TRACK_ROW
      if( ROW_NODE == node && ROW == row && ROW_SYS == system_type &&
            ( (block_type == LINKING_CONS_BLOCK && ROW_IS_LINK) || (!ROW_IS_LINK)) && (my_rank == 0 || node != -1) )
      {
         std::cout << "TRACKING ROW: activities of tracked row have been updated: " << std::endl;
         if( ROW_BLOCK == LINKING_CONS_BLOCK)
         {
            std::cout << "\tnow linking row with " << ( switch_upperlower ? "upper" : "lower" ) << " activity " << act_part[row] << " + " << act_part_chgs[row] << "(changes) and " <<
                  (switch_upperlower ? "upper" : "lower" ) << " unbounded counters " << act_ubndd[row] << " + " << act_ubndd_chgs[row] << "(changes)" << std::endl;
         }
         else
         {
            std::cout << "\tnow non-linking row with " << ( switch_upperlower ? "upper" : "lower" ) << " activity " << act_part[row] << " and " << (switch_upperlower ? "upper" : "lower" ) <<
                  " unbounded counters " << act_ubndd[row] << std::endl;
         }
      }
#endif
   }
}

/// returns activity and unbounded counters for specified row
/// +/- infinity if there are two or more unbouded entries in a row
void PresolveData::getRowActivities(SystemType system_type, int node, BlockType block_type, int row, double& max_act,
      double& min_act, int& max_ubndd, int& min_ubndd) const
{
   assert(!nodeIsDummy(node, system_type));

   double row_max_counters = (system_type == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmax_eq_ubndd, node, block_type)[row]
         : getSimpleVecRowFromStochVec(*actmax_ineq_ubndd, node, block_type)[row];
   double row_min_counters = (system_type == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmin_eq_ubndd, node, block_type)[row]
         : getSimpleVecRowFromStochVec(*actmin_ineq_ubndd, node, block_type)[row];

   if(block_type == LINKING_CONS_BLOCK)
   {
      row_max_counters += (system_type == EQUALITY_SYSTEM) ? (*actmax_eq_ubndd_chgs)[row]
            : (*actmax_ineq_ubndd_chgs)[row];
      row_min_counters += (system_type == EQUALITY_SYSTEM) ? (*actmin_eq_ubndd_chgs)[row]
            : (*actmin_ineq_ubndd_chgs)[row];
   }

   max_ubndd = static_cast<int>(row_max_counters);
   min_ubndd = static_cast<int>(row_min_counters);

#ifndef NDEBUG
   /// asserting that cast from double to int is correct
   double int_part_max, int_part_min;
   std::modf(row_max_counters, &int_part_max);
   std::modf(row_min_counters, &int_part_min);
   assert(int_part_max == max_ubndd);
   assert(int_part_min == min_ubndd);
   assert(max_ubndd >= 0);
   assert(min_ubndd >= 0);
#endif

   max_act = (system_type == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmax_eq_part, node, block_type)[row]
         : getSimpleVecRowFromStochVec(*actmax_ineq_part, node, block_type)[row];
   min_act = (system_type == EQUALITY_SYSTEM) ? getSimpleVecRowFromStochVec(*actmin_eq_part, node, block_type)[row]
         : getSimpleVecRowFromStochVec(*actmin_ineq_part, node, block_type)[row];

   if(block_type == LINKING_CONS_BLOCK)
   {
      max_act += (system_type == EQUALITY_SYSTEM) ? (*actmax_eq_chgs)[row] : (*actmax_ineq_chgs)[row];
      min_act += (system_type == EQUALITY_SYSTEM) ? (*actmin_eq_chgs)[row] : (*actmin_ineq_chgs)[row];
   }

   if( row_max_counters >= 2)
      assert( max_act == std::numeric_limits<double>::infinity() );
   else
   {
      if( block_type != LINKING_CONS_BLOCK)
      {
         if( !(max_act < std::numeric_limits<double>::infinity()))
            std::cout << row << " " << block_type << " " << node << " " << system_type << std::endl;
         assert(max_act < std::numeric_limits<double>::infinity());
      }
   }

   if( row_min_counters >= 2)
      assert( min_act == -std::numeric_limits<double>::infinity() );
   else
   {
      if( block_type != LINKING_CONS_BLOCK)
      {
         assert(-std::numeric_limits<double>::infinity() < min_act);
      }
   }
}

SimpleVector& PresolveData::getSimpleVecRowFromStochVec(const StochVector& stochvec, int node, BlockType block_type) const
{
   assert(-1 <= node && node < nChildren);

   if(node == -1)
   {
      if(block_type == LINKING_CONS_BLOCK)
      {
         assert(stochvec.vecl);
         return dynamic_cast<SimpleVector&>(*(stochvec.vecl));
      }
      else
      {
         assert(stochvec.vec);
         return dynamic_cast<SimpleVector&>(*(stochvec.vec));
      }
   }
   else
   {
      if(block_type == CHILD_BLOCK || block_type == LINKING_VARS_BLOCK)
      {
         assert(stochvec.children[node]->vec);
         return dynamic_cast<SimpleVector&>(*(stochvec.children[node]->vec));
      }
      else
      {
         assert(stochvec.vecl);
         return dynamic_cast<SimpleVector&>(*(stochvec.vecl));
      }
   }
}

SimpleVector& PresolveData::getSimpleVecColFromStochVec(const StochVector& stochvec, int node) const
{
   assert(-1 <= node && node < nChildren);

   if(node == -1)
   {
      assert(stochvec.vecl == NULL);
      assert(stochvec.vec);
      return dynamic_cast<SimpleVector&>(*(stochvec.vec));
   }
   else
   {
      assert(stochvec.children[node]->vecl == NULL);
      assert(stochvec.children[node]->vec);
      return dynamic_cast<SimpleVector&>(*(stochvec.children[node]->vec));
   }
}

SparseGenMatrix* PresolveData::getSparseGenMatrix(SystemType system_type, int node, BlockType block_type) const
{
   assert( -1 <= node && node <= nChildren );
   assert(!nodeIsDummy(node, system_type));

   StochGenMatrix& sMat = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<StochGenMatrix&>(*presProb->A) : dynamic_cast<StochGenMatrix&>(*presProb->C);

   if(node == -1)
      return (block_type == LINKING_CONS_BLOCK) ? sMat.Blmat : sMat.Bmat;
   else
   {
      if(block_type == LINKING_VARS_BLOCK)
         return sMat.children[node]->Amat;
      else if(block_type == CHILD_BLOCK)
         return sMat.children[node]->Bmat;
      else if(block_type == LINKING_CONS_BLOCK)
         return sMat.children[node]->Blmat;
   }
   return NULL;
}

/* only prints the part of a linking constraint the current process knows about */
// todo does not yet print linking constraints
void PresolveData::writeRowLocalToStreamDense(std::ostream& out, SystemType system_type, int node, BlockType block_type, int row) const
{
   if(nodeIsDummy(node, system_type))
      return;

   if(node == -1 && block_type != LINKING_CONS_BLOCK && my_rank != 0)
      return;

   out << "SystemType: " << system_type << "\tnode: " << node << "\tLinkingCons: " << (block_type == LINKING_CONS_BLOCK) << "\trow: " << row << std::endl;

   if(system_type == EQUALITY_SYSTEM)
   {
      out << getSimpleVecRowFromStochVec(*presProb->bA, node, block_type)[row] << " = ";
   }
   else
   {
      double iclow = getSimpleVecRowFromStochVec(*presProb->iclow, node, block_type)[row];
      double clow = (iclow == 1.0) ? getSimpleVecRowFromStochVec(*presProb->bl, node, block_type)[row] : -std::numeric_limits<double>::infinity();

      out << clow << " <= ";
   }

   if(node != -1 && block_type != LINKING_CONS_BLOCK)
   {
      block_type = CHILD_BLOCK;
      writeMatrixRowToStreamDense(out, *getSparseGenMatrix(system_type, node, LINKING_VARS_BLOCK), node, row, getSimpleVecColFromStochVec(*presProb->ixupp, -1),
            getSimpleVecColFromStochVec(*presProb->bux, -1), getSimpleVecColFromStochVec(*presProb->ixlow, -1), getSimpleVecColFromStochVec(*presProb->blx, -1));
   }

   writeMatrixRowToStreamDense(out, *getSparseGenMatrix(system_type, node, block_type), node, row, getSimpleVecColFromStochVec(*presProb->ixupp, node),
         getSimpleVecColFromStochVec(*presProb->bux, node),getSimpleVecColFromStochVec(*presProb->ixlow, node),getSimpleVecColFromStochVec(*presProb->blx, node));


   if(system_type == INEQUALITY_SYSTEM)
   {
      double icupp = getSimpleVecRowFromStochVec(*presProb->icupp, node, block_type)[row];
      double cupp = (icupp == 1.0) ? getSimpleVecRowFromStochVec(*presProb->bu, node, block_type)[row] : std::numeric_limits<double>::infinity();

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

      out << " + " << val << " * x_" << node << "_" << col << " € [" << ( (ixlow[col] == 0.0) ? -std::numeric_limits<double>::infinity() : xlow[col]) << ", "
            << ( (ixupp[col] == 0.0) ? std::numeric_limits<double>::infinity() : xupp[col]) << "]";
   }
}

/* stupid slow - do not use this for things.. */
StochVectorHandle PresolveData::getRowAsStochVector(SystemType system_type, int node, int row, bool linking_row)
{
   // create new StochVector and set pattern according to row   
   StochVector* svec;
   if(!linking_row)
   {
      svec = new StochVector( nnzs_col->vec->length(), MPI_COMM_WORLD, -1);
      const SparseStorageDynamic& amat = getSparseGenMatrix(system_type, node, LINKING_VARS_BLOCK)->getStorageDynamicRef();
      // todo assert size is correct
      
      /* copy Amat */
      for(int i = amat.getRowPtr(row).start; i < amat.getRowPtr(row).end; ++i)
      {
        dynamic_cast<SimpleVector&>(*svec->vec)[amat.getJcolM(i)] = amat.getMat(i);
      }
      
      /* copy bmat and add dummy children */
      for(int child = 0; child < nChildren; ++child)
      {
        if(child != node)
        {
           svec->AddChild( new StochDummyVector() );
        }
        else
        {
          StochVector* child_vec = new StochVector( nnzs_col->children[node]->length(), MPI_COMM_WORLD, -1);
          const SparseStorageDynamic& bmat = getSparseGenMatrix(system_type, node, CHILD_BLOCK)->getStorageDynamicRef();
          // todo assert size is correct
          
          for(int i = bmat.getRowPtr(row).start; i < bmat.getRowPtr(row).end; ++i)
          {
            dynamic_cast<SimpleVector&>(*child_vec->vec)[bmat.getJcolM(i)] = bmat.getMat(i);
          }
          svec->AddChild(child_vec);
        }
      }
   }
   else
   {
     svec = nnzs_col->clone();
     /* copy all blmat */
     for(int child = -1; child < nChildren; ++child)
     {
       if( !nodeIsDummy(child, system_type) )
       {
          const SparseStorageDynamic& blmat = getSparseGenMatrix(system_type, child, LINKING_CONS_BLOCK)->getStorageDynamicRef();
          // todo assert size is correct
          
          for(int i = blmat.getRowPtr(row).start; i < blmat.getRowPtr(row).end; ++i)
          {
            getSimpleVecColFromStochVec(*svec, child)[blmat.getJcolM(i)] = blmat.getMat(i);
          }
       }
     }
   }
   return StochVectorHandle(svec);
}

StochVectorHandle PresolveData::getColAsStochVector(SystemType system_type, int node, int col)
{
  // create new StochVector and set pattern according to row
  // todo: we need two vectors to display one column
  assert(false);
  return NULL;
}

void PresolveData::printVarBoundStatistics(std::ostream& out) const
{
   const StochVector& xupp = dynamic_cast<const StochVector&>(*presProb->bux);
   const StochVector& xlow = dynamic_cast<const StochVector&>(*presProb->blx);
   const StochVector& ixupp = dynamic_cast<const StochVector&>(*presProb->ixupp);
   const StochVector& ixlow = dynamic_cast<const StochVector&>(*presProb->ixlow);

   StochVectorHandle xlow_def = StochVectorHandle(xlow.cloneFull());
   StochVectorHandle xupp_def = StochVectorHandle(xupp.cloneFull());

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
