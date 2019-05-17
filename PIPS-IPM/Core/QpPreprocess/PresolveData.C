/*
 * PresolveData.C
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

#include "PresolveData.h"
#include "StochGenMatrix.h"
#include "DoubleMatrixTypes.h"
#include "pipsdef.h"

#include <stdexcept>
#include <limits>
#include <algorithm>
#include <string>
#include <cmath>

PresolveData::PresolveData(const sData* sorigprob, StochPostsolver* postsolver) :
      postsolver(postsolver),
      outdated_activities(true),
      outdated_lhsrhs(false),
      outdated_nnzs(false),
      outdated_linking_var_bounds(false),
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
      nChildren(nnzs_col->children.size()),
      objOffset(0.0), obj_offset_chgs(0.0),
      elements_deleted(0), elements_deleted_transposed(0)
{
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

   double* array_act_unbounded_chgs = new double[lenght_array_act_chgs]();

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
   initNnzCounter();
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
      if(vec_ixlow.elements()[i] == 0.0)
         vec_xlow.elements()[i] = -value;

      if(vec_ixupp.elements()[i] == 0.0)
         vec_xupp.elements()[i] = value;
   }
}

sData* PresolveData::finalize()
{

#ifndef NDEBUG
   if(distributed)
   {
      // todo make one array?
      MPI_Allreduce(MPI_IN_PLACE, &outdated_activities, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_lhsrhs, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_nnzs, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_linking_var_bounds, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
   }
   assert(!outdated_activities && !outdated_lhsrhs && !outdated_nnzs && !outdated_linking_var_bounds);
#endif
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

void PresolveData::allreduceLinkingVarBounds()
{
   MPI_Allreduce(MPI_IN_PLACE, &outdated_linking_var_bounds, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);

   if(!outdated_linking_var_bounds)
      return;

   if(distributed)
   {
      StochVector& xlow = dynamic_cast<StochVector&>(*presProb->blx);
      StochVector& xupp = dynamic_cast<StochVector&>(*presProb->bux);
      StochVector& ixlow = dynamic_cast<StochVector&>(*presProb->ixlow);
      StochVector& ixupp = dynamic_cast<StochVector&>(*presProb->ixupp);

      /* allreduce root node bounds */
      SimpleVector& vec_xlow = dynamic_cast<SimpleVector&>(*xlow.vec);
      SimpleVector& vec_ixlow = dynamic_cast<SimpleVector&>(*ixlow.vec);
      SimpleVector& vec_xupp = dynamic_cast<SimpleVector&>(*xupp.vec);
      SimpleVector& vec_ixupp = dynamic_cast<SimpleVector&>(*ixupp.vec);

      MPI_Allreduce(MPI_IN_PLACE, vec_xlow.elements(), vec_xlow.length(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
      MPI_Allreduce(MPI_IN_PLACE, vec_ixlow.elements(), vec_ixlow.length(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
      MPI_Allreduce(MPI_IN_PLACE, vec_xupp.elements(), vec_xupp.length(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
      MPI_Allreduce(MPI_IN_PLACE, vec_ixupp.elements(), vec_ixupp.length(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
   }

   outdated_linking_var_bounds = false;
}

/** allreduces changes in the activities of the linking rows and updates the linking row activities */
void PresolveData::allreduceAndApplyLinkingRowActivities()
{
   MPI_Allreduce(MPI_IN_PLACE, &outdated_activities, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);

   if(!outdated_activities)
      return;

   if(distributed)
   {
      MPI_Allreduce(MPI_IN_PLACE, array_act_chgs, lenght_array_act_chgs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, array_act_unbounded_chgs, lenght_array_act_chgs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }

   dynamic_cast<SimpleVector&>(*actmin_eq_part->vecl).axpy( 1.0, *actmin_eq_chgs);
   dynamic_cast<SimpleVector&>(*actmax_eq_part->vecl).axpy( 1.0, *actmax_eq_chgs);
   dynamic_cast<SimpleVector&>(*actmin_ineq_part->vecl).axpy( 1.0, *actmin_ineq_chgs);
   dynamic_cast<SimpleVector&>(*actmax_ineq_part->vecl).axpy( 1.0, *actmax_ineq_chgs);

   dynamic_cast<SimpleVector&>(*actmin_eq_ubndd->vecl).axpy( 1.0, *actmin_eq_ubndd_chgs);
   dynamic_cast<SimpleVector&>(*actmax_eq_ubndd->vecl).axpy( 1.0, *actmax_eq_ubndd_chgs);
   dynamic_cast<SimpleVector&>(*actmin_ineq_ubndd->vecl).axpy( 1.0, *actmin_ineq_ubndd_chgs);
   dynamic_cast<SimpleVector&>(*actmax_ineq_ubndd->vecl).axpy( 1.0, *actmax_ineq_ubndd_chgs);

   actmin_eq_chgs->setToZero();
   actmax_eq_chgs->setToZero();
   actmin_ineq_chgs->setToZero();
   actmax_ineq_chgs->setToZero();

   actmin_eq_ubndd_chgs->setToZero();
   actmax_eq_ubndd_chgs->setToZero();
   actmin_ineq_ubndd_chgs->setToZero();
   actmax_ineq_ubndd_chgs->setToZero();

   outdated_activities = false;
}

void PresolveData::allreduceAndApplyNnzChanges()
{
   MPI_Allreduce(MPI_IN_PLACE, &outdated_nnzs, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);

   if(!outdated_nnzs)
      return;

   if( distributed )
      MPI_Allreduce(MPI_IN_PLACE, array_nnz_chgs, length_array_nnz_chgs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   /* update local nnzCounters */
   nnzs_col->vec->axpy(-1.0, *nnzs_col_chgs);
   nnzs_row_A->vecl->axpy(-1.0, *nnzs_row_A_chgs);
   nnzs_row_C->vecl->axpy(-1.0, *nnzs_row_C_chgs);

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
   MPI_Allreduce(MPI_IN_PLACE, &outdated_lhsrhs, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);

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

/** Computes minimal and maximal activity of all rows in given matrix. Adds activities to min/max_activities accordingly. */
void PresolveData::addActivityOfBlock( const SparseStorageDynamic& matrix, SimpleVector& min_partact, SimpleVector& unbounded_min, SimpleVector& max_partact,
      SimpleVector& unbounded_max, const SimpleVector& xlow, const SimpleVector& ixlow, const SimpleVector& xupp, const SimpleVector& ixupp) const
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
            // add entry * lower_bound to infRow
            if( ixlow[col] != 0.0)
               min_partact[row] += entry * xlow[col];
            else
               ++unbounded_min[row];
            // add entry * upper_bound to supRow
            if( ixupp[col] != 0.0 )
               max_partact[row] += entry * xupp[col];
            else
               ++unbounded_max[row];
         }
         else
         {
            // add entry * upper_bound to infRow
            if( ixupp[col] != 0.0 )
               min_partact[row] += entry * xupp[col];
            else
               ++unbounded_min[row];
            // add entry * lower_bound to supRow
            if( ixlow[col] != 0.0 )
               max_partact[row] += entry * xlow[col];
            else
               ++unbounded_max[row];
         }
      }
   }
}

void PresolveData::initNnzCounter()
{
   StochGenMatrix& A = dynamic_cast<StochGenMatrix&>(*(presProb->A));
   StochGenMatrix& C = dynamic_cast<StochGenMatrix&>(*(presProb->C));

   StochVectorHandle colClone(dynamic_cast<StochVector*>(nnzs_col->clone()));

   A.getNnzPerRow(*nnzs_row_A);
   C.getNnzPerRow(*nnzs_row_C);
   A.getNnzPerCol(*nnzs_col);
   C.getNnzPerCol(*colClone);

   nnzs_col->axpy(1.0, *colClone);
}

void PresolveData::initSingletons()
{
   /* rows of A */
   /* B0 */
   for(int i = 0; i < nnzs_row_A->vec->n; ++i)
      if( dynamic_cast<SimpleVector&>(*nnzs_row_A->vec)[i] == 1.0)
      {
         singleton_rows.push_back( sROWINDEX(-1, i, EQUALITY_SYSTEM));
      }

   /* Bl0 */
   for(int i = 0; i < nnzs_row_A->vec->n; ++i)
      if(nnzs_row_A->vecl != NULL)
         if( dynamic_cast<SimpleVector&>(*nnzs_row_A->vecl)[i] == 1.0)
            singleton_rows.push_back( sROWINDEX(-2, i, EQUALITY_SYSTEM));

   /* children An + Bn */
   for(int i = 0; i < nChildren; ++i)
      if(nnzs_row_A->children[i]->vec != NULL)
         for(int j = 0; j < nnzs_row_A->children[i]->vec->n; ++j)
            if( dynamic_cast<SimpleVector&>(*nnzs_row_A->children[i]->vec)[j] == 1.0)
               singleton_rows.push_back( sROWINDEX(i, j, EQUALITY_SYSTEM));

   /* rows of C */
   /* B0 */
   for(int i = 0; i < nnzs_row_C->vec->n; ++i)
      if( dynamic_cast<SimpleVector&>(*nnzs_row_C->vec)[i] == 1.0)
         singleton_rows.push_back( sROWINDEX(-1, i, EQUALITY_SYSTEM));

   /* Bl0 */
   for(int i = 0; i < nnzs_row_C->vec->n; ++i)
      if(nnzs_row_C->vecl != NULL)
         if( dynamic_cast<SimpleVector&>(*nnzs_row_C->vecl)[i] == 1.0)
            singleton_rows.push_back( sROWINDEX(-2, i, EQUALITY_SYSTEM));

   /* children An + Bn */
   for(int i = 0; i < nChildren; ++i)
      if(nnzs_row_C->children[i]->vec != NULL)
         for(int j = 0; j < nnzs_row_C->children[i]->vec->n; ++j)
            if( dynamic_cast<SimpleVector&>(*nnzs_row_C->children[i]->vec)[j] == 1.0)
               singleton_rows.push_back( sROWINDEX(i, j, EQUALITY_SYSTEM));

   /* columns */
   for(int i = 0; i < nnzs_col->vec->n; ++i)
      if( dynamic_cast<SimpleVector&>(*nnzs_col->vec)[i] == 1.0)
         singleton_cols.push_back( sCOLINDEX(-1, i));

   for(int i = 0; i < nChildren; ++i)
      if(nnzs_col->children[i]->vec != NULL)
         for(int j = 0; j < nnzs_col->children[i]->vec->n; ++j)
            if( dynamic_cast<SimpleVector&>(*nnzs_col->children[i]->vec)[j] == 1.0)
               singleton_cols.push_back( sCOLINDEX(i, j));
}

bool PresolveData::reductionsEmpty()
{
   if(distributed)
   {
      MPI_Allreduce(MPI_IN_PLACE, &outdated_activities, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_lhsrhs, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_nnzs, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_linking_var_bounds, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
   }
   return !outdated_activities && !outdated_lhsrhs && !outdated_linking_var_bounds && !outdated_nnzs;
}

double PresolveData::addObjOffset(double addOffset)
{
   objOffset += addOffset;
   return objOffset;
}
void PresolveData::setObjOffset(double offset)
{
   objOffset = offset;
}

// todo : if small entry was removed from system no postsolve is necessary - if coefficient was removed because impact of changes in variable are small
// also rhs lhs will be adjusted - this has to be reversed later - there will also be a problem with the reduced costs in than particular row ?
void PresolveData::deleteEntry(SystemType system_type, int node, BlockType block_type, int row_index,
      int& index_k, int& row_end)
{
   assert(-1 <= node && node < nChildren);
   SparseStorageDynamic* storage = getSparseGenMatrix(system_type, node , block_type)->getStorageDynamic();
   assert(0 <= row_index && row_index <= storage->m);

   SimpleVector& xlower = (node == -1) ? dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->blx).vec)
         : dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->blx).children[node]->vec);

   if(block_type == LINKING_CONS_BLOCK || block_type == LINKING_VARS_BLOCK || node == -1)
      outdated_nnzs = true;

   /* adjust rhs and lhs */
   adjustMatrixBoundsBy(system_type, node, block_type, row_index, -storage->M[index_k] * xlower[storage->jcolM[index_k]]);

   /* adjust activity */

   /* adjust nnz counters */
   removeIndexRow(system_type, node, block_type, row_index, 1);
   removeIndexColumn(node, block_type, storage->jcolM[index_k], 1);

   std::swap(storage->M[index_k], storage->M[row_end - 1]);
   std::swap(storage->jcolM[index_k], storage->jcolM[row_end - 1]);
   storage->rowptr[row_index].end--;
   row_end = storage->rowptr[row_index].end;
   index_k--;

   ++elements_deleted;

}

void PresolveData::fixColumn(int node, int col, double value)
{
   assert(-1 <= node && node < nChildren);
   assert(0 <= col);
   postsolver->notifyFixedColumn(node, col, value);

   removeColumn(node, col, value);
}

bool PresolveData::rowPropagatedBounds( SystemType system_type, int node, BlockType block_type, int row, int col, double ubx, double lbx)
{
   assert( -1 <= node && node < nChildren );

   const SimpleVector& ixlow = (node == -1) ? dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->ixlow).vec)
         : dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->ixlow).children[node]->vec);
   const SimpleVector& xlow = (node == -1) ? dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->blx).vec)
         : dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->blx).children[node]->vec);
   const SimpleVector& ixupp = (node == -1) ? dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->ixupp).vec)
         : dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->ixupp).children[node]->vec);
   const SimpleVector& xupp = (node == -1) ? dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->bux).vec)
         : dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->bux).children[node]->vec);

   assert(0 <= col && col < ixlow.n);

   if( ( ixlow[col] != 0.0 && PIPSisLT(ubx, xlow[col]) )
         || (ixupp[col] != 0.0 && PIPSisLT(xupp[col], lbx) )
         || (lbx > ubx))
      abortInfeasible(MPI_COMM_WORLD, "Row Propagation detected infeasible new bounds!", "PresolveData.C", "rowPropagetedBounds");

   bool bounds_changed = false;

   // we do not tighten bounds if impact is too low or bound is bigger than 10e8 // todo : maybe different limit
   // set lower bound
   if( fabs(lbx) < 1e8 && (ixlow[col] == 0.0  || feastol * 1e3 <= fabs(xlow[col] - lbx) ) )
      bounds_changed = bounds_changed || updateUpperBoundVariable(system_type, block_type, node, col, ubx);
   if( fabs(ubx) < 1e8 && (ixupp[col] == 0.0  || feastol * 1e3 <= fabs(xupp[col] - ubx) ) )
      bounds_changed = bounds_changed || updateLowerBoundVariable(system_type, block_type, node, col, lbx);

   /// every process should have the same root node data thus all of them should propagate it's rows similarly
   if( bounds_changed && (block_type == LINKING_VARS_BLOCK || (node == -1 && block_type == LINKING_CONS_BLOCK) ) )
      assert(outdated_linking_var_bounds == true);

   // todo : how to undo propagations from linking constraint rows..
// todo : in case a linking row propagated we'll have to store the whole linking row
//   SparseGenMatrix* mat = getSparseGenMatrix(system_type, node, block_type);
//   assert(row < mat->getStorageDynamic()->m );
//
//   int row_start = mat->getStorageDynamic()->rowptr[row].start;
//   int row_end = mat->getStorageDynamic()->rowptr[row].end;
//
//   assert(row_start < row_end);
// if bounds_changed
//   postsolver->notifyRowPropagated(system_type, node, row, (block_type == LINKING_CONS_BLOCK), col, lbx, ubx, mat->getStorageDynamic()->M + row_start,
//         mat->getStorageDynamic()->jcolM + row_start, row_end - row_start );

   return bounds_changed;
}

// todo : postsolve if bounds adjusted because of deleted matrix entry simply reverse the adjustment - no changes in multipliers - row stays active / inactive ? not true..
void PresolveData::adjustMatrixBoundsBy(SystemType system_type, int node, BlockType block_type, int row_index, double value)
{
   assert( -1 <= node && node <= nChildren);

   if( PIPSisEQ(value, 0.0) )
      return;

   if(block_type == LINKING_CONS_BLOCK || node == -1)
      outdated_lhsrhs = true;

   if(system_type == EQUALITY_SYSTEM)
   {
      if(node == -1)
      {
         (block_type == LINKING_CONS_BLOCK && my_rank == 0) ? (*bound_chgs_A)[row_index] += value :
               dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->bA).vec)[row_index] += value;
      }
      else
      {
         dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->bA).children[node]->vec)[row_index] += value;
      }
   }
   else
   {
      if(node == -1)
      {
         if(block_type == LINKING_CONS_BLOCK && my_rank == 0)
            (*bound_chgs_C)[row_index] += value;
         else
         {
            assert(dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->icupp).vec)[row_index] == 1.0);
            assert(dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->iclow).vec)[row_index] == 1.0);

            dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->bu).vec)[row_index] += value;
            dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->bl).vec)[row_index] += value;
         }
      }
      else
      {
         if( dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->icupp).children[node]->vec)[row_index] == 1.0)
            dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->bu).children[node]->vec)[row_index] += value;

         if( dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->iclow).children[node]->vec)[row_index] == 1.0 )
            dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->bl).children[node]->vec)[row_index] += value;
      }

   }
}

void PresolveData::updateTransposedSubmatrix( SystemType system_type, int node, BlockType block_type, std::vector<std::pair<int, int> >& elements)
{
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
         (system_type == EQUALITY_SYSTEM) ? ((*nnzs_row_A_chgs)[row_index] += amount) : ((*nnzs_row_C_chgs)[row_index] += amount);
         outdated_nnzs = true;
      }
   }
   else
   {
      if(system_type == EQUALITY_SYSTEM)
      {
         (node == -1) ? (dynamic_cast<SimpleVector&>(*nnzs_row_A->vec)[row_index] -= amount) :
               (dynamic_cast<SimpleVector&>(*nnzs_row_A->children[node]->vec)[row_index] -= amount);
         (node == -1) ? assert(0 <= dynamic_cast<SimpleVector&>(*nnzs_row_A->vec)[row_index]) :
               assert(0 <= dynamic_cast<SimpleVector&>(*nnzs_row_A->children[node]->vec)[row_index]);
      }
      else
      {
         (node == -1) ? (dynamic_cast<SimpleVector&>(*nnzs_row_C->vec)[row_index] -= amount) :
               (dynamic_cast<SimpleVector&>(*nnzs_row_C->children[node]->vec)[row_index] -= amount);
         (node == -1) ? assert(0 <= dynamic_cast<SimpleVector&>(*nnzs_row_C->vec)[row_index]) :
               assert(0 <= dynamic_cast<SimpleVector&>(*nnzs_row_C->children[node]->vec)[row_index]);
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
         (*nnzs_col_chgs)[col_index] += amount;
         outdated_nnzs = true;
      }
   }
   else
   {
      dynamic_cast<SimpleVector&>(*nnzs_col->children[node]->vec)[col_index] -= amount;

      assert(0 <= dynamic_cast<SimpleVector&>(*nnzs_col->children[node]->vec)[col_index]);
   }
}

/// removes a column from the whole system A, C by fixing x to fixation
/// updates non-zero counters, rhs, lhs, objective offset and activities
void PresolveData::removeColumn(int node, int col, double fixation)
{
   assert( -1 <= node && node < nChildren );

   if(node == -1)
   {
      removeColumnFromMatrix(EQUALITY_SYSTEM, node, LINKING_VARS_BLOCK, col, fixation);
      removeColumnFromMatrix(INEQUALITY_SYSTEM, node, LINKING_VARS_BLOCK, col, fixation);

      if(hasLinking(EQUALITY_SYSTEM))
         removeColumnFromMatrix(EQUALITY_SYSTEM, node, LINKING_CONS_BLOCK, col, fixation);
      if( hasLinking(INEQUALITY_SYSTEM) )
         removeColumnFromMatrix(INEQUALITY_SYSTEM, node, LINKING_CONS_BLOCK, col, fixation);

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
      double objective_factor = ( ( node == -1 ) ? dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->g).vec)[col]
         : dynamic_cast<SimpleVector&>(*(dynamic_cast<StochVector&>(*presProb->g).children[node]->vec))[col] );
      obj_offset_chgs += objective_factor * fixation;
   }
}

/** remove column - adjust lhs, rhs and activity as well as nnz_counters */
void PresolveData::removeColumnFromMatrix(SystemType system_type, int node, BlockType block_type, int col, double fixation)
{
   bool removed = false;
   SparseGenMatrix* mat = getSparseGenMatrix(system_type, node, block_type);

   SparseStorageDynamic& matrix = mat->getStorageDynamicRef();
   SparseStorageDynamic& matrix_transp = mat->getStorageDynamicTransposedRef();

   assert(0 <= col && col < matrix_transp.m);

   /// remove all entries in column from the sparse storage dynamic
   for( int j = matrix_transp.rowptr[col].start; j < matrix_transp.rowptr[col].end; j++ )
   {
      removed = true;
      const int row = matrix_transp.jcolM[j];
      const double coeff = matrix_transp.M[j];

      removeEntryInDynamicStorage(matrix, row, col);

      removeIndexRow(system_type, node, block_type, row, 1);

      adjustMatrixBoundsBy(system_type, node, block_type, row, - coeff * fixation);
      adjustRowActivityFromDeletion(system_type, node, block_type, row, col, coeff);
   }

   removeIndexColumn(node, block_type, col, matrix_transp.rowptr[col].end - matrix_transp.rowptr[col].start);

   /// update the transposed
   matrix_transp.rowptr[col].end = matrix_transp.rowptr[col].start;
   // todo assert(transposed and normal matrix are in sync)

   if( removed && (node == -1 || block_type == LINKING_CONS_BLOCK) )
   {
      assert(outdated_nnzs == true);
      assert(outdated_activities == true);
      assert(outdated_lhsrhs == true);
   }
}

void PresolveData::removeParallelRow(SystemType system_type, int node, int row, bool linking)
{
   throw std::runtime_error("Not yet implemented");
//   if(postsolver)
//      postsolver->notifyParallelRow()

   removeRow(system_type, node, row, linking);
}
void PresolveData::removeRedundantRow(SystemType system_type, int node, int row, bool linking)
{
   if(postsolver)
      postsolver->notifyRedundantRow(system_type, node, row, linking);

   removeRow(system_type, node, row, linking);
}

void PresolveData::removeRow(SystemType system_type, int node, int row, bool linking)
{
   // todo set activity to zero
   // todo set lhs rhs to zero
   assert(!nodeIsDummy(node, system_type));
   if(linking)
   {
      assert(node == -1);

      /* Bl0 */
      removeRowFromMatrix(system_type, -1, LINKING_CONS_BLOCK, row);

      /* linking rows Bli */
      for(int child = 0; child < nChildren; ++child)
         if(!nodeIsDummy(child, system_type))
            removeRowFromMatrix(system_type, child, LINKING_CONS_BLOCK, row);
   }
   else
   {
      /* Amat */
      removeRowFromMatrix(system_type, node, LINKING_VARS_BLOCK, row);

      /* Bmat */
      if(node != -1)
         removeRowFromMatrix(system_type, node, CHILD_BLOCK, row);
   }
}


// todo : should rhs and lhs be set to zero too? -nah? activity !
void PresolveData::removeRowFromMatrix(SystemType system_type, int node, BlockType block_type, int row)
{
   assert(!nodeIsDummy(node, system_type));
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
      removeIndexColumn(node, block_type, col_idx, 1);
   }

   removeIndexRow(system_type, node, block_type, row, mat_storage->rowptr[row].end - mat_storage->rowptr[row].start);
   mat_storage->rowptr[row].end = mat_storage->rowptr[row].start;
}

void PresolveData::removeEntryInDynamicStorage(SparseStorageDynamic& storage, int rowIdx, int colIdx) const
{
   int i = -1;
   int end = storage.rowptr[rowIdx].end;
   int start = storage.rowptr[rowIdx].start;

   for( i = start; i < end; i++)
   {
      if( storage.jcolM[i] == colIdx )
         break;
   }
   assert( storage.jcolM[i] == colIdx);

   std::swap(storage.M[i], storage.M[end-1]);
   std::swap(storage.jcolM[i], storage.jcolM[end-1]);
   storage.rowptr[rowIdx].end--;
}

/** Verifies if the nnzCounters are still correct. */
// todo make const!
bool PresolveData::verifyNnzcounters()
{
   allreduceAndApplyNnzChanges();

   bool nnzCorrect = true;
   StochVectorHandle nnzColOrig(dynamic_cast<StochVector*>(nnzs_col->cloneFull()));
   StochVectorHandle nnzRowAOrig(dynamic_cast<StochVector*>(nnzs_row_A->cloneFull()));
   StochVectorHandle nnzRowCOrig(dynamic_cast<StochVector*>(nnzs_row_C->cloneFull()));

   nnzs_col->setToZero();
   nnzs_row_A->setToZero();
   nnzs_row_C->setToZero();

   initNnzCounter();

   // linking variables:
   SimpleVector* nColOrigSimple = dynamic_cast<SimpleVector*>(nnzColOrig->vec);
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
      nColOrigSimple = dynamic_cast<SimpleVector*>(nnzColOrig->children[it]->vec);
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
   SimpleVector* nRowAOrigSimple = dynamic_cast<SimpleVector*>(nnzRowAOrig->vec);
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
      nRowAOrigSimple = dynamic_cast<SimpleVector*>(nnzRowAOrig->children[it]->vec);
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
   if(nnzRowAOrig->vecl) // linking rows:
   {
      nRowAOrigSimple = dynamic_cast<SimpleVector*>(nnzRowAOrig->vecl);
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
   SimpleVector* nRowCOrigSimple = dynamic_cast<SimpleVector*>(nnzRowCOrig->vec);
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
      nRowCOrigSimple = dynamic_cast<SimpleVector*>(nnzRowCOrig->children[it]->vec);
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
   if(nnzRowCOrig->vecl) // linking rows:
   {
      nRowCOrigSimple = dynamic_cast<SimpleVector*>(nnzRowCOrig->vecl);
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

SparseGenMatrix* PresolveData::getSparseGenMatrix(SystemType system_type, int node, BlockType block_type)
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

void PresolveData::adjustRowActivityFromDeletion(SystemType system_type, int node, BlockType block_type, int row, int col, double coeff)
{
   assert(-1 <= node && node < nChildren);
   assert(0 <= col);
   assert(0 <= row);
   assert(!PIPSisZero(coeff));

   // get upper and lower bound on variable
   const SimpleVector& ixlow = (node == -1) ? dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->ixlow).vec)
         : dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->ixlow).children[node]->vec);
   const SimpleVector& xlow = (node == -1) ? dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->blx).vec)
         : dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->blx).children[node]->vec);
   const SimpleVector& ixupp = (node == -1) ? dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->ixupp).vec)
         : dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->ixupp).children[node]->vec);
   const SimpleVector& xupp = (node == -1) ? dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->bux).vec)
         : dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*presProb->bux).children[node]->vec);

   StochVector& actmax_part_stoch = (system_type == EQUALITY_SYSTEM) ? *actmax_eq_part : *actmax_ineq_part;
   StochVector& actmin_part_stoch = (system_type == EQUALITY_SYSTEM) ? *actmin_eq_part : *actmin_ineq_part;

   StochVector& actmax_ubndd_stoch = (system_type == EQUALITY_SYSTEM) ? *actmax_eq_ubndd : *actmax_ineq_ubndd;
   StochVector& actmin_ubndd_stoch = (system_type == EQUALITY_SYSTEM) ? *actmin_eq_ubndd : *actmin_ineq_ubndd;

   SimpleVector* actmin_chgs = (system_type == EQUALITY_SYSTEM) ? actmin_eq_chgs : actmin_ineq_chgs;
   SimpleVector* actmax_chgs = (system_type == EQUALITY_SYSTEM) ? actmax_eq_chgs : actmax_ineq_chgs;

   SimpleVector* actmin_ubndd_chgs = (system_type == EQUALITY_SYSTEM) ? actmin_eq_ubndd_chgs : actmin_ineq_ubndd_chgs;
   SimpleVector* actmax_ubndd_chgs = (system_type == EQUALITY_SYSTEM) ? actmax_eq_ubndd_chgs : actmax_ineq_ubndd_chgs;

   SimpleVector* actmin_part;
   SimpleVector* actmax_part;

   SimpleVector* actmin_ubbnd;
   SimpleVector* actmax_ubbnd;

   if( block_type == LINKING_CONS_BLOCK )
   {
      actmin_part = actmin_chgs;
      actmax_part = actmax_chgs;

      actmin_ubbnd = actmin_ubndd_chgs;
      actmax_ubbnd = actmax_ubndd_chgs;
   }
   else
   {
      actmin_part = (node == -1) ? dynamic_cast<SimpleVector*>(actmin_part_stoch.vec) : dynamic_cast<SimpleVector*>(actmin_part_stoch.children[node]->vec);
      actmax_part = (node == -1) ? dynamic_cast<SimpleVector*>(actmax_part_stoch.vec) : dynamic_cast<SimpleVector*>(actmax_part_stoch.children[node]->vec);

      actmin_ubbnd = (node == -1) ? dynamic_cast<SimpleVector*>(actmin_ubndd_stoch.vec) : dynamic_cast<SimpleVector*>(actmin_ubndd_stoch.children[node]->vec);
      actmax_ubbnd = (node == -1) ? dynamic_cast<SimpleVector*>(actmax_ubndd_stoch.vec) : dynamic_cast<SimpleVector*>(actmax_ubndd_stoch.children[node]->vec);
   }

   assert(actmin_part);
   assert(actmax_part);
   assert(actmin_ubbnd);
   assert(actmax_ubbnd);
   assert(col < ixlow.n);
   assert(row < actmin_part->n);
   assert(row < actmax_part->n);
   assert(row < actmin_ubbnd->n);
   assert(row < actmax_ubbnd->n);

   /* depending on the sign of coeff add / substract lbx * coeff/ ubx * coeff from actmin and actmax */
   if( PIPSisLT(0.0, coeff) )
   {
      if(ixupp[col] == 1.0 && xupp[col] != 0.0)
         (*actmax_part)[row] -= coeff*xupp[col];
      else
         --(*actmax_ubbnd)[row];

      if(ixlow[col] == 1.0 && xlow[col] != 0.0)
         (*actmin_part)[row] -= coeff*xlow[col];
      else
         --(*actmin_ubbnd)[row];
   }
   else
   {
      if(ixlow[col] == 1.0 && xlow[col] != 0.0)
         (*actmax_part)[row] -= coeff*xlow[col];
      else
         --(*actmax_ubbnd)[row];

      if(ixupp[col] == 1.0 && xupp[col] != 0.0)
         (*actmin_part)[row] -= coeff*xupp[col];
      else
         --(*actmin_ubbnd)[row];
   }

   if( block_type == LINKING_CONS_BLOCK)
      outdated_activities = true;
}

// todo update activites as well!
bool PresolveData::updateBoundsVariable(SystemType system_type, BlockType block_type, int node, int col, double ubx, double lbx)
{
   SimpleVector& ixlow = (node == -1) ? dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->ixlow).vec)
         : dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->ixlow).children[node]->vec);
   SimpleVector& xlow = (node == -1) ? dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->blx).vec)
         : dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->blx).children[node]->vec);
   SimpleVector& ixupp = (node == -1) ? dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->ixupp).vec)
         : dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->ixupp).children[node]->vec);
   SimpleVector& xupp = (node == -1) ? dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->bux).vec)
         : dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->bux).children[node]->vec);

   if( ( ixlow[col] != 0.0 && PIPSisLT(ubx, xlow[col]) )
         || (ixupp[col] != 0.0 && PIPSisLT(xupp[col], lbx) )
         || (lbx > ubx))
      abortInfeasible(MPI_COMM_WORLD, "Varbounds update detected infeasible new bounds!", "PresolveData.C", "updateBoundsVariable");

   bool updated = false;

   if( ubx < std::numeric_limits<double>::max() && (ixupp[col] == 0.0 || ubx < xupp[col]) )
   {
      updated = true;
      ixupp[col] = 1.0;
      xupp[col] = ubx;
   }

   if( -std::numeric_limits<double>::max() < lbx && (ixlow[col] == 0.0 || xlow[col] < lbx) )
   {
      updated = true;
      ixlow[col] = 1.0;
      xlow[col] = lbx;
   }

   if( updated && (block_type == LINKING_VARS_BLOCK || (node == -1 && block_type == LINKING_CONS_BLOCK) ) )
      outdated_linking_var_bounds = true;

   return updated;
}

//bool PresolveData::isRootInfoInSync()
//{
//
//}
