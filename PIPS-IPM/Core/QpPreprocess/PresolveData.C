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
#include <limits>

//      SimpleVectorHandle rhs_A_chgs;
//      SimpleVectorHandle lhs_C_chgs;
//      SimpleVectorHandle rhs_C_chgs;

PresolveData::PresolveData(const sData* sorigprob) :
      outdated_activities(false),
      outdated_bounds(false),
      outdated_nnzs(false),
      nnzs_row_A(dynamic_cast<StochVector*>(sorigprob->bA->clone())),
      nnzs_row_C(dynamic_cast<StochVector*>(sorigprob->icupp->clone())),
      nnzs_col(dynamic_cast<StochVector*>(sorigprob->g->clone())),
      actmax_eq(nnzs_row_A->clone()),
      actmin_eq(nnzs_row_A->clone()),
      actmax_ineq(nnzs_row_C->clone()),
      actmin_ineq(nnzs_row_C->clone()),
      nChildren(nnzs_col->children.size()),
      objOffset(0.0), obj_offset_chgs(0.0),
      elements_deleted(0), elements_deleted_transposed(0)
{
   presProb = sorigprob->cloneFull(true);

   int n_linking_vars = (nnzs_col->vec) ? nnzs_col->vec->n : 0;
   nnzs_col_chgs = SimpleVectorHandle(new SimpleVector(n_linking_vars));

   long long n_linking_A = (nnzs_row_A->vecl) ? nnzs_row_A->vecl->n : 0;
   long long n_linking_C = (nnzs_row_C->vecl) ? nnzs_row_C->vecl->n : 0;
   assert( 2 * n_linking_A + 2 * n_linking_C < std::numeric_limits<int>::max() );

   nnzs_row_A_chgs = SimpleVectorHandle(new SimpleVector(n_linking_A));
   nnzs_row_C_chgs = SimpleVectorHandle(new SimpleVector(n_linking_C));

   lenght_array_act_chgs = n_linking_A * 2 + n_linking_C * 2;
   array_act_chgs = new double[lenght_array_act_chgs]();

   actmax_eq_chgs = SimpleVectorHandle( new SimpleVector(array_act_chgs, n_linking_A));
   actmin_eq_chgs = SimpleVectorHandle( new SimpleVector(array_act_chgs + n_linking_A, n_linking_A));
   actmax_ineq_chgs = SimpleVectorHandle( new SimpleVector(array_act_chgs + 2 * n_linking_A, n_linking_C));
   actmin_ineq_chgs = SimpleVectorHandle( new SimpleVector(array_act_chgs + 2 * n_linking_A + n_linking_C, n_linking_C));

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
}

PresolveData::~PresolveData()
{
   delete[] array_act_chgs;
   delete[] array_bound_chgs;
}


sData* PresolveData::finalize()
{

#ifndef NDEBUG
   if(distributed)
   {
      MPI_Allreduce(MPI_IN_PLACE, &outdated_activities, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_bounds, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &outdated_nnzs, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
   }
   assert(!outdated_activities && !outdated_bounds && !outdated_nnzs);
#endif

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
   if( distributed )
      MPI_Allreduce(MPI_IN_PLACE, &outdated_activities, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);

   if(!outdated_activities)
      return;

   const StochGenMatrix& mat_A = dynamic_cast<const StochGenMatrix&>(*presProb->A);
   const StochGenMatrix& mat_C = dynamic_cast<const StochGenMatrix&>(*presProb->C);

   const StochVector& xupp = dynamic_cast<StochVector&>(*presProb->bux);
   const StochVector& ixupp = dynamic_cast<StochVector&>(*presProb->ixupp);
   const StochVector& xlow = dynamic_cast<StochVector&>(*presProb->blx);
   const StochVector& ixlow = dynamic_cast<StochVector&>(*presProb->ixlow);

   /* reset vectors keeping track of activities */
   if(!linking_only)
   {
      actmin_eq->setToZero();
      actmax_eq->setToZero();
      actmin_ineq->setToZero();
      actmax_ineq->setToZero();
   }
   else
   {
      actmin_eq->vecl->setToZero();
      actmax_eq->vecl->setToZero();
      actmin_ineq->vecl->setToZero();
      actmax_ineq->vecl->setToZero();
   }

   actmax_eq_chgs->setToZero();
   actmin_eq_chgs->setToZero();
   actmax_ineq_chgs->setToZero();
   actmin_ineq_chgs->setToZero();

   /* compute activities at root node */
   const SimpleVector& xupp_root = dynamic_cast<const SimpleVector&>(*xupp.vec);
   const SimpleVector& ixupp_root = dynamic_cast<const SimpleVector&>(*ixupp.vec);
   const SimpleVector& xlow_root = dynamic_cast<const SimpleVector&>(*xlow.vec);
   const SimpleVector& ixlow_root = dynamic_cast<const SimpleVector&>(*ixlow.vec);

   /* A0/B0 */
   if(!linking_only)
   {
      SimpleVector& actmin_eq_root = dynamic_cast<SimpleVector&>(*actmin_eq->vec);
      SimpleVector& actmax_eq_root = dynamic_cast<SimpleVector&>(*actmax_eq->vec);
      SimpleVector& actmin_ineq_root = dynamic_cast<SimpleVector&>(*actmin_ineq->vec);
      SimpleVector& actmax_ineq_root = dynamic_cast<SimpleVector&>(*actmax_ineq->vec);

      addActivityOfBlock(mat_A.Bmat->getStorageDynamicRef(), actmin_eq_root, actmax_eq_root, xlow_root, ixlow_root, xupp_root, ixupp_root);

      addActivityOfBlock(mat_C.Bmat->getStorageDynamicRef(), actmin_ineq_root, actmax_ineq_root, xlow_root, ixlow_root, xupp_root, ixupp_root);
   }

   SimpleVector& actmin_eq_link = dynamic_cast<SimpleVector&>(*actmin_eq->vecl);
   SimpleVector& actmax_eq_link = dynamic_cast<SimpleVector&>(*actmax_eq->vecl);
   SimpleVector& actmin_ineq_link = dynamic_cast<SimpleVector&>(*actmin_ineq->vecl);
   SimpleVector& actmax_ineq_link = dynamic_cast<SimpleVector&>(*actmax_ineq->vecl);
   /* Bl0 */
   if(my_rank == 0)
   {
      addActivityOfBlock(mat_A.Blmat->getStorageDynamicRef(), actmin_eq_link, actmax_eq_link, xlow_root, ixlow_root, xupp_root, ixupp_root);

      addActivityOfBlock(mat_C.Blmat->getStorageDynamicRef(), actmin_ineq_link, actmax_ineq_link, xlow_root, ixlow_root, xupp_root, ixupp_root);
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
            SimpleVector& actmin_eq_child = dynamic_cast<SimpleVector&>(*actmin_eq->children[node]->vec);
            SimpleVector& actmax_eq_child = dynamic_cast<SimpleVector&>(*actmax_eq->children[node]->vec);

            /* Ai */
            addActivityOfBlock(mat_A.children[node]->Amat->getStorageDynamicRef(), actmin_eq_child, actmax_eq_child, xlow_root, ixlow_root, xupp_root, ixupp_root);

            /* Bi */
            addActivityOfBlock(mat_A.children[node]->Bmat->getStorageDynamicRef(), actmin_eq_child, actmax_eq_child, xlow_child, ixlow_child, xupp_child, ixupp_child);
         }

         /* Bli */
         addActivityOfBlock(mat_A.children[node]->Blmat->getStorageDynamicRef(), actmin_eq_link, actmax_eq_link, xlow_child, ixlow_child, xupp_child, ixupp_child);


      }
      /* INEQUALITY_SYSTEM */
      if( !mat_C.children[node]->isKindOf(kStochGenDummyMatrix) )
      {
         if( !linking_only )
         {
            SimpleVector& actmin_ineq_child = dynamic_cast<SimpleVector&>(*actmin_ineq->children[node]->vec);
            SimpleVector& actmax_ineq_child = dynamic_cast<SimpleVector&>(*actmax_ineq->children[node]->vec);

            /* Ai */
            addActivityOfBlock(mat_C.children[node]->Amat->getStorageDynamicRef(), actmin_ineq_child, actmax_ineq_child, xlow_root, ixlow_root, xupp_root, ixupp_root);

            /* Bi */
            addActivityOfBlock(mat_C.children[node]->Bmat->getStorageDynamicRef(), actmin_ineq_child, actmax_ineq_child, xlow_child, ixlow_child, xupp_child, ixupp_child);
         }

         /* Bli */
         addActivityOfBlock(mat_C.children[node]->Blmat->getStorageDynamicRef(), actmin_ineq_link, actmax_ineq_link, xlow_child, ixlow_child, xupp_child, ixupp_child);
      }

   }

   /* allreduce linking constraint activities */
   if( distributed )
   {
      // todo is copying and then allreducing once cheaper than allreducing 4 times ?

      MPI_Allreduce(MPI_IN_PLACE, dynamic_cast<SimpleVector*>(actmin_eq->vecl)->elements(), actmin_eq->vecl->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, dynamic_cast<SimpleVector*>(actmax_eq->vecl)->elements(), actmax_eq->vecl->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, dynamic_cast<SimpleVector*>(actmin_ineq->vecl)->elements(), actmin_ineq->vecl->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, dynamic_cast<SimpleVector*>(actmax_ineq->vecl)->elements(), actmax_ineq->vecl->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }

   outdated_activities = false;
}

/** allreduces changes in the activities of the linking rows and updates the linking row activities */
void PresolveData::allreduceAndApplyLinkingRowActivities()
{
   if(distributed)
   {
      MPI_Allreduce(MPI_IN_PLACE, array_act_chgs, lenght_array_act_chgs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }

   dynamic_cast<SimpleVector&>(*actmin_eq->vecl).axpy( 1.0, *actmin_eq_chgs);
   dynamic_cast<SimpleVector&>(*actmax_eq->vecl).axpy( 1.0, *actmax_eq_chgs);
   dynamic_cast<SimpleVector&>(*actmin_ineq->vecl).axpy( 1.0, *actmin_ineq_chgs);
   dynamic_cast<SimpleVector&>(*actmax_ineq->vecl).axpy( 1.0, *actmax_ineq_chgs);

   actmin_eq_chgs->setToZero();
   actmax_eq_chgs->setToZero();
   actmin_ineq_chgs->setToZero();
   actmax_ineq_chgs->setToZero();
}

/** Computes minimal and maximal activity of all rows in given matrix. Adds activities to min/max_activities accordingly. */
void PresolveData::addActivityOfBlock( const SparseStorageDynamic& matrix, SimpleVector& min_activities, SimpleVector& max_activities,
      const SimpleVector& xlow, const SimpleVector& ixlow, const SimpleVector& xupp, const SimpleVector& ixupp) const
{
   assert( xlow.n == matrix.n && ixlow.n == matrix.n && xupp.n == matrix.n && ixupp.n == matrix.n );
   assert( max_activities.n == matrix.m && min_activities.n == matrix.m);

   for( int row = 0; row < matrix.m; ++row)
   {
      for( int j = matrix.rowptr[row].start; j < matrix.rowptr[row].end; j++)
      {
         const int col = matrix.jcolM[j];
         const double entry = matrix.M[j];

         if( entry > 0)
         {
            // add entry * lower_bound to infRow
            if( ixlow[col] != 0.0)
               min_activities[row] += entry * xlow[col];
            else
               min_activities[row] = -std::numeric_limits<double>::infinity();
            // add entry * upper_bound to supRow
            if( ixupp[col] != 0.0 )
               max_activities[row] += entry * xupp[col];
            else
               max_activities[row] = std::numeric_limits<double>::infinity();
         }
         else
         {
            // add entry * upper_bound to infRow
            if( ixupp[col] != 0.0 )
               min_activities[row] += entry * xupp[col];
            else
               min_activities[row] = -std::numeric_limits<double>::infinity();
            // add entry * lower_bound to supRow
            if( ixlow[col] != 0.0 )
               max_activities[row] += entry * xlow[col];
            else
               max_activities[row] = std::numeric_limits<double>::infinity();
         }

         /* stop computing the activity of a row if it is already +/-infinity */
         if( max_activities[row] == std::numeric_limits<double>::infinity() && min_activities[row] == -std::numeric_limits<double>::infinity() )
            continue;
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

//bool PresolveData::combineColAdaptParent()
//{
//   int myRank, world_size;
//   bool iAmDistrib = false;
//   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
//   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//
//   if( world_size > 1)
//      iAmDistrib = true;
//
//   if( iAmDistrib )
//   {
//      // allgather the length of each colAdaptParent
//      const int mylen = getNumberColAdParent();
//      int* recvcounts = new int[world_size];
//
//      MPI_Allgather(&mylen, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
//
//      // allgatherv the actual colAdaptParents
//      // First, extract the colIdx and val into int* and double* arrays:
//      int* colIndicesLocal = new int[mylen];
//      double* valuesLocal = new double[mylen];
//      for(int i=0; i<mylen; i++)
//      {
//         colIndicesLocal[i] = getColAdaptParent(i).colIdx;
//         valuesLocal[i] = getColAdaptParent(i).val;
//      }
//      // Second, prepare the receive buffers:
//      int lenghtGlobal = recvcounts[0];
//      int* displs = new int[world_size];
//      displs[0] = 0;
//      for(int i = 1; i < world_size; i++)
//      {
//         lenghtGlobal += recvcounts[i];
//         displs[i] = displs[i-1] + recvcounts[i-1];
//      }
//      int* colIndicesGlobal = new int[lenghtGlobal];
//      double* valuesGlobal = new double[lenghtGlobal];
//      // Then, do the actual MPI communication:
//      MPI_Allgatherv(colIndicesLocal, mylen, MPI_INT, colIndicesGlobal, recvcounts, displs , MPI_INT, MPI_COMM_WORLD);
//      MPI_Allgatherv(valuesLocal, mylen, MPI_DOUBLE, valuesGlobal, recvcounts, displs , MPI_DOUBLE, MPI_COMM_WORLD);
//
//      // Reconstruct a colAdaptParent:
//      clearColAdaptParent();
//      for(int i=0; i<lenghtGlobal; i++)
//      {
//         COLUMNFORDELETION colWithVal = {colIndicesGlobal[i], valuesGlobal[i]};
//         addColToAdaptParent(colWithVal);
//      }
//
//      delete[] displs;
//      delete[] recvcounts;
//      delete[] colIndicesLocal;
//      delete[] valuesLocal;
//      delete[] colIndicesGlobal;
//      delete[] valuesGlobal;
//   }
//
//   // Sort colIndicesGlobal (and valuesGlobal accordingly), remove duplicates and find infeasibilities
//   std::sort(linkingVariablesMarkedForDeletion.begin(), linkingVariablesMarkedForDeletion.end(), col_is_smaller());
//   for(int i = 1; i < getNumberColAdParent(); i++)
//      assert( getColAdaptParent(i-1).colIdx <= getColAdaptParent(i).colIdx);
//
//   if(getNumberColAdParent() > 0)
//   {
//      int colIdxCurrent = getColAdaptParent(0).colIdx;
//      double valCurrent = getColAdaptParent(0).val;
//      for(int i = 1; i < getNumberColAdParent(); i++)
//      {
//         if( getColAdaptParent(i).colIdx == colIdxCurrent )
//         {
//            if( getColAdaptParent(i).val != valCurrent )
//            {
//               std::cout << "Detected infeasibility (in variable) " << colIdxCurrent << std::endl;
//               return false;
//            }
//            else
//            {
//               linkingVariablesMarkedForDeletion.erase(linkingVariablesMarkedForDeletion.begin()+i);   //todo: implement more efficiently
//               i--;
//            }
//         }
//         else{
//            colIdxCurrent = getColAdaptParent(i).colIdx;
//            valCurrent = getColAdaptParent(i).val;
//         }
//      }
//   }
//   assert( getNumberColAdParent() <= nnzs_col->vec->n );
//
//   return true;
//}

bool PresolveData::reductionsEmpty()
{
   return nnzs_row_A_chgs->isZero() && nnzs_row_C_chgs->isZero() && nnzs_col_chgs->isZero();
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

//COLUMNFORDELETION PresolveData::getColAdaptParent(int i) const
//{
//   assert( i<getNumberColAdParent() );
//   return linkingVariablesMarkedForDeletion[i];
//}
//int PresolveData::getNumberColAdParent() const
//{
//   return (int)linkingVariablesMarkedForDeletion.size();
//}
//void PresolveData::addColToAdaptParent(COLUMNFORDELETION colToAdapt)
//{
//   linkingVariablesMarkedForDeletion.push_back(colToAdapt);
//}
//void PresolveData::clearColAdaptParent()
//{
//   linkingVariablesMarkedForDeletion.clear();
//}

// todo make one vector?
void PresolveData::allreduceAndApplyNnzChanges()
{
   // todo check if necessary
   if( distributed )
   {
      /* allreduce the linking variable columns */
      MPI_Allreduce(MPI_IN_PLACE, nnzs_col_chgs->elements(), nnzs_col_chgs->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

      /* allreduce the linking cons rows */
      MPI_Allreduce(MPI_IN_PLACE, nnzs_row_A_chgs->elements(), nnzs_row_A_chgs->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      MPI_Allreduce(MPI_IN_PLACE, nnzs_row_C_chgs->elements(), nnzs_row_C_chgs->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   }

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
}

void PresolveData::allreduceAndApplyBoundChanges()
{
   if(distributed)
   {
      MPI_Allreduce(MPI_IN_PLACE, array_bound_chgs, lenght_array_bound_chgs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }

   dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*presProb->bA).vecl)->axpy( 1.0, *bound_chgs_A);
   dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*presProb->bl).vecl)->axpy( 1.0, *bound_chgs_C);
   dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*presProb->bu).vecl)->axpy( 1.0, *bound_chgs_C);

   bound_chgs_A->setToZero();
   bound_chgs_C->setToZero();
}

// todo : no postsolve required?
void PresolveData::deleteEntry(SystemType system_type, int node, BlockType block_type, SparseStorageDynamic* storage, int row_index,
      int& index_k, int& row_end)
{
   /* adjust nnz counters */
   removeIndexRow(system_type, node, block_type, row_index);
   removeIndexColumn(node, block_type, storage->jcolM[index_k]);

   std::swap(storage->M[index_k], storage->M[row_end - 1]);
   std::swap(storage->jcolM[index_k], storage->jcolM[row_end - 1]);
   storage->rowptr[row_index].end--;
   row_end = storage->rowptr[row_index].end;
   index_k--;

   ++elements_deleted;
}

// todo : no postsolve necessary?
void PresolveData::adjustMatrixBoundsBy(SystemType system_type, int node, BlockType block_type, int row_index, double value)
{
   assert( -1 <= node && node <= nChildren);

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
         assert(dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->icupp).children[node]->vec)[row_index] == 1.0);
         assert(dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->iclow).children[node]->vec)[row_index] == 1.0);

         dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->bu).children[node]->vec)[row_index] += value;
         dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*presProb->bl).children[node]->vec)[row_index] += value;
      }

   }
}

void PresolveData::updateTransposedSubmatrix( SparseStorageDynamic* transposed, std::vector<std::pair<int, int> >& elements)
{
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
}

void PresolveData::removeIndexRow(SystemType system_type, int node, BlockType block_type, int row_index)
{
   assert(-1 <= node && node <= nChildren);

   /* linking constraints get stored */
   if(block_type == LINKING_CONS_BLOCK)
   {
      if(my_rank == 0 || node != -1)
         (system_type == EQUALITY_SYSTEM) ? ++(*nnzs_row_A_chgs)[row_index] : ++(*nnzs_row_C_chgs)[row_index];
   }
   else
   {
      if(system_type == EQUALITY_SYSTEM)
      {
         (node == -1) ? --dynamic_cast<SimpleVector&>(*nnzs_row_A->vec)[row_index] : --dynamic_cast<SimpleVector&>(*nnzs_row_A->children[node]->vec)[row_index];
         assert(0 <= dynamic_cast<SimpleVector&>(*nnzs_row_A->vec)[row_index]);
         assert(0 <= dynamic_cast<SimpleVector&>(*nnzs_row_A->children[node]->vec)[row_index]);
      }
      else
      {
         (node == -1) ? --dynamic_cast<SimpleVector&>(*nnzs_row_C->vec)[row_index] : --dynamic_cast<SimpleVector&>(*nnzs_row_C->children[node]->vec)[row_index];
         assert(0 <= dynamic_cast<SimpleVector&>(*nnzs_row_C->vec)[row_index]);
         assert(0 <= dynamic_cast<SimpleVector&>(*nnzs_row_C->children[node]->vec)[row_index]);
      }
   }
}

void PresolveData::removeIndexColumn(int node, BlockType block_type, int col_index)
{
   assert(-1 <= node && node <= nChildren);

   /* linking constraints get stored */
   if(node == -1 || block_type == LINKING_VARS_BLOCK)
   {
      if(my_rank == 0 || node != -1)
         ++(*nnzs_col_chgs)[col_index];
   }
   else
   {
      --dynamic_cast<SimpleVector&>(*nnzs_col->children[node]->vec)[col_index];
      assert(0 <= dynamic_cast<SimpleVector&>(*nnzs_col->children[node]->vec)[col_index]);
   }
}

void PresolveData::deleteColumn()
{

}
void PresolveData::deleteRow(SystemType system_type, int node, int idx, bool linking)
{

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

void PresolveData::getStorageDynamic(SystemType system_type, int node, BlockType block_type, SparseGenMatrix* mat)
{
//   if(system_type == EQUALITY_SYSTEM)
//   {
//      if(node == -1)
//      {
//         mat = (block_type == LINKING_CONS_BLOCK) ? dynamic_cast<StochGenMatrixHandle>(presProb->A)->Blmat :
//               dynamic_cast<StochGenMatrixHandle>(presProb->A)->Bmat;
//      }
//      else
//      {
//         mat = ()
//      }
//   }
//   else
//   {
//      if(node == -1)
//      {
//
//      }
//      else
//      {
//
//      }
//   }
}

