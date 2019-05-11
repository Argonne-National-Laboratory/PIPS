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

PresolveData::PresolveData(const sData* sorigprob) :
      presProb(sorigprob->cloneFull(true)),
      nRowElemsA(StochVectorHandle(dynamic_cast<StochVector*>(sorigprob->bA->clone()))),
      nRowElemsC(StochVectorHandle(dynamic_cast<StochVector*>(sorigprob->icupp->clone()))),
      nColElems(StochVectorHandle(dynamic_cast<StochVector*>(sorigprob->g->clone()))),
      redRowA(nRowElemsA->clone()),
      redRowC(nRowElemsC->clone()),
      redCol(nColElems->clone()),
      outdated_activities(false),
      actmax_eq(nRowElemsA->clone()),
      actmin_eq(nRowElemsA->clone()),
      actmax_ineq(nRowElemsC->clone()),
      actmin_ineq(nRowElemsC->clone()),
      lenght_array_act_chgs(nRowElemsA->vecl->n * 2 + nRowElemsC->vecl->n * 2),
      array_act_chgs( new double[lenght_array_act_chgs]),
      actmax_eq_chgs(array_act_chgs, nRowElemsA->vecl->n),
      actmin_eq_chgs(array_act_chgs + nRowElemsA->vecl->n, nRowElemsA->vecl->n),
      actmax_ineq_chgs(array_act_chgs + 2 * nRowElemsA->vecl->n, nRowElemsC->vecl->n),
      actmin_ineq_chgs(array_act_chgs + 2 * nRowElemsA->vecl->n + nRowElemsC->vecl->n, nRowElemsC->vecl->n),
      nChildren(nColElems->children.size()),
      objOffset(0.0)
{
   getRankDistributed(MPI_COMM_WORLD, my_rank, distributed);

   // initialize all dynamic transposed sub matrices
   dynamic_cast<StochGenMatrix&>(*presProb->A).initTransposed(true);
   dynamic_cast<StochGenMatrix&>(*presProb->C).initTransposed(true);

   recomputeActivities();
   initNnzCounter();
}

PresolveData::~PresolveData()
{
   delete[] array_act_chgs;
}


sData* PresolveData::finalize()
{
   // this removes all columns and rows that are now empty from the problem
   presProb->cleanUpPresolvedData(*nRowElemsA, *nRowElemsC, *nColElems);

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

   actmax_eq_chgs.setToZero();
   actmin_eq_chgs.setToZero();
   actmax_ineq_chgs.setToZero();
   actmin_ineq_chgs.setToZero();

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
}

/** allreduces changes in the activities of the linking rows and updates the linking row activities */
void PresolveData::updateLinkingRowActivities()
{
   if(!distributed)
      return;

   MPI_Allreduce(MPI_IN_PLACE, array_act_chgs, lenght_array_act_chgs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   SimpleVector& actmin_eq_link = dynamic_cast<SimpleVector&>(*actmin_eq->vecl);
   SimpleVector& actmax_eq_link = dynamic_cast<SimpleVector&>(*actmax_eq->vecl);
   SimpleVector& actmin_ineq_link = dynamic_cast<SimpleVector&>(*actmin_ineq->vecl);
   SimpleVector& actmax_ineq_link = dynamic_cast<SimpleVector&>(*actmax_ineq->vecl);

   actmin_eq_link.axpy( 1.0, actmin_eq_chgs);
   actmax_eq_link.axpy( 1.0, actmax_eq_chgs);
   actmin_ineq_link.axpy( 1.0, actmin_ineq_chgs);
   actmax_ineq_link.axpy( 1.0, actmax_ineq_chgs);

   actmin_eq_chgs.setToZero();
   actmax_eq_chgs.setToZero();
   actmin_ineq_chgs.setToZero();
   actmax_ineq_chgs.setToZero();
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

   StochVectorHandle colClone(dynamic_cast<StochVector*>(nColElems->clone()));

   A.getNnzPerRow(*nRowElemsA);
   C.getNnzPerRow(*nRowElemsC);
   A.getNnzPerCol(*nColElems);
   C.getNnzPerCol(*colClone);

   nColElems->axpy(1.0, *colClone);

#if 0
   int rank = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if( rank == 0 ) std::cout << "write Cols with all cols" << std::endl;
   nColElems->writeToStreamAll(std::cout);

   if( rank == 0 ) std::cout << "write Rows A" << std::endl;
   nRowElemsA->writeToStreamAll(std::cout);

   if( rank == 0 ) std::cout << "write Rows C " << std::endl;
   nRowElemsC->writeToStreamAll(std::cout);
#endif
}

bool PresolveData::combineColAdaptParent()
{
   int myRank, world_size;
   bool iAmDistrib = false;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);

   if( world_size > 1)
      iAmDistrib = true;

   if( iAmDistrib )
   {
      // allgather the length of each colAdaptParent
      const int mylen = getNumberColAdParent();
      int* recvcounts = new int[world_size];

      MPI_Allgather(&mylen, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

      // allgatherv the actual colAdaptParents
      // First, extract the colIdx and val into int* and double* arrays:
      int* colIndicesLocal = new int[mylen];
      double* valuesLocal = new double[mylen];
      for(int i=0; i<mylen; i++)
      {
         colIndicesLocal[i] = getColAdaptParent(i).colIdx;
         valuesLocal[i] = getColAdaptParent(i).val;
      }
      // Second, prepare the receive buffers:
      int lenghtGlobal = recvcounts[0];
      int* displs = new int[world_size];
      displs[0] = 0;
      for(int i = 1; i < world_size; i++)
      {
         lenghtGlobal += recvcounts[i];
         displs[i] = displs[i-1] + recvcounts[i-1];
      }
      int* colIndicesGlobal = new int[lenghtGlobal];
      double* valuesGlobal = new double[lenghtGlobal];
      // Then, do the actual MPI communication:
      MPI_Allgatherv(colIndicesLocal, mylen, MPI_INT, colIndicesGlobal, recvcounts, displs , MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(valuesLocal, mylen, MPI_DOUBLE, valuesGlobal, recvcounts, displs , MPI_DOUBLE, MPI_COMM_WORLD);

      // Reconstruct a colAdaptParent:
      clearColAdaptParent();
      for(int i=0; i<lenghtGlobal; i++)
      {
         COLUMNFORDELETION colWithVal = {colIndicesGlobal[i], valuesGlobal[i]};
         addColToAdaptParent(colWithVal);
      }

      delete[] displs;
      delete[] recvcounts;
      delete[] colIndicesLocal;
      delete[] valuesLocal;
      delete[] colIndicesGlobal;
      delete[] valuesGlobal;
   }

   // Sort colIndicesGlobal (and valuesGlobal accordingly), remove duplicates and find infeasibilities
   std::sort(linkingVariablesMarkedForDeletion.begin(), linkingVariablesMarkedForDeletion.end(), col_is_smaller());
   for(int i = 1; i < getNumberColAdParent(); i++)
      assert( getColAdaptParent(i-1).colIdx <= getColAdaptParent(i).colIdx);

   if(getNumberColAdParent() > 0)
   {
      int colIdxCurrent = getColAdaptParent(0).colIdx;
      double valCurrent = getColAdaptParent(0).val;
      for(int i = 1; i < getNumberColAdParent(); i++)
      {
         if( getColAdaptParent(i).colIdx == colIdxCurrent )
         {
            if( getColAdaptParent(i).val != valCurrent )
            {
               std::cout << "Detected infeasibility (in variable) " << colIdxCurrent << std::endl;
               return false;
            }
            else
            {
               linkingVariablesMarkedForDeletion.erase(linkingVariablesMarkedForDeletion.begin()+i);   //todo: implement more efficiently
               i--;
            }
         }
         else{
            colIdxCurrent = getColAdaptParent(i).colIdx;
            valCurrent = getColAdaptParent(i).val;
         }
      }
   }
   assert( getNumberColAdParent() <= nColElems->vec->n );

   return true;
}

bool PresolveData::reductionsEmpty()
{
   return redRowA->isZero() && redRowC->isZero() && redCol->isZero();
}

void PresolveData::resetRedCounters()
{
   redRowA->setToZero();
   redRowC->setToZero();
   redCol->setToZero();
}

int PresolveData::getNChildren() const
{
   return nChildren;
}

double PresolveData::getObjOffset() const
{
   return objOffset;
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

COLUMNFORDELETION PresolveData::getColAdaptParent(int i) const
{
   assert( i<getNumberColAdParent() );
   return linkingVariablesMarkedForDeletion[i];
}
int PresolveData::getNumberColAdParent() const
{
   return (int)linkingVariablesMarkedForDeletion.size();
}
void PresolveData::addColToAdaptParent(COLUMNFORDELETION colToAdapt)
{
   linkingVariablesMarkedForDeletion.push_back(colToAdapt);
}
void PresolveData::clearColAdaptParent()
{
   linkingVariablesMarkedForDeletion.clear();
}
