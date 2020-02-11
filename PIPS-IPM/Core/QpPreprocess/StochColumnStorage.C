/*
 * StochColumnStorage.C
 *
 *  Created on: 12.12.2019
 *      Author: Nils-Christian Kempke
 *
 *
 *
 *  stored_cols_eq/ineq:
 *      | 0   |
 *      | Bl1 | B1 |
 *      | Bl2 |    | B2 |
 *         .            | B3 |
 *         .                   ...
 *         .
 *      | BlN |                    | BN |
 *      | Bl0 | A1 | A2 | A3 | ... | AN |
 *
 *  b0_block_linking_cols_eq/ineq
 *      | B0  |
 *
 *
 */
// todo improve description

#include "StochColumnStorage.h"
#include "DoubleMatrixTypes.h"
#include "SystemType.h"
#include "StochMatrixUtilities.h"
#include "StochVectorUtilities.h"

StochColumnStorage::StochColumnStorage(const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part) :
 nChildren(matrix_eq_part.children.size() )
{
   assert(matrix_ineq_part.children.size() == nChildren );

   createStorageMatrix(EQUALITY_SYSTEM, matrix_eq_part);
   createStorageMatrix(INEQUALITY_SYSTEM, matrix_ineq_part);
}

StochColumnStorage::~StochColumnStorage()
{
}

void StochColumnStorage::createStorageMatrix(SystemType system_type, const StochGenMatrix& sys_matrix)
{
   SparseGenMatrixHandle& b0_block_storage = (system_type == EQUALITY_SYSTEM) ? B0_eq : B0_ineq;
   StochGenMatrixHandle& col_storage = (system_type == EQUALITY_SYSTEM) ? stored_cols_eq : stored_cols_ineq;

   /* extra storage for b0mat entries in linking variable column we want to store */
   assert(sys_matrix.Bmat);
   b0_block_storage = SparseGenMatrixHandle(sys_matrix.Bmat->cloneEmptyRowsTransposed(true));

   // todo : n, m are wrong here but the counters are broken anyways?
   col_storage = StochGenMatrixHandle(new StochGenMatrix(sys_matrix.id, sys_matrix.n, sys_matrix.m, sys_matrix.mpiComm));


   /* the rest of the matrix is going to get transposed - note that this would noramlly reverse the dummy non-dummy structure of the matrix
    * we will not reverse it here but rather permute the child array such that the dummy non dummy structure stays the same as before
    * this will make the implementation of multiplication etc easier
    * e.g. A1 will not be swapped with Bln but Bl1
    * Bi will stay at their position in the matrix and just get transposed
    */

   /* clone submatrices */
   /* Amat is empty in root node */
   col_storage->Amat = sys_matrix.Amat->cloneEmptyRows(true);
   /* B0mat will not be used and stay empty */
   col_storage->Bmat = sys_matrix.Blmat->cloneEmptyRowsTransposed(true);
   col_storage->Blmat = sys_matrix.Blmat->cloneEmptyRowsTransposed(true);

   for( size_t it = 0; it < sys_matrix.children.size(); it++ )
   {
      const StochGenMatrix& child = *sys_matrix.children[it];

      /* create child */
      StochGenMatrix* child_clone;
      if( child.isKindOf( kStochGenDummyMatrix ) )
         child_clone = new StochGenDummyMatrix(child.id);
      else
         child_clone = new StochGenMatrix(child.id, child.m, child.n, child.mpiComm);

      /* clone submatrices */
      child_clone->Amat = child.Blmat->cloneEmptyRowsTransposed(true);
      child_clone->Bmat = child.Bmat->cloneEmptyRowsTransposed(true);
      child_clone->Blmat = child.Amat->cloneEmptyRowsTransposed(true);

      col_storage->children.push_back(child_clone);
   }
}

int StochColumnStorage::storeCol( int node, int col, const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part)
{
   assert( matrix_eq_part.children.size() == matrix_ineq_part.children.size() );

   if( node == -1 )
      return storeLinkingCol(col, matrix_eq_part, matrix_ineq_part);
   else
      return storeLocalCol(node, col, matrix_eq_part, matrix_ineq_part);
}

int StochColumnStorage::storeLinkingCol(int col, const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part)
{
   assert( matrix_eq_part.children.size() == matrix_ineq_part.children.size() );
   assert( matrix_eq_part.children.size() == stored_cols_eq->children.size() );

   const SparseGenMatrix& B0 = *getSparseGenMatrixFromStochMat(matrix_eq_part, -1, B_MAT);
   const SparseGenMatrix& D0 = *getSparseGenMatrixFromStochMat(matrix_ineq_part, -1, B_MAT);
   const SparseGenMatrix& Bl0 = *getSparseGenMatrixFromStochMat(matrix_eq_part, -1, BL_MAT);
   const SparseGenMatrix& Dl0 = *getSparseGenMatrixFromStochMat(matrix_ineq_part, -1, BL_MAT);

   SparseGenMatrix& Bl0_storage = *getSparseGenMatrixFromStochMat(*stored_cols_eq, -1, BL_MAT);
   SparseGenMatrix& Dl0_storage = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, -1, BL_MAT);

   const int B0_index = B0_eq->appendRow(B0, col);
#ifndef NDEBUG
   const int D0_index = B0_ineq->appendRow(D0, col);
   const int Bl0_index = Bl0_storage.appendRow(Bl0, col);
   const int Dl0_index = Dl0_storage.appendRow(Dl0, col);
   assert(B0_index == D0_index);
   assert(Bl0_index == Dl0_index);
   assert(B0_index == Bl0_index);
#else
   B0_ineq->appendRow(D0, col);
   Bl0_storage.appendRow(Bl0, col);
   Dl0_storage.appendRow(Dl0, col);
#endif

   for(unsigned int i = 0; i < nChildren; ++i)
   {
      const SparseGenMatrix& Ai = *getSparseGenMatrixFromStochMat(matrix_eq_part, i, A_MAT);
      const SparseGenMatrix& Ci = *getSparseGenMatrixFromStochMat(matrix_ineq_part, i, A_MAT);

      SparseGenMatrix& Ai_storage = *getSparseGenMatrixFromStochMat(*stored_cols_eq, i, BL_MAT);
      SparseGenMatrix& Ci_storage = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, i, BL_MAT);

      if( Ai.isKindOf(kStochGenDummyMatrix) )
      {
         assert( Ci.isKindOf(kStochGenDummyMatrix) );
         assert( Ai_storage.isKindOf(kStochGenDummyMatrix) );
         assert( Ci_storage.isKindOf(kStochGenDummyMatrix) );
      }
      else
      {
         assert( !Ci.isKindOf(kStochGenDummyMatrix) );
         assert( !Ai_storage.isKindOf(kStochGenDummyMatrix) );
         assert( !Ci_storage.isKindOf(kStochGenDummyMatrix) );
      }

#ifndef NDEBUG
      const int Ai_index = Ai_storage.appendRow(Ai, col);
      const int Ci_index = Ci_storage.appendRow(Ci, col);
      assert(Ai_index == Ci_index);
      assert(Ai_index == B0_index);
#else
      Ai_storage.appendRow(Ai, col);
      Ci_storage.appendRow(Ci, col);
#endif
   }

   return B0_index;
}

int StochColumnStorage::storeLocalCol(int node, int col, const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part)
{
   assert( 0 <= node && node < static_cast<int>(matrix_eq_part.children.size()) );
   assert( matrix_eq_part.children.size() == matrix_ineq_part.children.size() );
   assert( matrix_eq_part.children.size() == stored_cols_eq->children.size() );

   const SparseGenMatrix& Bi = *getSparseGenMatrixFromStochMat(matrix_eq_part, node, B_MAT);
   const SparseGenMatrix& Di = *getSparseGenMatrixFromStochMat(matrix_ineq_part, node, BL_MAT);
   const SparseGenMatrix& Bli = *getSparseGenMatrixFromStochMat(matrix_eq_part, node, B_MAT);
   const SparseGenMatrix& Dli = *getSparseGenMatrixFromStochMat(matrix_ineq_part, node, BL_MAT);

   SparseGenMatrix& Bi_storage = *getSparseGenMatrixFromStochMat(*stored_cols_eq, node, B_MAT);
   SparseGenMatrix& Di_storage = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, node, A_MAT);
   SparseGenMatrix& Bli_storage = *getSparseGenMatrixFromStochMat(*stored_cols_eq, node, B_MAT);
   SparseGenMatrix& Dli_storage = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, node, A_MAT);

   const int Bi_index = Bi_storage.appendRow(Bi, col);
#ifndef NDEBUG
   const int Di_index = Di_storage.appendRow(Di, col);
   const int Bli_index = Bli_storage.appendRow(Bli, col);
   const int Dli_index = Dli_storage.appendRow(Dli, col);

   assert(Bi_index == Di_index);
   assert(Bli_index == Dli_index);
   assert(Bi_index == Bli_index);
#else
   Di_storage.appendRow(Di, col);
   Bli_storage.appendRow(Bli, col);
   Dli_storage.appendRow(Dli, col);
#endif

   return Bi_index;
}

double StochColumnStorage::multColTimesVec( int node, int col, const StochVector& vec_eq, const StochVector& vec_ineq ) const
{
   assert(-1 <= node && node < static_cast<int>(nChildren));
   assert(nChildren == vec_eq.children.size());
   assert(nChildren == vec_ineq.children.size());

   if(node == -1)
      return multiplyLinkingColTimesVec(col, vec_eq, vec_ineq);
   else
      return multiplyLocalColTimesVec(node, col, vec_eq, vec_ineq);
}

double StochColumnStorage::multColTimesVecWithoutRootNode( int node, int col, const StochVector& vec_eq, const StochVector& vec_ineq ) const
{
   assert(-1 <= node && node < static_cast<int>(nChildren));
   assert(nChildren == vec_eq.children.size());
   assert(nChildren == vec_ineq.children.size());

   if(node == -1)
      return multiplyLinkingColTimesVecWithoutRootNode(col, vec_eq, vec_ineq);
   else
      return multiplyLocalColTimesVec(node, col, vec_eq, vec_ineq);
}

double StochColumnStorage::multiplyLinkingColTimesVec(int col, const StochVector& vec_eq, const StochVector& vec_ineq) const
{
   double res = 0.0;

   /* equality system */
   const SparseGenMatrix& Bl0_mat = *getSparseGenMatrixFromStochMat(*stored_cols_eq, -1, BL_MAT);

   res += B0_eq->localRowTimesVec(getSimpleVecFromRowStochVec(vec_eq, -1, false), col);
   res += Bl0_mat.localRowTimesVec(getSimpleVecFromRowStochVec(vec_eq, -1, true), col);

   /* inequality system */
   const SparseGenMatrix& Dl0_mat = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, -1, BL_MAT);

   res += B0_ineq->localRowTimesVec(getSimpleVecFromRowStochVec(vec_ineq, -1, false), col);
   res += Dl0_mat.localRowTimesVec(getSimpleVecFromRowStochVec(vec_ineq, -1, true), col);

   /* Amat equality and inequality */
   for(unsigned int i = 0; i < nChildren; ++i)
   {
      if( !vec_eq.children[i]->isKindOf(kStochDummy) )
      {
         assert(!vec_ineq.children[i]->isKindOf(kStochDummy));

         const SparseGenMatrix& A_mat = *getSparseGenMatrixFromStochMat(*stored_cols_eq, -1, BL_MAT);
         assert(vec_eq.children[i]->vec);

         const SimpleVector& a_vec = dynamic_cast<const SimpleVector&>(*vec_eq.children[i]->vec);
         res += A_mat.localRowTimesVec(a_vec, col);

         const SparseGenMatrix& C_mat = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, -1, BL_MAT);
         assert(vec_ineq.children[i]->vec);

         const SimpleVector& c_vec = dynamic_cast<const SimpleVector&>(*vec_ineq.children[i]->vec);
         res += C_mat.localRowTimesVec(c_vec, col);
      }
   }

   return res;
}

double StochColumnStorage::multiplyLinkingColTimesVecWithoutRootNode(int col, const StochVector& vec_eq, const StochVector& vec_ineq) const
{
   double res = 0.0;

   /* Amat equality and inequality */
   for(unsigned int i = 0; i < nChildren; ++i)
   {
      if( !vec_eq.children[i]->isKindOf(kStochDummy) )
      {
         assert(!vec_ineq.children[i]->isKindOf(kStochDummy));

         const SparseGenMatrix& A_mat = *getSparseGenMatrixFromStochMat(*stored_cols_eq, -1, BL_MAT);
         assert(vec_eq.children[i]->vec);

         const SimpleVector& a_vec = dynamic_cast<const SimpleVector&>(*vec_eq.children[i]->vec);
         res += A_mat.localRowTimesVec(a_vec, col);

         const SparseGenMatrix& C_mat = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, -1, BL_MAT);
         assert(vec_ineq.children[i]->vec);

         const SimpleVector& c_vec = dynamic_cast<const SimpleVector&>(*vec_ineq.children[i]->vec);
         res += C_mat.localRowTimesVec(c_vec, col);
      }
   }

   return res;
}

double StochColumnStorage::multiplyLocalColTimesVec(int node, int col, const StochVector& vec_eq, const StochVector& vec_ineq) const
{
   double res = 0.0;
   assert(!vec_eq.isKindOf(kStochDummy) && !vec_ineq.isKindOf(kStochDummy));
   /* equality system */
   const SparseGenMatrix& Bi_mat = *getSparseGenMatrixFromStochMat(*stored_cols_eq, node, B_MAT);
   const SparseGenMatrix& Bli_mat = *getSparseGenMatrixFromStochMat(*stored_cols_eq, node, A_MAT);

   res += Bi_mat.localRowTimesVec(getSimpleVecFromRowStochVec(vec_eq, node, false), node);
   res += Bli_mat.localRowTimesVec(getSimpleVecFromRowStochVec(vec_eq, -1, true), node);

   /* inequality system */
   const SparseGenMatrix& Di_mat = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, node, B_MAT);
   const SparseGenMatrix& Dli_mat = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, node, A_MAT);

   res += Di_mat.localRowTimesVec(getSimpleVecFromRowStochVec(vec_ineq, node, false), node);
   res += Dli_mat.localRowTimesVec(getSimpleVecFromRowStochVec(vec_ineq, -1, true), node);

   return res;
}

