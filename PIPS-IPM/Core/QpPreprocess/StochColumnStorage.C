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

#include "StochColumnStorage.h"
#include "DoubleMatrixTypes.h"

StochColumnStorage::StochColumnStorage(const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part)
{
   createStorageMatrix(b0_block_linking_cols_eq, stored_cols_eq, matrix_eq_part);
   createStorageMatrix(b0_block_linking_cols_ineq, stored_cols_ineq, matrix_ineq_part);
}

StochColumnStorage::~StochColumnStorage()
{
}

void StochColumnStorage::createStorageMatrix(SparseGenMatrix* b0_block_storage, StochGenMatrix* col_storage, const StochGenMatrix& sys_matrix)
{
   /* extra storage for b0mat entries in linking variable column we want to store */
   assert(sys_matrix.Bmat);
   b0_block_storage = sys_matrix.Bmat->cloneEmptyRowsTransposed(true);

   // todo : n, m are wrong here but the counters are broken anyways?
   col_storage = new StochGenMatrix(sys_matrix.id, sys_matrix.n, sys_matrix.m, sys_matrix.mpiComm);


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
   if( node == -1 )
      return storeLinkingCol(col, matrix_eq_part, matrix_ineq_part);
   else
      return storeLocalCol(node, col, matrix_eq_part, matrix_ineq_part);
}

int StochColumnStorage::storeLinkingCol(int col, const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part)
{

   return 0;
}

int StochColumnStorage::storeLocalCol(int node, int col, const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part)
{

   return 0;
}

double StochColumnStorage::multColTimesVec( int node, int col, const StochVector& vec_eq, const StochVector& vec_ineq ) const
{
   return 0;
}

double StochColumnStorage::getColCoefficientAtRow( int node, int col, int row) const
{
   return 0;
}
