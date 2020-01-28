/*
 * StochRowStorage.C
 *
 *  Created on: 12.12.2019
 *      Author: Nils-Christian Kempke
 */

#include "StochRowStorage.h"
#include "StochMatrixUtilities.h"

StochRowStorage::StochRowStorage(const StochGenMatrix& system_matrix) :
   row_storage(system_matrix.cloneEmptyRows(true))
{
}


StochRowStorage::~StochRowStorage()
{
}

int StochRowStorage::storeRow( int node, int row, bool linking_row, const StochGenMatrix& matrix_row)
{
   const int stored_row_index = row_storage->appendRow(matrix_row, node, row, linking_row);
   return stored_row_index;
}

/* TODO : it would probably be nice to have something like a StochVectorView that one can create at this point and then use like a normal StochVec.. */
/** y = beta * y + alpha * stored */
void StochRowStorage::axpyAtRow(double beta, StochVector &y, double alpha, int node, int row, bool linking_row) const
{
   y.scale(beta);

   row_storage->axpyWithRowAt(beta, y, node, row, linking_row);
}

double StochRowStorage::multRowTimesVec( int node, int row, bool linking_row, const StochVector& vec ) const
{
   const double res = row_storage->localRowTimesVec(vec, node, row, linking_row);
   return res;
}

double StochRowStorage::multLinkingRowTimesVecWithoutBl0( int row, const StochVector& vec) const
{
   const double res_full = row_storage->localRowTimesVec(vec, -1, row, true);
   const double res_bl0 = row_storage->Blmat->localRowTimesVec(dynamic_cast<const SimpleVector&>(*vec.vec), row);

   return res_full - res_bl0;
}

double StochRowStorage::getRowCoefficientAtColumn( int node_row, int row, bool linking_row, int node_column, int column ) const
{
   BlockType block_col = B_MAT;
   if( linking_row )
      block_col = BL_MAT;
   if( node_column == -1 && node_row != 1)
      block_col = A_MAT;

   const SparseStorageDynamic& mat = *getSparseGenMatrixFromStochMat(*row_storage, node_row, block_col )->getStorageDynamic();

   const int row_start = mat.getRowPtr(row).start;
   const int row_end = mat.getRowPtr(row).end;

   for(int i = row_start; i < row_end; ++i)
   {
      if( mat.getJcolM(i) == column )
         return mat.getMat(i);
   }
   return 0.0;
}
