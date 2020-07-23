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

int StochRowStorage::storeRow( const INDEX& row, const StochGenMatrix& matrix_row)
{
   assert(row.isRow());
   const int stored_row_index = row_storage->appendRow(matrix_row, row.getNode(), row.getIndex(), row.getLinking());
   return stored_row_index;
}

/* TODO : it would probably be nice to have something like a StochVectorView that one can create at this point and then use like a normal StochVec.. */
/** y = beta * y + alpha * stored */
/** if y_linking is defined the result of beta * y + alpha * stored for linking variables will be put there */
void StochRowStorage::axpyAtRow(double beta, StochVector* y, SimpleVector* y_linking, double alpha, const INDEX& row) const
{
   assert( row.isRow() );
   assert( y );
   if( !PIPSisEQ(beta, 1.0) )
   {
      y->scale(beta);
      if( y_linking )
         y_linking->scale(beta);
   }

   row_storage->axpyWithRowAt(alpha, y, y_linking, row.getNode(), row.getIndex(), row.getLinking());
}

void StochRowStorage::axpyAtRowPosNeg(double beta, StochVector* y_pos, SimpleVector* y_link_pos, StochVector* y_neg, SimpleVector* y_link_neg, double alpha, const INDEX& row) const
{
   assert( row.isRow() );
   assert( (y_link_pos && y_link_neg) || (!y_link_pos && !y_link_neg) );

   if( !PIPSisEQ(beta, 1.0) )
   {
      y_pos->scale(beta);
      y_neg->scale(beta);

      if( y_link_pos )
      {
         y_link_pos->scale(beta);
         y_link_neg->scale(beta);
      }
   }

   row_storage->axpyWithRowAtPosNeg(alpha, y_pos, y_link_pos, y_neg, y_link_neg, row.getNode(), row.getIndex(), row.getLinking());
}

double StochRowStorage::multRowTimesVec( const INDEX& row, const StochVector& vec ) const
{
   assert(row.isRow());
   const double res = row_storage->localRowTimesVec(vec, row.getNode(), row.getIndex(), row.getLinking());
   return res;
}

double StochRowStorage::multLinkingRowTimesVecWithoutBl0( int row, const StochVector& vec) const
{
   const double res_full = row_storage->localRowTimesVec(vec, -1, row, true);
   const double res_bl0 = row_storage->Blmat->localRowTimesVec(dynamic_cast<const SimpleVector&>(*vec.vec), row);

   return res_full - res_bl0;
}

double StochRowStorage::getRowCoefficientAtColumn( const INDEX& row, const INDEX& col) const
{
   BlockType block_col = B_MAT;
   if( row.getLinking() )
      block_col = BL_MAT;
   if( col.getNode() == -1 && row.getNode() != -1)
      block_col = A_MAT;

   const SparseStorageDynamic& mat = *getSparseGenMatrixFromStochMat(*row_storage, block_col == BL_MAT ? col.getNode() : row.getNode(), block_col )->getStorageDynamic();

   const int row_start = mat.getRowPtr(row.getIndex()).start;
   const int row_end = mat.getRowPtr(row.getIndex()).end;

   for(int i = row_start; i < row_end; ++i)
   {
      if( mat.getJcolM(i) == col.getIndex() )
         return mat.getMat(i);
   }
   return 0.0;
}
