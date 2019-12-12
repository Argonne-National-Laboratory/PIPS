/*
 * StochRowStorage.C
 *
 *  Created on: 12.12.2019
 *      Author: Nils-Christian Kempke
 */

#include "StochRowStorage.h"


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


double StochRowStorage::multRowTimesVec( int node, int row, bool linking_row, const StochVector& vec ) const
{
   /// if we multiply a linking row but are not root - we do not add the bl0 matrix to the result
   return 0;
}

double StochRowStorage::getRowCoefficientAtColumn( int node, int row, bool linking_row, int colum ) const
{
   // todo
   return 0;
}
