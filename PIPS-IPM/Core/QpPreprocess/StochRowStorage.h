/*
 * StochRowStorage.h
 *
 *  Created on: 12.12.2019
 *      Author: Nils-Christian Kempke
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHROWSTORAGE_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHROWSTORAGE_H_

#include "StochGenMatrix.h"
#include "StochVector.h"

class StochRowStorage
{
public:
   StochRowStorage(const StochGenMatrix& system_matrix);
   ~StochRowStorage();

   int storeRow( int node, int row, bool linking_row, const StochGenMatrix& matrix_row);

   double multRowTimesVec( int node, int row, bool linking_row, const StochVector& vec ) const;
   double multLinkingRowTimesVecWithoutB0( int row, const StochVector& vec) const;

   double getRowCoefficientAtColumn( int node_row, int row, bool linking_row, int node_column, int colum ) const;

   // todo : deleteRowFromStorage
private:
   StochGenMatrixHandle row_storage;

};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHROWSTORAGE_H_ */
