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
#include "SystemType.h"

class StochRowStorage
{
public:
   StochRowStorage(const StochGenMatrix& system_matrix);
   ~StochRowStorage();

   int storeRow( const INDEX& row, const StochGenMatrix& matrix_row);

   /** y = beta * y + alpha * stored */
   void axpyAtRow(double beta, StochVector& y, double alpha, const INDEX& row ) const;

   double multRowTimesVec(const INDEX& row, const StochVector &vec) const;
   double multLinkingRowTimesVecWithoutBl0(int row, const StochVector &vec) const;

   double getRowCoefficientAtColumn( const INDEX& row, const INDEX& col ) const;

   // todo : deleteRowFromStorage
private:
   StochGenMatrixHandle row_storage;

};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHROWSTORAGE_H_ */
