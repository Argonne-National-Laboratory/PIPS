/*
 * StochColumnStorage.h
 *
 *  Created on: 12.12.2019
 *      Author: Nils-Christian Kempke
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHCOLUMNSTORAGE_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHCOLUMNSTORAGE_H_

#include "StochGenMatrix.h"
#include "StochVector.h"
#include "SystemType.h"

class StochColumnStorage
{
public:
   StochColumnStorage( const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part);
   ~StochColumnStorage();

   int storeCol( int node, int col, const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part);

   double multColTimesVec( int node, int col, const StochVector& vec_eq, const StochVector& vec_ineq ) const;
   double multColTimesVecWithoutRootNode( int node, int col, const StochVector& vec_eq, const StochVector& vec_ineq ) const;
   double getColCoefficientAtRow( SystemType system_type, int node, int col, int row) const;

   // todo: delete Column from storage
private:
   // todo : assert that transposed is initalized
   /// all columns get stored as rows
   /// for a linking column the additional block B0 is stored in the additional associated SparseGenMatrix
   SparseGenMatrixHandle B0_eq;
   StochGenMatrixHandle stored_cols_eq;

   SparseGenMatrixHandle B0_ineq;
   StochGenMatrixHandle stored_cols_ineq;

   unsigned int nChildren;

   int storeLinkingCol(int col, const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part);
   int storeLocalCol(int node, int col, const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part);

   double multiplyLocalColTimesVec(int node, int col, const StochVector& vec_eq, const StochVector& vec_ineq) const;
   double multiplyLinkingColTimesVec(int col, const StochVector& vec_eq, const StochVector& vec_ineq) const;
   double multiplyLinkingColTimesVecWithoutRootNode(int col, const StochVector& vec_eq, const StochVector& vec_ineq) const;

   void createStorageMatrix(SparseGenMatrix* b0_block_storage, StochGenMatrix* col_storage, const StochGenMatrix& sys_matrix);
};




#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHCOLUMNSTORAGE_H_ */
