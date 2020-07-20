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

   int storeCol( const INDEX& col, const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part);

   /** y = beta * vec + alpha * stored - either of the vectors can be nullptr and will then not be considered */
   /** there is the option of specifying a special vector to store the linking row results in - if no vector is given the corresponding link vector will be ignored */
   void axpyAtCol( double beta, StochVector* eq_vec, StochVector* ineq_vec, SimpleVector* eq_link, SimpleVector* ineq_link, double alpha, const INDEX& col ) const;

   double multColTimesVec( const INDEX& col, const StochVector& vec_eq, const StochVector& vec_ineq ) const;
   double multColTimesVecWithoutRootNode( const INDEX& col, const StochVector& vec_eq, const StochVector& vec_ineq ) const;

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
   int storeLocalCol( const INDEX& col, const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part);

   double multiplyLocalColTimesVec( const INDEX& col, const StochVector& vec_eq, const StochVector& vec_ineq) const;
   double multiplyLinkingColTimesVec(int col, const StochVector& vec_eq, const StochVector& vec_ineq) const;
   double multiplyLinkingColTimesVecWithoutRootNode(int col, const StochVector& vec_eq, const StochVector& vec_ineq) const;

   void createStorageMatrix(SystemType system_type, const StochGenMatrix& sys_matrix);

   /* calculate vec + alpha * stored for system */
   void axpyAtCol(StochVector& vec, SimpleVector* vec_link, double alpha, const INDEX& col, SystemType system_type) const;
};




#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHCOLUMNSTORAGE_H_ */
