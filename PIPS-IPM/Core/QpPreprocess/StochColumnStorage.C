/*
 * StochColumnStorage.C
 *
 *  Created on: 12.12.2019
 *      Author: Nils-Christian Kempke
 */

#include "StochColumnStorage.h"


StochColumnStorage::StochColumnStorage( const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part)
{
   // todo
}

StochColumnStorage::~StochColumnStorage()
{
}

int StochColumnStorage::storeCol( int node, int col, const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part)
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
