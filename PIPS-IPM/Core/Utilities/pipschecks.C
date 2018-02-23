/*
 * pipschecks.C
 *
 *  Created on: 23.02.2018
 *      Author: bzfrehfe
 */

#include "pipschecks.h"
#include <cassert>

// are the columns of the given sub-matrix ordered?
bool subMatrixIsOrdered(const int* rowptr, const int* colidx,
      int rowstart, int rowend)
{
   assert(rowptr);
   assert(colidx);
   assert(rowstart >= 0 && rowstart <= rowend);

   for( int r = rowstart; r < rowend; r++ )
      for( int ci = rowptr[r] + 1; ci < rowptr[r + 1]; ci++ )
         if( colidx[ci - 1] >= colidx[ci] )
            return false;

   return true;
}

