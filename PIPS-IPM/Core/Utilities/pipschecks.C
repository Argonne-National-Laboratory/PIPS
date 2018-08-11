/*
 * pipschecks.C
 *
 *  Created on: 23.02.2018
 *      Author: bzfrehfe
 */

#include "pipschecks.h"
#include <math.h>
#include <cassert>

bool permutationIsValid(const std::vector<unsigned int>& perm)
{
   size_t size = perm.size();
   std::vector<bool> permflag(size, false);

   for( size_t i = 0; i < size; ++i )
   {
      const unsigned int p = perm[i];
      if( p >= size )
         return false;
      permflag[p] = true;
   }

   for( size_t i = 0; i < size; ++i )
      if( !permflag[i] )
         return false;

   return true;
}


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

void computeFortranCSRMatResidualNorms(const int* rowptr, const int* colidx, const double* vals, /*const*/ SimpleVector& rhs,
      /*const*/ SimpleVector& x, double& res_norm2, double& res_nrmInf, double& sol_inf, double& mat_max)
{
   const int dim=rhs.length();
   double* tmp_resid = new double[dim];
   memcpy(tmp_resid, rhs.elements(), dim * sizeof(double));

   mat_max = 0.0, res_norm2 = 0.0, res_nrmInf = 0.0, sol_inf = 0.0;

   for( int i = 0; i < dim; i++ )
   {
      for( int p = rowptr[i]; p < rowptr[i + 1]; p++ )
      {
         const int j = colidx[p - 1] - 1;
         if( j + 1 <= dim )
         {
            //r[i] = r[i] + M(i,j)*x(j)
            tmp_resid[i] -= vals[p - 1] * x[j];

            if( fabs(vals[p - 1]) > mat_max )
               mat_max = fabs(vals[p - 1]);

            if( j != i )
            {
               //r[j] = r[j] + M(j,i)*x(i)
               tmp_resid[j] -= vals[p - 1] * x[i];
            }
         }
      }
   }

   for( int i = 0; i < dim; i++ )
   {
      res_norm2 += tmp_resid[i] * tmp_resid[i];
      if( res_nrmInf < fabs(tmp_resid[i]) )
         res_nrmInf = tmp_resid[i];
      if( fabs(x[i]) > sol_inf )
         sol_inf = fabs(x[i]);
   }
   res_norm2 = sqrt(res_norm2);

   delete[] tmp_resid;
}
