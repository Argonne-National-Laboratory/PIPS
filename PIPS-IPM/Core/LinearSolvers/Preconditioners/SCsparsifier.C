/*
 * SCsparsifier.C
 *
 *  Created on: 21.06.2019
 *      Author: bzfrehfe
 */


#include "SCsparsifier.h"
#include "pipsdef.h"
#include "pipsport.h"
#include "StochOptions.h"
#include "mpi.h"
#include <cassert>

extern double g_iterNumber;
extern int gOuterBiCGFails;
extern int gOuterBiCGIter;
extern int gInnerBiCGIter;
extern int gInnerBiCGFails;


SCsparsifier::SCsparsifier(MPI_Comm mpiComm_)
{
   mpiComm = mpiComm_;
   diagDomBoundsPosition = 0;
}


SCsparsifier::~SCsparsifier()
{

}


double
SCsparsifier::getDiagDomBound() const
{
	assert(diagDomBoundsPosition < sizeof(diagDomBounds) / sizeof(diagDomBounds[0]));
	return diagDomBounds[diagDomBoundsPosition];
}


void
SCsparsifier::decreaseDiagDomBound(bool& success)
{
	const size_t positionUpperBound = sizeof(diagDomBounds) / sizeof(diagDomBounds[0]);

	assert(positionUpperBound > 0);

	if( diagDomBoundsPosition < positionUpperBound - 1 )
	{
	    int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

		diagDomBoundsPosition++;
		success = true;

		if( myRank == 0 )
			printf("\n SCsparsifier switched to sparsify factor %f (less aggressive) \n", diagDomBounds[diagDomBoundsPosition]);
	}
	else
	{
		success = false;
	}
}


void
SCsparsifier::increaseDiagDomBound(bool& success)
{
	if( diagDomBoundsPosition > 0 )
	{
	    int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

		diagDomBoundsPosition--;
		success = true;

		if( myRank == 0 )
			printf("\n SCsparsifier switched to sparsify factor %f (more aggressive)  \n", diagDomBounds[diagDomBoundsPosition]);
	}
	else
	{
		success = false;
	}
}


void
SCsparsifier::unmarkDominatedSCdistLocals(const sData& prob,
      SparseSymMatrix& sc) const
{
   const std::vector<bool>& rowIsMyLocal = prob.getSCrowMarkerMyLocal();

   const int* const krowM = sc.krowM();
   int* const jcolM = sc.jcolM();
   const double* const M = sc.M();
   const int sizeSC = sc.size();

   std::vector<double> diag = getDomDiagDist(prob, sc);

   /* loop over all rows */
   for( int r = 0; r < sizeSC; r++ )
   {
      const int rowEnd = krowM[r + 1];

      if( rowIsMyLocal[r] )
      {
         for( int j = krowM[r]; j < rowEnd; j++ )
         {
            const int col = jcolM[j];
            const double absM = fabs(M[j]);

            assert(col >= 0);

            if( (absM < epsilonZero && col != r) || (absM < diag[r] && absM < diag[col]) )
            {
               assert(col != r);
               jcolM[j] = -jcolM[j] - 1;
            }

         }
         continue;
      }

      for( int j = krowM[r]; j < rowEnd; j++ )
      {
         const int col = jcolM[j];
         assert(col >= 0);

         if( rowIsMyLocal[col] )
         {
            const double absM = fabs(M[j]);

            if( (absM < epsilonZero && col != r) || (absM < diag[r] && absM < diag[col]) )
            {
               assert(col != r);
               jcolM[j] = -jcolM[j] - 1;
            }
         }
      }
   } /* loop over all rows */
}


void
SCsparsifier::resetSCdistEntries(SparseSymMatrix& sc) const
{
   const int nnz = sc.numberOfNonZeros();
   int* const jcolM = sc.jcolM();

   for( int i = 0; i < nnz; i++ )
      if( jcolM[i] < 0 )
         jcolM[i] = -jcolM[i] - 1;

   assert(sc.getStorageRef().isValid(false));
}


void
SCsparsifier::getSparsifiedSC_fortran(const sData& prob,
      SparseSymMatrix& sc)
{
   assert(!sc.getStorageRef().fortranIndexed());

   int* const krowM = sc.krowM();
   int* const jcolM = sc.jcolM();
   double* const M = sc.M();
   const int sizeSC = sc.size();
   const std::vector<bool>& rowIsLocal = prob.getSCrowMarkerLocal();

   std::vector<double> diag(sizeSC);

// assert
   const double t = diagDomBounds[diagDomBoundsPosition];

   assert(t > 0.0 && t < 1.0);

   for( int r = 0; r < sizeSC; r++ )
   {
      const int j = krowM[r];
      assert(jcolM[j] == r);

      diag[r] = fabs(M[j]) * t;
   }

   int nnznew = 0;
   int rowStart = 0;
   krowM[0] = 1;

   for( int r = 0; r < sizeSC; r++ )
   {
      const int rowEnd = krowM[r + 1];

      if( rowIsLocal[r] )
      {
         for( int j = rowStart; j < rowEnd; ++j )
         {
            assert(fabs(M[j]) >= diag[r] || fabs(M[j]) >= diag[jcolM[j]]);
            assert(nnznew <= j);

            jcolM[nnznew] = jcolM[j] + 1;
            M[nnznew++] = M[j];
         }
      }
      else
      {
         for( int j = rowStart; j < rowEnd; ++j )
         {
            const double absM = fabs(M[j]);
            if( absM != 0.0 || jcolM[j] == r )
            {
               if( absM >= diag[r] || absM >= diag[jcolM[j]] )
               {
                  assert(nnznew <= j);

                  jcolM[nnznew] = jcolM[j] + 1;
                  M[nnznew++] = M[j];
               }
               else
               {
                  assert(jcolM[j] != r);
               }
            }
         }
      }

      rowStart = krowM[r + 1];
      krowM[r + 1] = nnznew + 1;
   }

   sc.getStorageRef().len = nnznew;

   sc.getStorageRef().set2FortranIndexed();
}


// todo should be done properly and needs to be tested....
void SCsparsifier::updateDiagDomBound()
{
   bool wasIncreased = false;
   const int nIter = std::max(gOuterBiCGIter, gInnerBiCGIter);
   const int nFails = std::max(gOuterBiCGFails, gInnerBiCGFails);

   assert(nIter >= 0);
   assert(nFails >= 0);

   if( nIter >= 5 && static_cast<int>(g_iterNumber) > 0 )
   {
      if( diagDomBoundsPosition == 0 )
      {
    	 decreaseDiagDomBound(wasIncreased);
    	 assert(wasIncreased);
      }
   }

   if( nFails >= 3 )
   {
      if( diagDomBoundsPosition == 1 )
      {
    	  decreaseDiagDomBound(wasIncreased);
    	  assert(wasIncreased);
      }
   }

   if( nFails >= 20 )
   {
	   if( diagDomBoundsPosition == 2 )
	   {
		   decreaseDiagDomBound(wasIncreased);
	  	   assert(wasIncreased);
	   }
   }

   if( nFails >= 60 )
   {
	   if( diagDomBoundsPosition == 3 )
	   {
		   decreaseDiagDomBound(wasIncreased);
		   assert(wasIncreased);
	   }
   }
}


std::vector<double> SCsparsifier::getDomDiagDist(const sData& prob, SparseSymMatrix& sc) const
{
   int* const krowM = sc.krowM();
   double* const M = sc.M();
   const int sizeSC = sc.size();
   const std::vector<bool>& rowIsLocal = prob.getSCrowMarkerLocal();
   const std::vector<bool>& rowIsMyLocal = prob.getSCrowMarkerMyLocal();

   std::vector<double> diag(sizeSC, 0.0);

   for( int r = 0; r < sizeSC; r++ )
   {
      if( krowM[r] == krowM[r + 1] )
         continue;

      assert(sc.jcolM()[krowM[r]] == r);

      if( rowIsLocal[r] && !rowIsMyLocal[r] )
         continue;

      diag[r] = M[krowM[r]];
   }

   MPI_Allreduce(MPI_IN_PLACE, &diag[0], sizeSC, MPI_DOUBLE, MPI_SUM, mpiComm);

   const double diagDomBound = getDiagDomBound();

   for( size_t i = 0; i < diag.size(); ++i )
      diag[i] = fabs(diag[i]) * diagDomBound;

   return diag;
}

