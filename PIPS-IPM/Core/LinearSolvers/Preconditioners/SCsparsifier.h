/*
 * SCsparsifier.h
 *
 *  Created on: 21.06.2019
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_LINEARSOLVERS_PRECONDITIONERS_SCSPARSIFIER_H_
#define PIPS_IPM_CORE_LINEARSOLVERS_PRECONDITIONERS_SCSPARSIFIER_H_

#include "pipsport.h"
#include "sData.h"
#include <vector>

class sData;


/** Schur complement sparsifier (based comparison with diagonal entries) */
class SCsparsifier
{
   public:
      static constexpr double diagDomBoundDefault = 0.0001;
      static constexpr double epsilonZero = 1e-15;

      SCsparsifier();

      // for any diagDomBound with value <= 0.0 the default value will be taken
      SCsparsifier(double diagDomBound, MPI_Comm mpiComm);
      ~SCsparsifier();

      double getDiagDomBound() const { return diagDomBound; };

      // set columns col of dominated Schur complement (distributed) entries to -col
      void unmarkDominatedSCdistEntries(const sData& prob, SparseSymMatrix& sc) const;

      // resets unmarkDominatedSCdistEntries actions
      void resetSCdistEntries(SparseSymMatrix& sc) const;

      // deletes dominated Schur complement entries and converts matrix to Fortran format
      void getSparsifiedSC_fortran(const sData& prob, SparseSymMatrix& sc) const;

   private:

      double diagDomBound;
      MPI_Comm mpiComm;

      std::vector<double> getDomDiagDist(SparseSymMatrix& sc) const;

};


#endif /* PIPS_IPM_CORE_LINEARSOLVERS_PRECONDITIONERS_SCSPARSIFIER_H_ */
