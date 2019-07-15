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
      constexpr static double diagDomBoundAggressive = 0.001;
      constexpr static double diagDomBoundNormal = 0.0002;
      constexpr static double diagDomBoundConservative = 0.000025;
      constexpr static double diagDomBoundDefault = diagDomBoundAggressive;


      constexpr static double epsilonZero = 1e-15;

      SCsparsifier();

      // for any diagDomBound with value <= 0.0 the default value will be taken
      SCsparsifier(double diagDomBound, MPI_Comm mpiComm);
      ~SCsparsifier();

      double getDiagDomBound() const { return diagDomBound; };

      // set columns col of dominated local Schur complement (distributed) entries to -col
      void unmarkDominatedSCdistLocals(const sData& prob, SparseSymMatrix& sc) const;

      // resets unmarkDominatedSCdistEntries actions
      void resetSCdistEntries(SparseSymMatrix& sc) const;

      // deletes dominated Schur complement entries and converts matrix to Fortran format
      void getSparsifiedSC_fortran(const sData& prob, SparseSymMatrix& sc);

   private:

      double diagDomBound;
      MPI_Comm mpiComm;

      void updateDiagDomBound();
      std::vector<double> getDomDiagDist(const sData& prob, SparseSymMatrix& sc) const;

};


#endif /* PIPS_IPM_CORE_LINEARSOLVERS_PRECONDITIONERS_SCSPARSIFIER_H_ */
