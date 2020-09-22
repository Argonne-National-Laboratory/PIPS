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

#define SCSPARSIFIER_SAVE_STATS

class sData;

const static double diagDomBoundsLeaf[] = {   0.002,
										      0.001,
										      0.0003,
										      0.000025,
										      0.000005,
										      0.000001   };

const static double diagDomBounds[] = {   0.001,
										  0.0005,
										  0.0002,
										  0.000025,
										  0.000005,
										  0.000001   };

/** Schur complement sparsifier (based on comparison with diagonal entries) */
class SCsparsifier
{
   public:

      unsigned diagDomBoundsPosition;

      constexpr static double epsilonZero = 1e-15;

      SCsparsifier(MPI_Comm mpiComm = MPI_COMM_NULL);
      ~SCsparsifier();

      // returns sparsification bound
      double getDiagDomBound() const;

      double getDiagDomBoundLeaf() const;

      // increases sparsification (more aggressive)
      void increaseDiagDomBound(bool& success);

      // decreases sparsification (less aggressive)
      void decreaseDiagDomBound(bool& success);

      void updateStats();

      // sets CSR column marker col of dominated local Schur complement (distributed) entries to -col
      void unmarkDominatedSCdistLocals(const sData& prob, SparseSymMatrix& sc);

      // resets unmarkDominatedSCdistEntries actions
      void resetSCdistEntries(SparseSymMatrix& sc) const;

      // deletes dominated Schur complement entries and converts matrix to Fortran format
      void getSparsifiedSC_fortran(const sData& prob, SparseSymMatrix& sc);

      // updates the bound according the convergence history of BICGStab
      // todo this method should be removed and increaseDiagDomBound/decreaseDiagDomBound should be used instead
      // by the solver
      void updateDiagDomBound();


   private:
#ifdef SCSPARSIFIER_SAVE_STATS
      std::vector<double> allratios;
      double ratioAvg;
      int nEntriesLocal;
      int nDeletedLocal;
#endif
      double diagDomBound;
      double diagDomBoundLeaf;


      MPI_Comm mpiComm;

      std::vector<double> getDomDiagDist(const sData& prob, SparseSymMatrix& sc, bool isLeaf = false) const;

};


#endif /* PIPS_IPM_CORE_LINEARSOLVERS_PRECONDITIONERS_SCSPARSIFIER_H_ */
