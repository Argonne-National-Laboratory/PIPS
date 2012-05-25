#ifndef OOQPINTERFACE_HPP
#define OOQPINTERFACE_HPP

#include "stochasticInput.hpp"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "SimpleVector.h"
#include "Status.h"
#include "QpGenData.h"
#include "OoqpVersion.h"
#include <boost/scoped_ptr.hpp>
#include "SparseGenMatrix.h"
#include "SparseSymMatrix.h"

using boost::scoped_ptr;


template<typename SOLVER, typename FORMULATION>
class OOQPInterface {
public:
	OOQPInterface(stochasticInput &);

	void go();


protected:
	scoped_ptr<FORMULATION> qp;
	scoped_ptr<QpGenData> prob;
	scoped_ptr<QpGenVars> vars;
	scoped_ptr<Residuals> resid;
	scoped_ptr<SOLVER> s;

};


namespace{

void scanRows(stochasticInput &input, int &my, int &mz, int &nnzA, int &nnzC) {
	
	my = mz = nnzA = nnzC = 0;
	{
		CoinPackedMatrix Arow;
		Arow.reverseOrderedCopyOf(input.getFirstStageConstraints());
		int nrow1 = input.nFirstStageCons();
		std::vector<double> const &l = input.getFirstStageRowLB(), &u = input.getFirstStageRowUB();
		for (int i = 0; i < nrow1; i++) {
			if (l[i] == u[i]) {
				my++;
				nnzA += Arow.getVectorSize(i);
			} else {
				mz++;
				nnzC += Arow.getVectorSize(i);
			}
		}
	}

	int nscen = input.nScenarios();
	for (int i = 0; i < nscen; i++) {
		CoinPackedMatrix Wrow, Trow;
		Wrow.reverseOrderedCopyOf(input.getSecondStageConstraints(i));
		Trow.reverseOrderedCopyOf(input.getLinkingConstraints(i));
		std::vector<double> const &l = input.getSecondStageRowLB(i), &u = input.getSecondStageRowUB(i);
		int nrow2 = input.nSecondStageCons(i);
		for (int k = 0; k < nrow2; k++) {
			if (l[k] == u[k]) {
				my++;
				nnzA += Wrow.getVectorSize(k) + Trow.getVectorSize(k);
			} else {
				mz++;
				nnzC += Wrow.getVectorSize(k) + Trow.getVectorSize(k);
			}
		}
	}
	
}


void separateRows(stochasticInput &input, SparseGenMatrix &A, SparseGenMatrix &C) {

	int nscen = input.nScenarios();
	int *rowptrA = A.krowM(), *rowptrC = C.krowM();
	int *colidxA = A.jcolM(), *colidxC = C.jcolM();
	double *eltA = A.M(), *eltC = C.M();

	int offset1 = 0;
	for (int i = 0; i < nscen; i++) offset1 += input.nSecondStageVars(i);

	int nnzA = 0, nnzC = 0;
	int nrowA = 0, nrowC = 0;
	
	int *rowptr, *colidx, *nnz, *nr;
	double *elt;

	// [ W T ]
	int offset2 = 0;
	for (int s = 0; s < nscen; s++) {
		CoinPackedMatrix Wrow, Trow;
		Wrow.reverseOrderedCopyOf(input.getSecondStageConstraints(s));
		Trow.reverseOrderedCopyOf(input.getLinkingConstraints(s));
		
		std::vector<double> const &l = input.getSecondStageRowLB(s), &u = input.getSecondStageRowUB(s);
		int nrow = input.nSecondStageCons(s);
		for (int i = 0; i < nrow; i++) {
			if (l[i] == u[i]) {
				rowptr = rowptrA;
				colidx = colidxA;
				elt = eltA;
				nr = &nrowA;
				nnz = &nnzA;
			} else {
				rowptr = rowptrC;
				colidx = colidxC;
				elt = eltC;
				nr = &nrowC;
				nnz = &nnzC;
			}
			rowptr[(*nr)++] = *nnz;
			for (CoinBigIndex k = Wrow.getVectorFirst(i); k < Wrow.getVectorLast(i); k++) {
				colidx[*nnz] = Wrow.getIndices()[k]+offset2;
				elt[(*nnz)++] = Wrow.getElements()[k];
			}
			for (CoinBigIndex k = Trow.getVectorFirst(i); k < Trow.getVectorLast(i); k++) {
				colidx[*nnz] = Trow.getIndices()[k]+offset1;
				elt[(*nnz)++] = Trow.getElements()[k];
			}
		}
		offset2 += input.nSecondStageVars(s);
	}

	CoinPackedMatrix Arow;
	Arow.reverseOrderedCopyOf(input.getFirstStageConstraints());
	int nrow = input.nFirstStageCons();
	std::vector<double> const &l = input.getFirstStageRowLB(), &u = input.getFirstStageRowUB();
	for (int i = 0; i < nrow; i++) {
		if (l[i] == u[i]) {
			rowptr = rowptrA;
			colidx = colidxA;
			elt = eltA;
			nr = &nrowA;
			nnz = &nnzA;
		} else {
			rowptr = rowptrC;
			colidx = colidxC;
			elt = eltC;
			nr = &nrowC;
			nnz = &nnzC;
		}
		rowptr[(*nr)++] = *nnz;
		for (CoinBigIndex k = Arow.getVectorFirst(i); k < Arow.getVectorLast(i); k++) {
			colidx[*nnz] = Arow.getIndices()[k]+offset1;
			elt[(*nnz)++] = Arow.getElements()[k];
		}
	}
	rowptrA[nrowA] = nnzA;
	rowptrC[nrowC] = nnzC;

}


void formQ(stochasticInput &input, SparseSymMatrix &Q) {
	/*
	[ Q_1
	       Q_2
	            ...
	                Q_N
	 Q'_1 Q'_2     Q'_N Q ]
	 Cross terms given by input are transposed forms of Q'_i
	 */

	int nscen = input.nScenarios();
	int *rowQ = Q.krowM(), *colidx = Q.jcolM();
	double *elt = Q.M();
	
	int nnzQ = 0;
	int nrowQ = 0;
	int offset = 0;
	for (int s = 0; s < nscen; s++) {
		CoinPackedMatrix Qrow;
		Qrow.reverseOrderedCopyOf(input.getSecondStageHessian(s));
		int nrow = input.nSecondStageVars(s);
		for (int i = 0; i < nrow; i++) {
			rowQ[nrowQ++] = nnzQ;
			for (CoinBigIndex k = Qrow.getVectorFirst(i); k < Qrow.getVectorLast(k); k++) {
				colidx[nnzQ] = Qrow.getIndices()[k]+offset;
				elt[nnzQ++] = Qrow.getElements()[k];
			}
		}
		offset += nrow;
	}
	int nvar1 = input.nFirstStageVars();

	CoinPackedMatrix Qrow;
	Qrow.reverseOrderedCopyOf(input.getFirstStageHessian());
	for (int i = 0; i < nvar1; i++) {
		rowQ[nrowQ++] = nnzQ;
		offset = 0;
		for (int s = 0; s < nscen; s++) {
			CoinPackedMatrix const& Qcross = input.getSecondStageCrossHessian(s); // this will make a new copy for each row, oops
			for (CoinBigIndex k = Qcross.getVectorFirst(i); k < Qcross.getVectorLast(i); k++) {
				colidx[nnzQ] = Qcross.getIndices()[k]+offset;
				elt[nnzQ++] = Qcross.getElements()[k];
			}
			offset += input.nSecondStageVars(s);
		}
		for (CoinBigIndex k = Qrow.getVectorFirst(i); k < Qrow.getVectorLast(i); k++) {
			colidx[nnzQ] = Qrow.getIndices()[k] + offset;
			elt[nnzQ++] = Qrow.getElements()[k];
		}
	}
	rowQ[nrowQ] = nnzQ;

}

}


template<typename SOLVER, typename FORMULATION>
OOQPInterface<SOLVER,FORMULATION>::OOQPInterface(stochasticInput &input) {

	// we put the first-stage variables at the end to help with the sparse reordering
	int nvar1 = input.nFirstStageVars();
	int nscen = input.nScenarios();
	int totalVars = nvar1;
	for (int i = 0; i < nscen; i++) {
		totalVars += input.nSecondStageVars(i);
	}

	int my, mz, nnzA, nnzC;

	scanRows(input, my, mz, nnzA, nnzC);

	int nnzQ = input.getFirstStageHessian().getNumElements();
	for (int i = 0; i < nscen; i++) {
		nnzQ += input.getSecondStageHessian(i).getNumElements();
		nnzQ += input.getSecondStageCrossHessian(i).getNumElements();
	}

	qp.reset(new FORMULATION(totalVars,my,mz,nnzQ,nnzA,nnzC));
	prob.reset(dynamic_cast<QpGenData*>(qp->makeData()));

	separateRows(input,dynamic_cast<SparseGenMatrix&>(*prob->A),dynamic_cast<SparseGenMatrix&>(*prob->C));
	formQ(input,dynamic_cast<SparseSymMatrix&>(*prob->Q));

	// the beauty of OOP...
	double *bA = dynamic_cast<SimpleVector&>(*prob->bA).elements();
	double *bl = dynamic_cast<SimpleVector&>(*prob->bl).elements();
	double *bu = dynamic_cast<SimpleVector&>(*prob->bu).elements();
	double *iclow = dynamic_cast<SimpleVector&>(*prob->iclow).elements();
	double *icupp = dynamic_cast<SimpleVector&>(*prob->icupp).elements();
	double *g = dynamic_cast<SimpleVector&>(*prob->g).elements();
	double *blx = dynamic_cast<SimpleVector&>(*prob->blx).elements();
	double *bux = dynamic_cast<SimpleVector&>(*prob->bux).elements();
	double *ixupp = dynamic_cast<SimpleVector&>(*prob->ixupp).elements();
	double *ixlow = dynamic_cast<SimpleVector&>(*prob->ixlow).elements();

	int eq_idx = 0, ineq_idx = 0;
	for (int i = 0; i < nscen; i++) {
		std::vector<double> const &l = input.getSecondStageRowLB(i), &u = input.getSecondStageRowUB(i);
		int nrow = input.nSecondStageCons(i);
		for (int k = 0; k < nrow; k++) {
			if (l[k] == u[k]) {
				bA[eq_idx++] = l[k];
			} else {
				if (l[k] > -1e20) {
					bl[ineq_idx] = l[k];
					iclow[ineq_idx] = 1.;
				} else {
					bl[ineq_idx] = 0.;
					iclow[ineq_idx] = 0.;
				}
				if (u[k] < 1e20) {
					bu[ineq_idx] = u[k];
					icupp[ineq_idx] = 1.;
				} else {
					bu[ineq_idx] = 0.;
					icupp[ineq_idx] = 0.;
				}
				ineq_idx++;
			}
		}
	}
	{
		std::vector<double> const &l = input.getFirstStageRowLB(), &u = input.getFirstStageRowUB();
		int nrow = input.nFirstStageCons();
		for (int k = 0; k < nrow; k++) {
			if (l[k] == u[k]) {
				bA[eq_idx++] = l[k];
			} else {
				if (l[k] > -1e20) {
					bl[ineq_idx] = l[k];
					iclow[ineq_idx] = 1.;
				} else {
					bl[ineq_idx] = 0.;
					iclow[ineq_idx] = 0.;
				}
				if (u[k] < 1e20) {
					bu[ineq_idx] = u[k];
					icupp[ineq_idx] = 1.;
				} else {
					bu[ineq_idx] = 0.;
					icupp[ineq_idx] = 0.;
				}
				ineq_idx++;
			}
		}
	}

	int idx = 0;
	for (int i = 0; i < nscen; i++) {
		std::vector<double> const &l = input.getSecondStageColLB(i), &u = input.getSecondStageColUB(i),
			&c = input.getSecondStageObj(i);
		int ncol = input.nSecondStageVars(i);
		for (int k = 0; k < ncol; k++) {
			g[idx] = c[k];
			if (l[k] > -1e20) {
				blx[idx] = l[k];
				ixlow[idx] = 1.;
			} else {
				blx[idx] = 0.;
				ixlow[idx] = 0.;
			}
			if (u[k] < 1e20) {
				bux[idx] = u[k];
				ixupp[idx] = 1.;
			} else {
				bux[idx] = 0.;
				ixupp[idx] = 0.;
			}
			idx++;
		}
	}
	{
		std::vector<double> const &l = input.getFirstStageColLB(), &u = input.getFirstStageColUB(),
			&c = input.getFirstStageObj();
		int ncol = input.nFirstStageVars();
		for (int k = 0; k < ncol; k++) {
			g[idx] = c[k];
			if (l[k] > -1e20) {
				blx[idx] = l[k];
				ixlow[idx] = 1.;
			} else {
				blx[idx] = 0.;
				ixlow[idx] = 0.;
			}
			if (u[k] < 1e20) {
				bux[idx] = u[k];
				ixupp[idx] = 1.;
			} else {
				bux[idx] = 0.;
				ixupp[idx] = 0.;
			}
			idx++;
		}


	}



	vars.reset(dynamic_cast<QpGenVars*>(qp->makeVariables(prob.get())));
	resid.reset(qp->makeResiduals(prob.get()));
	s.reset(new SOLVER(qp.get(),prob.get()));

}


template<typename SOLVER, typename FORMULATION>
void OOQPInterface<SOLVER,FORMULATION>::go() {

	s->monitorSelf();
	int result = s->solve(prob.get(),vars.get(),resid.get());

	if ( 0 == result ) {
	 double objective = prob->objectiveValue(vars.get());
      
      cout << " " << prob->nx << " variables, " 
	   << prob->my  << " equality constraints, " 
	   << prob->mz  << " inequality constraints.\n";
      
      cout << " Iterates: " << s->iter
	   <<",    Optimal Solution:  " << objective << endl;
      }
	


}

#endif
