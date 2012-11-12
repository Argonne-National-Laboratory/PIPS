#ifndef OOQPRECOURSEINTERFACE_HPP
#define OOQPRECOURSEINTERFACE_HPP

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
class OOQPRecourseInterface {
public:
  OOQPRecourseInterface(stochasticInput &, int scenarioNumber, const std::vector<double> &firstStageSolution);

  void go();
  double getObjective() const;
  
  std::vector<double> getPrimalColSolution() const;
  std::vector<double> getDualRowSolution() const;
  
protected:
  scoped_ptr<FORMULATION> qp;
  scoped_ptr<QpGenData> prob;
  scoped_ptr<QpGenVars> vars;
  scoped_ptr<Residuals> resid;
  scoped_ptr<SOLVER> s;
  vector<int> primalOffsets;
  vector<int> rowToEqIneqIdx;
  vector<int> rowOffsets;
};


namespace{

  void scanRows(CoinPackedMatrix& Wrow, 
		std::vector<double>& lb, std::vector<double>& ub, 		
		int &my, int &mz, int &nnzA, int &nnzC, 
		vector<int> &rowToEqIneqIdx) 
  {
    
    my = mz = nnzA = nnzC = 0;
    
    int nrow = lb.size(), rowi=0;
    for (int i = 0; i < nrow; i++) {
      if (lb[i] == ub[i]) {
	my++;
	nnzA += Wrow.getVectorSize(i);
	rowToEqIneqIdx[rowi++]=my;
      } else {
	mz++;
	nnzC += Wrow.getVectorSize(i);
	rowToEqIneqIdx[rowi++]=-mz-1;
      }
    }
  }
  
  
  void separateRows(CoinPackedMatrix& Wrow, 
		    std::vector<double>& lb, std::vector<double>& ub,
		    SparseGenMatrix &A, SparseGenMatrix &C) 
  {
    int *rowptrA = A.krowM(), *rowptrC = C.krowM();
    int *colidxA = A.jcolM(), *colidxC = C.jcolM();
    double *eltA = A.M(), *eltC = C.M();
    
    int nnzA = 0, nnzC = 0;
    int nrowA = 0, nrowC = 0;
    
    int *rowptr, *colidx, *nnz, *nr;
    double *elt;
    
    int nrow = lb.size();
    for (int i = 0; i < nrow; i++) {
      if (lb[i] == ub[i]) {
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
	colidx[*nnz] = Wrow.getIndices()[k];
	elt[(*nnz)++] = Wrow.getElements()[k];
      }
    }
    rowptrA[nrowA] = nnzA;
    rowptrC[nrowC] = nnzC;
    
  }
  
  //we do not deal with cross terms here
  void formQ(stochasticInput &input, int scenNumber,
	     const CoinPackedMatrix& Qrow,
	     SparseSymMatrix &Q) 
  {
    int *rowQ = Q.krowM(), *colidx = Q.jcolM();
    double *elt = Q.M();
    
    int nnzQ = 0;
    int nrowQ = 0;
    //int offset = 0;
    
    int nrow = input.nSecondStageVars(scenNumber);
    for (int i = 0; i < nrow; i++) {
      rowQ[nrowQ++] = nnzQ;
      for (CoinBigIndex k = Qrow.getVectorFirst(i); k < Qrow.getVectorLast(i); k++) {
	//colidx[nnzQ] = Qrow.getIndices()[k]+offset;
	colidx[nnzQ] = Qrow.getIndices()[k];
	elt[nnzQ++] = Qrow.getElements()[k];
      }
    }
    rowQ[nrowQ] = nnzQ;
    //offset += nrow;
  }
}

template<typename SOLVER, typename FORMULATION>
OOQPRecourseInterface<SOLVER,FORMULATION>::OOQPRecourseInterface(stochasticInput &input, int scenNumber, const std::vector<double> &firstStageSolution)
{

  //Qrow.reverseOrderedCopyOf(input.getSecondStageHessian(s));
  
  int nvar1 = input.nFirstStageVars();
  int nvar2 = input.nSecondStageVars(scenNumber);
  int ncons2 = input.nSecondStageCons(scenNumber);

  CoinPackedMatrix Wrow, Trow;
  Wrow.reverseOrderedCopyOf(input.getSecondStageConstraints(scenNumber));
  Trow.reverseOrderedCopyOf(input.getLinkingConstraints(scenNumber));

  vector<double> Tx(ncons2);
  Trow.times(&firstStageSolution[0],&Tx[0]);

  vector<double> rowlb = input.getSecondStageRowLB(scenNumber),
    rowub = input.getSecondStageRowUB(scenNumber);
  
  for (int k = 0; k < ncons2; k++) {
    if (rowub[k] < 1e20) {
      rowub[k] -= Tx[k];
    }
    if (rowlb[k] >-1e20) {
      rowlb[k] -= Tx[k];
    }
  }

  rowToEqIneqIdx.resize(ncons2);

  int my, mz, nnzA, nnzC;
  scanRows(Wrow, rowlb, rowub, my, mz, nnzA, nnzC, rowToEqIneqIdx);


  int nnzQ = input.getSecondStageHessian(scenNumber).getNumElements();
  CoinPackedMatrix Qrow;
  Qrow.reverseOrderedCopyOf(input.getSecondStageHessian(scenNumber));


  qp.reset(new FORMULATION(nvar2,my,mz,nnzQ,nnzA,nnzC));
  prob.reset(dynamic_cast<QpGenData*>(qp->makeData()));
  
  separateRows(Wrow, rowlb, rowub, 
	       dynamic_cast<SparseGenMatrix&>(*prob->A),
	       dynamic_cast<SparseGenMatrix&>(*prob->C));

  formQ(input, scenNumber, Qrow,
	dynamic_cast<SparseSymMatrix&>(*prob->Q));

  
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
  //std::vector<double> const &l = input.getSecondStageRowLB(scenNumber);
  //std::vector<double> const &u = input.getSecondStageRowUB(scenNumber);
  int nrow = input.nSecondStageCons(scenNumber);
  for (int k = 0; k < nrow; k++) {
    if (rowlb[k] == rowub[k]) {
      bA[eq_idx++] = rowlb[k];
    } else {
      if (rowlb[k] > -1e20) {
	bl[ineq_idx] = rowlb[k];
	iclow[ineq_idx] = 1.;
      } else {
	bl[ineq_idx] = 0.;
	iclow[ineq_idx] = 0.;
      }
      if (rowub[k] < 1e20) {
	bu[ineq_idx] = rowub[k];
	icupp[ineq_idx] = 1.;
      } else {
	bu[ineq_idx] = 0.;
	icupp[ineq_idx] = 0.;
      }
      ineq_idx++;
    }
  }
  

  double RESCALE=1.0;
  int idx = 0;
  std::vector<double> const &l = input.getSecondStageColLB(scenNumber);
  std::vector<double> const &u = input.getSecondStageColUB(scenNumber);
  std::vector<double> const &c = input.getSecondStageObj(scenNumber);
  int ncol = input.nSecondStageVars(scenNumber);
  for (int k = 0; k < ncol; k++) {
    //g[idx] = RESCALE*c[k];
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

  vars.reset(dynamic_cast<QpGenVars*>(qp->makeVariables(prob.get())));
  resid.reset(qp->makeResiduals(prob.get()));
  s.reset(new SOLVER(qp.get(),prob.get()));

}


template<typename S, typename F>
void OOQPRecourseInterface<S,F>::go() 
{
  //s->monitorSelf();
  int result = s->solve(prob.get(),vars.get(),resid.get());
  
  // if ( 0 == result ) {
  //   double objective = prob->objectiveValue(vars.get());
    
  //   cout << " " << prob->nx << " variables, " 
  // 	 << prob->my  << " equality constraints, " 
  // 	 << prob->mz  << " inequality constraints.\n";
    
  //   cout << " Iterates: " << s->iter
  // 	 <<",    Optimal Solution:  " << objective << endl;
  // }
}

template<typename S, typename F>
double OOQPRecourseInterface<S,F>::getObjective() const {
	return prob->objectiveValue(vars.get());
}

template<typename S, typename F>
std::vector<double> OOQPRecourseInterface<S,F>::getPrimalColSolution() const
{
  double const *sol = &dynamic_cast<SimpleVector const&>(*vars->x)[0];
  return std::vector<double>(sol, vars->x->length());
}

template<typename S, typename F>
std::vector<double> OOQPRecourseInterface<S,F>::getDualRowSolution() const
{
  std::vector<double> out; out.reserve(vars->y->length()+vars->z->length());
  double const *eqsol = &dynamic_cast<SimpleVector const&>(*vars->y)[0];
  double const *ineqsol = &dynamic_cast<SimpleVector const&>(*vars->z)[0];

  for(int idx=0; idx<vars->y->length()+vars->z->length(); idx++) {
    int r=rowToEqIneqIdx[idx];
    if(r>=0) out.push_back(eqsol[r]);
    else     out.push_back(ineqsol[-r-1]);
  }
}

#endif
