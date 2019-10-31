/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#ifndef SROOTLINSYSCOMM2
#define SROOTLINSYSCOMM2

#include "sLinsysRoot.h"

class sFactory;
class sData;

// DEBUG only
//#include "ScaDenSymMatrix.h"

/**
 * ROOT (= NON-leaf) linear system
 */
class sLinsysRootComm2 : public sLinsysRoot {
 protected:
  sLinsysRootComm2() {};
 public:
  sLinsysRootComm2(sFactory * factory_, sData * prob_);
  sLinsysRootComm2(sFactory* factory,
	      sData* prob_,
	      OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_,
	      OoqpVector* rhs_);

  virtual SymMatrix*   createKKT     (sData* prob) = 0;
  virtual DoubleLinearSolver*
                       createSolver  (sData* prob,
				      SymMatrix* kktmat) = 0;


  virtual void factor2(sData *prob, Variables *vars);
  /* Atoms methods of FACTOR2 for a non-leaf linear system */
  using sLinsysRoot::reduceKKT;
  virtual void reduceKKT();

  using sLinsysRoot::factorizeKKT;
  virtual void factorizeKKT();

  virtual void finalizeKKT(sData* prob, Variables* vars)=0;

  virtual void Dsolve ( sData *prob, OoqpVector& x);
  virtual void Lsolve ( sData *prob, OoqpVector& x );
  virtual void Ltsolve( sData *prob, OoqpVector& x );
  virtual void solveReduced( sData *prob, SimpleVector& b)=0;
 public:
  virtual ~sLinsysRootComm2();

  void submatrixReduce(DenseSymMatrix* A,
			  int row, int col, int drow, int dcol,
			  MPI_Comm comm);

};

#endif
