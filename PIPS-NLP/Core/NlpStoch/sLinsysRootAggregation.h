/* PIPS-NLP                                                           	*
 * Author:  Nai-Yuan Chiang                                       	*
 * (C) 2015 Argonne National Laboratory. 				*/

#ifndef SROOTLINSYSAGG
#define SROOTLINSYSAGG

#include "sLinsysRoot.h"

class sFactory;
class sData;

// DEBUG only
//#include "ScaDenSymMatrix.h"

/** 
 * ROOT (= NON-leaf) linear system
 */
class sLinsysRootAggregation : public sLinsysRoot {
 protected:
  sLinsysRootAggregation() {};
 public:
  sLinsysRootAggregation(sFactory * factory_, sData * prob_);
  sLinsysRootAggregation(sFactory* factory,
	      sData* prob_,				    
	      OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_,
	      OoqpVector* rhs_, OoqpVector* additiveDiag);
  virtual ~sLinsysRootAggregation();

  virtual SymMatrix*   createKKT     (sData* prob);
  virtual DoubleLinearSolver* 
                       createSolver  (sData* prob, 
				      SymMatrix* kktmat);


  virtual int factor2(sData *prob, Variables *vars);
  /* Atoms methods of FACTOR2 for a non-leaf linear system */
  
  virtual void initializeKKT(sData* prob, Variables* vars){};
  virtual void reduceKKT();
  virtual int factorizeKKT(); 
  virtual void finalizeKKT(sData* prob, Variables* vars){};

  virtual void Dsolve ( sData *prob, OoqpVector& x);
  virtual void Lsolve ( sData *prob, OoqpVector& x );
  virtual void Ltsolve( sData *prob, OoqpVector& x );
  virtual void solveReduced( sData *prob, SimpleVector& b);


  void submatrixReduce(DenseSymMatrix* A, 
			  int row, int col, int drow, int dcol,
			  MPI_Comm comm);


  SparseSymMatrix* MatQ;
  SparseGenMatrix* MatJ;
  
  SimpleVector* redRhs;
  SimpleVector* temp_nrow;
  SimpleVector* temp_ncol;

  int redDim;
  int redNnz;
  
  int QNnz;
  int JNnz;

  int n_col;
  int n_row;


  // dummy dense mat
  DenseSymMatrix* CtDC;

  SimpleVector* temp_ixlow;
  SimpleVector* temp_ixupp;

  // vectors for the reduced system
  SimpleVector *rQ;
  SimpleVector *rA;
  SimpleVector *rv;
  SimpleVector *rw;
  SimpleVector *rgamma;
  SimpleVector *rphi;

  SimpleVector* temp_x;
  SimpleVector* temp_y;
  SimpleVector* temp_gamma;
  SimpleVector* temp_phi;
  SimpleVector* temp_v;
  SimpleVector* temp_w;

  SimpleVector* temp_rhs_x;
  SimpleVector* temp_rhs_y;

  SimpleVector* temp_xl;
  SimpleVector* temp_xu; 

  bool setIX;

  protected: 
	void calcPreCondKKTResids(Data *prob_in);	
	void computeReducedRhs();

  private:
	void _joinRedRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in, OoqpVector& rhs2_in);
	void _separateVars( OoqpVector& vars_in, OoqpVector& x_in, OoqpVector& y_in );
	void _set1stStVar(SimpleVector& redSol,SimpleVector& sol0);
};

#endif

