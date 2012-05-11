/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#ifndef SROOTLINSYS
#define SROOTLINSYS

#include "sLinsys.h"
class StochTree;
class sFactory;
class QpGenStochData;

// DEBUG only
//#include "ScaDenSymMatrix.h"

/** 
 * ROOT (= NON-leaf) linear system
 */
class sLinsysRoot : public sLinsys {
 protected:
  sLinsysRoot() {};

  virtual void         createChildren(QpGenStochData* prob);
  virtual void         deleteChildren();

  virtual SymMatrix*   createKKT     (QpGenStochData* prob) = 0;
  virtual DoubleLinearSolver* 
                       createSolver  (QpGenStochData* prob, 
				      SymMatrix* kktmat) = 0;
 public:
  std::vector<sLinsys*> children;
  int iAmDistrib;
 public:

  sLinsysRoot(sFactory * factory_, QpGenStochData * prob_);
  sLinsysRoot(sFactory* factory,
	      QpGenStochData* prob_,				    
	      OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_,
	      OoqpVector* rhs_);

  virtual void factor2(QpGenStochData *prob, Variables *vars);
  /* Atoms methods of FACTOR2 for a non-leaf linear system */
  virtual void initializeKKT(QpGenStochData* prob, Variables* vars);
  virtual void reduceKKT();
  virtual void factorizeKKT(); 
  virtual void finalizeKKT(QpGenStochData* prob, Variables* vars)=0;


  virtual void Lsolve ( QpGenStochData *prob, OoqpVector& x );
  virtual void Dsolve ( QpGenStochData *prob, OoqpVector& x );
  virtual void Ltsolve( QpGenStochData *prob, OoqpVector& x );

  virtual void Ltsolve2( QpGenStochData *prob, StochVector& x, SimpleVector& xp);

  virtual void solveReduced( QpGenStochData *prob, SimpleVector& b)=0;

  virtual void putXDiagonal( OoqpVector& xdiag_ );
  virtual void putZDiagonal( OoqpVector& zdiag );
 
  virtual void AddChild(sLinsys* child);

  void sync();
 public:
  virtual ~sLinsysRoot();

 public: //utilities
  void myAtPutZeros(DenseSymMatrix* mat);
  void myAtPutZeros(DenseSymMatrix* mat, 
		    int row, int col, 
		    int rowExtent, int colExtent);
  void submatrixAllReduce(DenseSymMatrix* A, 
			  int row, int col, int drow, int dcol,
			  MPI_Comm comm);
 protected: //buffers

  OoqpVector* zDiag;
  OoqpVector* xDiag;

#ifdef STOCH_TESTING
 protected: 
  static void dumpRhs(int proc, const char* nameToken,  SimpleVector& rhs);
  static void dumpMatrix(int scen, int proc, const char* nameToken, DenseSymMatrix& M);
#endif
#ifdef TIMING
 protected:
  void afterFactor();
#endif
};

#endif

