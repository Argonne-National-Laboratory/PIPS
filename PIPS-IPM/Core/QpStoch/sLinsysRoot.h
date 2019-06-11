/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#ifndef SROOTLINSYS
#define SROOTLINSYS

#include "sLinsys.h"
#include "StochGenMatrix.h"


class sFactory;
class sData;

// DEBUG only
//#include "ScaDenSymMatrix.h"

/** 
 * ROOT (= NON-leaf) linear system
 */
class sLinsysRoot : public sLinsys {
 protected:
  sLinsysRoot() {};

  virtual void         createChildren(sData* prob);
  virtual void         deleteChildren();

  virtual SymMatrix*   createKKT     (sData* prob) = 0;
  virtual DoubleLinearSolver* 
                       createSolver  (sData* prob, 
				      SymMatrix* kktmat) = 0;
 public:
  std::vector<sLinsys*> children;
  int iAmDistrib;

 public:

  sLinsysRoot(sFactory * factory_, sData * prob_);
  sLinsysRoot(sFactory* factory,
	      sData* prob_,				    
	      OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_,
	      OoqpVector* rhs_);

  virtual void factor2(sData *prob, Variables *vars);
  /* Atoms methods of FACTOR2 for a non-leaf linear system */
  virtual void initializeKKT(sData* prob, Variables* vars);
  virtual void reduceKKT();
  virtual void factorizeKKT(); 
  virtual void finalizeKKT(sData* prob, Variables* vars)=0;


  virtual void Lsolve ( sData *prob, OoqpVector& x );
  virtual void Dsolve ( sData *prob, OoqpVector& x );
  virtual void Ltsolve( sData *prob, OoqpVector& x );

  virtual void Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp);

  virtual void solveReduced( sData *prob, SimpleVector& b)=0;
  virtual void solveReducedLinkCons( sData *prob, SimpleVector& b) {assert("not implemented here \n" && 0);};

  virtual void putXDiagonal( OoqpVector& xdiag_ );
  virtual void putZDiagonal( OoqpVector& zdiag );
 
  virtual void AddChild(sLinsys* child);

  virtual bool usingSparseKkt() {return hasSparseKkt;};

  void sync();
 public:
  virtual ~sLinsysRoot();

 public: //utilities
  void myAtPutZeros(DenseSymMatrix* mat);
  void myAtPutZeros(DenseSymMatrix* mat, 
		    int row, int col, 
		    int rowExtent, int colExtent);

  // all_reduces specified submatrix (in chunks)
  void submatrixAllReduce(DenseSymMatrix* A, 
			  int startRow, int startCol, int nRows, int nCols,
			  MPI_Comm comm);

  // all_reduces specified submatrix as a while
  void submatrixAllReduceFull(DenseSymMatrix* A,
           int startRow, int startCol, int nRows, int nCols,
           MPI_Comm comm);

  // all_reducees lower half (including diagonal) of specified submatrix
  void submatrixAllReduceDiagLower(DenseSymMatrix* A,
            int substart, int subsize,
            MPI_Comm comm);

 protected: //buffers

  OoqpVector* zDiag;
  OoqpVector* zDiagLinkCons;
  OoqpVector* xDiag;

  double* sparseKktBuffer;

  bool hasSparseKkt;
  bool usePrecondDist;

 private:
  void addTermToSchurCompl(sData* prob, size_t childindex);
  void reduceKKTdense();
  void reduceKKTsparse();

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

