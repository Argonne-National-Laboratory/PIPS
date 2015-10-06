/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef SROOTLINSYS_NLP
#define SROOTLINSYS_NLP

#include "sLinsys.h"

class sFactory;
class sData;
class RegularizationAlg;

// DEBUG only
//#include "ScaDenSymMatrix.h"

/** 
 * ROOT (= NON-leaf) linear system
 */
class sLinsysRoot : public sLinsys {
 protected:
  sLinsysRoot() { };

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

  sLinsysRoot(sFactory * factory_, sData * prob_, bool createChild = true);
  sLinsysRoot(sFactory* factory,
	      sData* prob_,				    
	      OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_,
	      OoqpVector* rhs_, OoqpVector* additiveDiag_, bool createChild = true);

  virtual int factor2(sData *prob, Variables *vars);
  	
  /* Atoms methods of FACTOR2 for a non-leaf linear system */
  virtual void initializeKKT(sData* prob, Variables* vars);
  virtual void reduceKKT();
  virtual int  factorizeKKT(); 
  virtual void finalizeKKT(sData* prob, Variables* vars)=0;


  virtual void Lsolve ( sData *prob, OoqpVector& x );
  virtual void Dsolve ( sData *prob, OoqpVector& x );
  virtual void Ltsolve( sData *prob, OoqpVector& x );

  virtual void Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp);

  virtual void solveReduced( sData *prob, SimpleVector& b)=0;

  virtual void putXDiagonal( OoqpVector& xdiag_ );
  virtual void putSDiagonal( OoqpVector& sdiag_ );
  virtual void putYDualDiagonal( OoqpVector& ydiag );
  virtual void putZDiagonal( OoqpVector& zdiag );

  virtual void setAdditiveDiagonal( );

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
  OoqpVector* yDiag;
  OoqpVector* sDiag;

#ifdef STOCH_TESTING
 protected: 
  static void dumpRhs(int proc, const char* nameToken,  SimpleVector& rhs);
  static void dumpMatrix(int scen, int proc, const char* nameToken, DenseSymMatrix& M);
#endif
#ifdef TIMING
 protected:
  void afterFactor();
#endif

  public: 
   virtual void UpdateMatrices( Data * prob_in,int const updateLevel=2);

   bool firstQUpdate,firstBUpdate;

   std::map<int,int> LocQMap;
   std::map<int,int> LocBMap;

   virtual void setXDiagonal( OoqpVector& xdiag );
   virtual void setSDiagonal( OoqpVector& sdiag );  
   virtual void setYDiagonal( OoqpVector& ydiag );  
   virtual void setZDiagonal( OoqpVector& zdiag );

};

#endif

