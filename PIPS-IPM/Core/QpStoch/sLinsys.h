/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#ifndef SLINSYS
#define SLINSYS

#include "QpGenLinsys.h"
#include "DoubleLinearSolver.h"
#include "OoqpVectorHandle.h"
#include "DenseSymMatrix.h"
#include "SparseSymMatrix.h"
#include "DenseGenMatrix.h"
#include "SimpleVector.h"
#include "StochVector.h"

#include <vector>

#include "mpi.h"

class sTree;
class sFactory;
class sData;

class sLinsys : public QpGenLinsys
{
 public:
  sLinsys(sFactory* factory, sData* prob);
  sLinsys(sFactory* factory,
		   sData* prob, 
		   OoqpVector* dd, 
		   OoqpVector* dq,
		   OoqpVector* nomegaInv,
		   OoqpVector* rhs);

  virtual ~sLinsys();

  virtual void factor (Data *prob, Variables *vars);
  virtual void factor2(sData *prob, Variables *vars)=0;

  virtual void Lsolve ( sData *prob, OoqpVector& x )=0;
  virtual void Dsolve ( sData *prob, OoqpVector& x )=0;
  virtual void Ltsolve( sData *prob, OoqpVector& x )=0;
  virtual void Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp)=0;

  virtual void putZDiagonal( OoqpVector& zdiag )=0;
  virtual void solveCompressed( OoqpVector& rhs );
  virtual void putXDiagonal( OoqpVector& xdiag_ )=0;

  void joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
		OoqpVector& rhs2_in, OoqpVector& rhs3_in );

  void separateVars( OoqpVector& x_in, OoqpVector& y_in,
		     OoqpVector& z_in, OoqpVector& vars_in );
  
  virtual void sync()=0;
  virtual void deleteChildren()=0;

  virtual bool isDummy() const { return false; };

 protected:
  sLinsys(){};

  SymMatrix* kkt;
  DoubleLinearSolver* solver;
  int locnx, locmy, locmyl, locmz, locmzl;
  sData* data;
  
 public:
  virtual void addLnizi(sData *prob, OoqpVector& z0, OoqpVector& zi);

  virtual void addLniziLinkCons(sData *prob, OoqpVector& z0, OoqpVector& zi, int parentmy, int parentmz);

  /** y += alpha * Lni^T * x */
  void LniTransMult(sData *prob, 
		    SimpleVector& y, 
		    double alpha, SimpleVector& x);

  /** Methods that use dense matrices U and V to compute the 
   *  terms from the Schur complement.
   */
  virtual void allocU(DenseGenMatrix ** Ut, int np);
  virtual void allocV (DenseGenMatrix ** V, int np);
  virtual void computeU_V(sData *prob, DenseGenMatrix* U, DenseGenMatrix* V);

  /** Method(s) that use a memory-friendly mechanism for computing
   *  the terms from the Schur Complement 
   */
  virtual void addTermToDenseSchurCompl(sData *prob, 
					DenseSymMatrix& SC);

  virtual void addTermToSchurComplBlocked(sData *prob, bool sparseSC,
               SymMatrix& SC);

  virtual void addTermToSparseSchurCompl(sData *prob,
               SparseSymMatrix& SC) { assert(0 && "not implemented here"); };
					
  virtual void addColsToDenseSchurCompl(sData *prob, 
					DenseGenMatrix& out, 
					int startcol, int endcol);

  virtual void symAddColsToDenseSchurCompl(sData *prob, 
				       double *out, 
				       int startcol, int endcol);
  /** Used in the iterative refinement for the dense Schur complement systems
   * Computes res += [0 A^T C^T ]*inv(KKT)*[0;A;C] x
   */
  virtual void addTermToSchurResidual(sData* prob, 
				      SimpleVector& res, 
				      SimpleVector& x);
 public:
  MPI_Comm mpiComm;
  sTree* stochNode;

 protected:
  int nThreads;

  /* members for blockwise schur complement computation */
  int blocksizemax;
  double* colsBlockDense;
  int* colId;
  int* colSparsity;

  void multLeftSchurComplBlocked(/*const*/sData *prob, /*const*/double* colsBlockDense,
        const int* colId, int blocksize, bool sparseSC, SymMatrix& SC);

  void multLeftSparseSchurComplBlocked(/*const*/sData *prob, /*const*/double* colsBlockDense,
        const int* colId, int blocksize, SparseSymMatrix& SC);

  void multLeftDenseSchurComplBlocked(/*const*/sData *prob, /*const*/double* colsBlockDense,
        const int* colId, int blocksize, DenseSymMatrix& SC);

};

#endif
