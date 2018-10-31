/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef SLINSYS_NLP
#define SLINSYS_NLP

#include "NlpGenLinsys.h"
#include "DoubleLinearSolver.h"
#include "OoqpVectorHandle.h"
#include "DenseSymMatrix.h"
#include "SparseSymMatrixRowMajList.h"
#include "DenseGenMatrix.h"
#include "SimpleVector.h"
#include "StochVector.h"

#include <vector>

#include "mpi.h"

class sTree;
class sFactory;
class sData;

class sLinsys : public NlpGenLinsys
{
 public:
  sLinsys(sFactory* factory, sData* prob);
  sLinsys(sFactory* factory,
		   sData* prob, 
		   OoqpVector* dd, 
		   OoqpVector* dq,
		   OoqpVector* nomegaInv,
		   OoqpVector* rhs,
		   OoqpVector* additiveDiag_=NULL);

  virtual ~sLinsys();

  virtual void factor(Data *prob, Variables *vars,RegularizationAlg *RegInfo);

  virtual void factor(Data *prob, Variables *vars);

  virtual int factor2(sData *prob, Variables *vars)=0;

  virtual void Lsolve ( sData *prob, OoqpVector& x )=0;
  virtual void Dsolve ( sData *prob, OoqpVector& x )=0;
  virtual void Ltsolve( sData *prob, OoqpVector& x )=0;
  virtual void Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp)=0;


  virtual void solveCompressed( OoqpVector& rhs );
  virtual void putXDiagonal( OoqpVector& xdiag_ )=0;
  virtual void putSDiagonal( OoqpVector& sdiag_ )=0;
  virtual void putYDualDiagonal( OoqpVector& ydiag_ )=0;
  virtual void putZDiagonal( OoqpVector& zdiag )=0;

  virtual void setAdditiveDiagonal()=0;

  void joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
		OoqpVector& rhs2_in, OoqpVector& rhs3_in );

  void separateVars( OoqpVector& x_in, OoqpVector& y_in,
		     OoqpVector& z_in, OoqpVector& vars_in );

  void joinRHSXSYZ( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			OoqpVector& rhs2_in, OoqpVector& rhs3_in, OoqpVector& rhs4_in );
  
  void separateVarsXSYZ( OoqpVector& x_in, OoqpVector& s_in, 
			OoqpVector& y_in, OoqpVector& z_in, OoqpVector& vars_in);  
  
  
  virtual void sync()=0;
  virtual void deleteChildren()=0;

  int negEigVal;

 protected:
  sLinsys(){};

  int locnx, locmy, locmz;
  
  DoubleLinearSolver* solver;
  sData* data;
  
 public:
  SymMatrix* kkt; 	
  virtual void addLnizi(sData *prob, OoqpVector& z0, OoqpVector& zi);

  /** y += alpha * Lni^T * x */
  void LniTransMult(sData *prob, 
		    OoqpVector& y, 
		    double alpha, OoqpVector& x);


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
  virtual void addTermToDenseSchurCompl(sData *prob, 
					SparseSymMatrixRowMajList& SC);
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

  virtual int GetNegEigVal(){return solver->negEigVal;};

  virtual void _backSolve(sData *prob, OoqpVector& ParSol_ , OoqpVector& Vec_, StochVector* End_Par_Pos_);

  virtual void _addTargetParsLnizi(sData *prob, OoqpVector& ParSol_ , OoqpVector& Vec_, OoqpVector* goal_Par); 

  virtual void _setupColOfBordMat( sData *prob, OoqpVector* rhs_St, const int ColIDX, bool &allzero, 
										const int aimlevel);
  
  virtual void _assembleSC( sData *prob, OoqpVector* rhs_St_in, const int ColIDX, 
										const int aimlevel, DenseSymMatrix& SC);
  virtual void _setupYaddLniTx( sData *prob, OoqpVector& y_, 
			   double alpha, SimpleVector& x, const int aimlevel);

  virtual void initialize(sFactory* factory, sData* prob) = 0;
 public:
  MPI_Comm mpiComm;
  sTree* stochNode;

  bool isActive;
};

#endif
