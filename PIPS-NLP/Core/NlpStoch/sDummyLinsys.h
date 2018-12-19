/* PIPS-IPM                                                             
 * Author: Cosmin G. Petra
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
 
/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef SDUMMYLINSYS
#define SDUMMYLINSYS

#include "sLinsys.h"

/** 
 * DUMMY Linear system class
 */
class sDummyLinsys : public sLinsys
{
 public:
  sDummyLinsys(sFactory* factory, sData* prob)
    : sLinsys(factory, prob, NULL, NULL, NULL, NULL) 
    {
      mpiComm = MPI_COMM_NULL;
    };


  virtual ~sDummyLinsys(){};

  virtual int factor2( sData *prob, Variables *vars){return 0;};
  virtual void Lsolve ( sData *prob, OoqpVector& x ){};
  virtual void Dsolve ( sData *prob, OoqpVector& x ){};
  virtual void Ltsolve( sData *prob, OoqpVector& x ){};

  virtual void Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp){};


  virtual void solveCompressed( OoqpVector& rhs ){};
  virtual void putXDiagonal( OoqpVector& xdiag_ ){};
  virtual void putSDiagonal( OoqpVector& sdiag_ ){};  
  virtual void putYDualDiagonal( OoqpVector& ydiag_ ){};
  virtual void putZDiagonal( OoqpVector& zdiag ){};
  virtual void setAdditiveDiagonal(){};

  void joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
		OoqpVector& rhs2_in, OoqpVector& rhs3_in ){};

  void separateVars( OoqpVector& x_in, OoqpVector& y_in,
		     OoqpVector& z_in, OoqpVector& vars_in ){};
  
  void joinRHSXSYZ( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
		OoqpVector& rhs2_in, OoqpVector& rhs3_in,OoqpVector& rhs4_in ){};

  void separateVarsXSYZ( OoqpVector& x_in, OoqpVector& s_in,
		     OoqpVector& y_in, OoqpVector& z_in, OoqpVector& vars_in ){};  

  virtual void addLnizi(sData *prob, OoqpVector& z0, OoqpVector& zi){};

  /** y += alpha * Lni^T * x */
  void LniTransMult(sData *prob, 
		    OoqpVector& y, 
		    double alpha, OoqpVector& x){};

  void addTermToSchurResidual(sData* prob, 
			      SimpleVector& res, 
			      SimpleVector& x) {};


  virtual void allocU(DenseGenMatrix ** Ut, int np){};
  virtual void allocV (DenseGenMatrix ** V, int np){};
  virtual void computeU_V(sData *prob, DenseGenMatrix* U, DenseGenMatrix* V){};
  void sync(){};
  virtual void deleteChildren(){};
  virtual void UpdateMatrices( Data * prob_in, int const onlyRegPart){};

  virtual void setXDiagonal( OoqpVector& xdiag ){};
  virtual void setSDiagonal( OoqpVector& sdiag ){};  
  virtual void setYDiagonal( OoqpVector& ydiag ){};  
  virtual void setZDiagonal( OoqpVector& zdiag ){};  

  virtual void initialize(sFactory* factory, sData* prob) {};

virtual void _backSolve(sData *prob, OoqpVector& ParSol_ , OoqpVector& Vec_, StochVector* End_Par_Pos_){};

virtual void _addTargetParsLnizi(sData *prob, OoqpVector& ParSol_ , OoqpVector& Vec_, OoqpVector* goal_Par){}; 

virtual void _setupColOfBordMat( sData *prob, OoqpVector* rhs_St, const int ColIDX, bool &allzero, 
									  const int aimlevel){};

virtual void _assembleSC( sData *prob, OoqpVector* rhs_St_in, const int ColIDX, 
									  const int aimlevel, DenseSymMatrix& SC){};


}; 

#endif
