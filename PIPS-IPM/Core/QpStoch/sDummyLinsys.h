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

  virtual void factor2( sData *prob, Variables *vars){};
  virtual void Lsolve ( sData *prob, OoqpVector& x ){};
  virtual void Dsolve ( sData *prob, OoqpVector& x ){};
  virtual void Ltsolve( sData *prob, OoqpVector& x ){};

  virtual void Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp){};

  virtual void putZDiagonal( OoqpVector& zdiag ){};
  virtual void solveCompressed( OoqpVector& rhs ){};
  virtual void putXDiagonal( OoqpVector& xdiag_ ){};

  void joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
		OoqpVector& rhs2_in, OoqpVector& rhs3_in ){};

  void separateVars( OoqpVector& x_in, OoqpVector& y_in,
		     OoqpVector& z_in, OoqpVector& vars_in ){};
  

  virtual void addLnizi(sData *prob, OoqpVector& z0, OoqpVector& zi){};
  virtual void addLniziLinkCons(sData *prob, OoqpVector& z0, OoqpVector& zi, int my0, int mz0){};

  /** y += alpha * Lni^T * x */
  void LniTransMult(sData *prob, 
		    SimpleVector& y, 
		    double alpha, SimpleVector& x){};

  void addTermToSchurResidual(sData* prob, 
			      SimpleVector& res, 
			      SimpleVector& x) {};


  virtual void allocU(DenseGenMatrix ** Ut, int np){};
  virtual void allocV (DenseGenMatrix ** V, int np){};
  virtual void computeU_V(sData *prob, DenseGenMatrix* U, DenseGenMatrix* V){};
  void sync(){};
  virtual void deleteChildren(){};
}; 

#endif
