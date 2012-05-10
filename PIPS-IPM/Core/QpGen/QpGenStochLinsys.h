#ifndef QPGENSTOCHLINSYS
#define QPGENSTOCHLINSYS

#include "QpGenLinsys.h"
#include "DoubleLinearSolver.h"
#include "OoqpVectorHandle.h"
#include "DenseSymMatrix.h"
#include "DenseGenMatrix.h"
#include "SimpleVector.h"
#include "StochVector.h"

#include <vector>

#include "mpi.h"

class StochTree;
class QpGenStoch;
class QpGenStochData;

class QpGenStochLinsys : public QpGenLinsys
{
 public:
  QpGenStochLinsys(QpGenStoch* factory, QpGenStochData* prob);
  QpGenStochLinsys(QpGenStoch* factory,
		   QpGenStochData* prob, 
		   OoqpVector* dd, 
		   OoqpVector* dq,
		   OoqpVector* nomegaInv,
		   OoqpVector* rhs);

  virtual ~QpGenStochLinsys();

  virtual void factor (Data *prob, Variables *vars);
  virtual void factor2(QpGenStochData *prob, Variables *vars)=0;

  virtual void Lsolve ( QpGenStochData *prob, OoqpVector& x )=0;
  virtual void Dsolve ( QpGenStochData *prob, OoqpVector& x )=0;
  virtual void Ltsolve( QpGenStochData *prob, OoqpVector& x )=0;
  virtual void Ltsolve2( QpGenStochData *prob, StochVector& x, SimpleVector& xp)=0;

  virtual void putZDiagonal( OoqpVector& zdiag )=0;
  virtual void solveCompressed( OoqpVector& rhs );
  virtual void putXDiagonal( OoqpVector& xdiag_ )=0;

  void joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
		OoqpVector& rhs2_in, OoqpVector& rhs3_in );

  void separateVars( OoqpVector& x_in, OoqpVector& y_in,
		     OoqpVector& z_in, OoqpVector& vars_in );
  
  virtual void sync()=0;
  virtual void deleteChildren()=0;
 protected:
  QpGenStochLinsys(){};
  SymMatrix* kkt;
  DoubleLinearSolver* solver;
  int locnx, locmy, locmz;
  QpGenStochData* data;
  
 public:
  virtual void addLnizi(QpGenStochData *prob, OoqpVector& z0, OoqpVector& zi);

  /** y += alpha * Lni^T * x */
  void LniTransMult(QpGenStochData *prob, 
		    SimpleVector& y, 
		    double alpha, SimpleVector& x);

  virtual void allocU(DenseGenMatrix ** Ut, int np);
  virtual void allocV (DenseGenMatrix ** V, int np);
  virtual void computeU_V(QpGenStochData *prob, DenseGenMatrix* U, DenseGenMatrix* V);

 public:
  MPI_Comm mpiComm;
  StochTree* stochNode;
};


/** 
 * DUMMY Linear system class
 */
class QpGenStochDummyLinsys : public QpGenStochLinsys
{
 public:
  QpGenStochDummyLinsys(QpGenStoch* factory, QpGenStochData* prob)
    : QpGenStochLinsys(factory, prob, NULL, NULL, NULL, NULL) 
    {
      mpiComm = MPI_COMM_NULL;
    };


  virtual ~QpGenStochDummyLinsys(){};

  virtual void factor2( QpGenStochData *prob, Variables *vars){};
  virtual void Lsolve ( QpGenStochData *prob, OoqpVector& x ){};
  virtual void Dsolve ( QpGenStochData *prob, OoqpVector& x ){};
  virtual void Ltsolve( QpGenStochData *prob, OoqpVector& x ){};

  virtual void Ltsolve2( QpGenStochData *prob, StochVector& x, SimpleVector& xp){};

  virtual void putZDiagonal( OoqpVector& zdiag ){};
  virtual void solveCompressed( OoqpVector& rhs ){};
  virtual void putXDiagonal( OoqpVector& xdiag_ ){};

  void joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
		OoqpVector& rhs2_in, OoqpVector& rhs3_in ){};

  void separateVars( OoqpVector& x_in, OoqpVector& y_in,
		     OoqpVector& z_in, OoqpVector& vars_in ){};
  

  virtual void addLnizi(QpGenStochData *prob, OoqpVector& z0, OoqpVector& zi){};

  /** y += alpha * Lni^T * x */
  void LniTransMult(QpGenStochData *prob, 
		    SimpleVector& y, 
		    double alpha, SimpleVector& x){};

  virtual void allocU(DenseGenMatrix ** Ut, int np){};
  virtual void allocV (DenseGenMatrix ** V, int np){};
  virtual void computeU_V(QpGenStochData *prob, DenseGenMatrix* U, DenseGenMatrix* V){};
  void sync(){};
  virtual void deleteChildren(){};
}; 


#endif
