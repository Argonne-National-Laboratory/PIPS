/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#ifndef SLINSYS
#define SLINSYS

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
class sFactory;
class QpGenStochData;

class sLinsys : public QpGenLinsys
{
 public:
  sLinsys(sFactory* factory, QpGenStochData* prob);
  sLinsys(sFactory* factory,
		   QpGenStochData* prob, 
		   OoqpVector* dd, 
		   OoqpVector* dq,
		   OoqpVector* nomegaInv,
		   OoqpVector* rhs);

  virtual ~sLinsys();

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
  sLinsys(){};
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


  /** Methods that use dense matrices U and V to compute the 
   *  terms from the Schur complement.
   */
  virtual void allocU(DenseGenMatrix ** Ut, int np);
  virtual void allocV (DenseGenMatrix ** V, int np);
  virtual void computeU_V(QpGenStochData *prob, DenseGenMatrix* U, DenseGenMatrix* V);

  /** Method(s) that use a memory-friendly mechanism for computing
   *  the terms from the Schur Complement 
   */
  virtual void addTermToDenseSchurCompl(QpGenStochData *prob, 
					DenseSymMatrix& SC);
					
  virtual void addColsToDenseSchurCompl(QpGenStochData *prob, 
					DenseGenMatrix& out, 
					int startcol, int endcol);

  virtual void symAddColsToDenseSchurCompl(QpGenStochData *prob, 
				       double *out, 
				       int startcol, int endcol);

 public:
  MPI_Comm mpiComm;
  StochTree* stochNode;
};

#endif
