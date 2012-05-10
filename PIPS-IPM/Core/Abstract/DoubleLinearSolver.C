#include "DoubleLinearSolver.h"
#include "SimpleVector.h"

DoubleIterativeLinearSolver::
DoubleIterativeLinearSolver( MatTimesVec* Ain, MatTimesVec* M1in, MatTimesVec* M2in )
: A(Ain), ML(M1in), MR(M2in)
{

}

DoubleIterativeLinearSolver::DoubleIterativeLinearSolver()
  : A(NULL), ML(NULL), MR(NULL)
{ }

DoubleIterativeLinearSolver::~DoubleIterativeLinearSolver()
{ }

void 
DoubleIterativeLinearSolver::diagonalChanged( int idiag, int extent )
{ }

void
DoubleIterativeLinearSolver::matrixChanged()
{ }

void DoubleIterativeLinearSolver::applyA (double beta, OoqpVector& res, 
					  double alpha, OoqpVector& x)
{
  A->doIt(beta,res,alpha,x);
}


void DoubleIterativeLinearSolver::applyM1(double beta, OoqpVector& res, 
					  double alpha, OoqpVector& x)
{
  ML->doIt(beta,res,alpha,x);
}

void DoubleIterativeLinearSolver::applyM2(double beta,  OoqpVector& res, 
					  double alpha, OoqpVector& x)
{
  if(NULL==MR) {
    //just a identity precond
    if(beta==0.0) {
      if(alpha==1.0) {res.copyFrom(x);return;}
      //alpha not zero
      res.setToZero();
      
    } else if(beta!=1.0) res.scale(beta);
    
    //beta not 0.0 and alpha not 1.0
    if(alpha!=0) res.axpy(alpha, x);
    
  } else         MR->doIt(beta,res,alpha,x);
}



/**********************************************************************
 * StoredMatTimesVec implementation
 **********************************************************************/

StoredMatTimesVec::StoredMatTimesVec(DoubleMatrix* mat)
  : mMat(mat)
{ };


void StoredMatTimesVec::doIt(double beta, OoqpVector& y_, 
			     double alpha, OoqpVector& x_)
{
  SimpleVector& x = dynamic_cast<SimpleVector&>(x_);
  SimpleVector& y = dynamic_cast<SimpleVector&>(y_);

  int m,n; mMat->getSize(m,n);
  assert(x.length() == n);
  assert(y.length() == m);
  
  mMat->mult(beta, y, alpha, x);
}
/**********************************************************************
 * StoredMatTransTimesVec implementation
 **********************************************************************/

StoredMatTransTimesVec::StoredMatTransTimesVec(DoubleMatrix* mat)
  : mMat(mat)
{ };


void StoredMatTransTimesVec::doIt(double beta, OoqpVector& y_, 
				  double alpha, OoqpVector& x_)
{
  SimpleVector& x = dynamic_cast<SimpleVector&>(x_);
  SimpleVector& y = dynamic_cast<SimpleVector&>(y_);

  int m,n; mMat->getSize(m,n);
  assert(x.length() == m);
  assert(y.length() == n);
  
  mMat->transMult(beta, y, alpha, x);
}
