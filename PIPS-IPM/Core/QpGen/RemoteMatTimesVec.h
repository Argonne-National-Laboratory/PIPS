#ifndef REMOTEMATVEC
#define REMOTEMATVEC

#include "DoubleLinearSolver.h"
#include "StochTreePrecond.h"

#define PREC_EXTRA_DOUBLES 1

class OoqpVector;

class RemoteMatTimesVec : public MatTimesVec {
 public:
  RemoteMatTimesVec(StochTreePrecond* tree);
  virtual ~RemoteMatTimesVec();
  /** y = beta * y + alpha * A * x */
  virtual void doIt(double beta, OoqpVector& y, double alpha, OoqpVector& x);

  void setTmpVec1(double* buf);
 protected:
  StochTreePrecond* stochNode;
  double* tmpVec1; int delTmpVec1;
 private:
  RemoteMatTimesVec();
};


#endif
