/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef QPGENRESIDUALS
#define QPGENRESIDUALS

#include <iostream>
#include <fstream>

#include "Residuals.h"
#include "OoqpVectorHandle.h"

class QpGen;
class QpGenData;
class Variables;
class LinearAlgebraPackage;

/** 
 * Residuals for the general QP formulation 
 *
 * @ingroup QpGen
 */

class QpGenResiduals : public Residuals {
protected:
  long long nx, my, mz;

  double nxupp;
  OoqpVectorHandle ixupp;

  double nxlow;
  OoqpVectorHandle ixlow;

  double mcupp;
  OoqpVectorHandle icupp;

  double mclow;
  OoqpVectorHandle iclow;

  QpGenResiduals() {};

public:
  OoqpVectorHandle rQ;
  OoqpVectorHandle rA;
  OoqpVectorHandle rC;
  OoqpVectorHandle rz;
  OoqpVectorHandle rv;
  OoqpVectorHandle rw;
  OoqpVectorHandle rt;
  OoqpVectorHandle ru;
  OoqpVectorHandle rgamma;
  OoqpVectorHandle rphi;
  OoqpVectorHandle rlambda;
  OoqpVectorHandle rpi;

  QpGenResiduals( LinearAlgebraPackage * la,
		  long long nx, long long my, long long mz,
		  OoqpVector * ixlow, OoqpVector * ixupp,
		  OoqpVector * iclow, OoqpVector * icupp );


  virtual void calcresids(Data *problem, Variables *vars);

  virtual void add_r3_xz_alpha(Variables *vars, double alpha);

  virtual void set_r3_xz_alpha(Variables *vars, double alpha);
  
  virtual void clear_r3();
  
  virtual void clear_r1r2();

  virtual void project_r3(double rmin, double rmax);

  virtual int  validNonZeroPattern();
  
  virtual ~QpGenResiduals();

  virtual void writeToStream(std::ostream& out);
};

#endif





