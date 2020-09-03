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

  long long nxupp;
  OoqpVectorHandle ixupp;

  long long nxlow;
  OoqpVectorHandle ixlow;

  long long mcupp;
  OoqpVectorHandle icupp;

  long long mclow;
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

  QpGenResiduals( const QpGenResiduals& res);

  const long long& getNxupp() { return nxupp; };
  const long long& getNxlow() { return nxlow; };
  const long long& getMcupp() { return mcupp; };
  const long long& getMclow() { return mclow; };

  virtual ~QpGenResiduals();
  
  void calcresids(Data *problem, Variables *vars, bool print_resids = false) override;

  void add_r3_xz_alpha(const Variables *vars, double alpha) override;

  double recomputeResidualNorm() override;

  void set_r3_xz_alpha(const Variables *vars, double alpha) override;
  
  void clear_r3() override;
  
  void clear_r1r2() override;

  void project_r3(double rmin, double rmax) override;

  virtual int  validNonZeroPattern();
  
  virtual void writeToStream(std::ostream& out);
};

#endif





