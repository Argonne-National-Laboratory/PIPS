/*
 * QpScaler.h
 *
 *  Created on: 19.12.2017
 *      Author: Daniel Rehfeldt
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_QPSCALER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_QPSCALER_H_

#include "Scaler.h"
#include "OoqpVector.h"
#include "DoubleMatrix.h"

/**  * @defgroup QpPreprocess
 *
 * QP scaler
 * @{
 */

/**
 * Abstract base class for QP scalers.
 */
class QpScaler : public Scaler
{
protected:
  static void invertAndRound(bool round, OoqpVector& vector)
  {
     vector.invertSave(1.0);
     if( round )
        vector.roundToPow2();
  }

  // scaling vector
  OoqpVector* vec_rowscaleQ;
  OoqpVector* vec_rowscaleA;
  OoqpVector* vec_rowscaleC;
  OoqpVector* vec_colscale;

  // problem data
  SymMatrixHandle Q;
  GenMatrixHandle A;
  GenMatrixHandle C;
  OoqpVectorHandle obj;
  OoqpVectorHandle bA;
  OoqpVectorHandle bux;
  OoqpVectorHandle blx;
  OoqpVectorHandle rhsC;
  OoqpVectorHandle lhsC;

  // scaling factor for objective
  double factor_objscale;

  virtual void applyScaling();
  virtual void unscaleVars( Variables& vars ) const;
  virtual void unscaleResids( Residuals& vars ) const;

  virtual void doObjScaling() = 0;

  /** get maximum absolute row ratio and write maximum row entries into vectors */
  virtual double maxRowRatio(OoqpVector& maxvecA, OoqpVector& maxvecC, OoqpVector& minvecA, OoqpVector& minvecC, const OoqpVector* colScalevec);

  /** get maximum absolute column ratio and write maximum column entries into vectors */
  virtual double maxColRatio(OoqpVector& maxvec, OoqpVector& minvec, const OoqpVector* rowScaleVecA,  const OoqpVector* rowScaleVecC);

  void scaleObjVector(double scaling_factor);
public:

  QpScaler(Data* prob, bool bitshifting = false);
  virtual ~QpScaler();

  /** scale */
  virtual void scale() = 0;

  virtual double getOrigObj(double objval) const;
  virtual Variables* getUnscaledVariables(const Variables& vars) const;
  virtual Residuals* getUnscaledResiduals(const Residuals& resids) const;
  virtual OoqpVector* getOrigPrimal(const OoqpVector& solprimal) const;
  virtual OoqpVector* getOrigDualEq(const OoqpVector& soldual) const;
  virtual OoqpVector* getOrigDualIneq(const OoqpVector& soldual) const;
  virtual OoqpVector* getOrigDualVarBoundsUpp(const OoqpVector& soldual) const;
  virtual OoqpVector* getOrigDualVarBoundsLow(const OoqpVector& soldual) const;
};

//@}


#endif /* PIPS_IPM_CORE_QPPREPROCESS_QPSCALER_H_ */
