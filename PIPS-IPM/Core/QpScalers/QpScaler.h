/*
 * QpScaler.h
 *
 *  Created on: 19.12.2017
 *      Author: Daniel Rehfeldt
 */

#ifndef PIPS_IPM_CORE_QPSCALERS_QPSCALER_H_
#define PIPS_IPM_CORE_QPSCALERS_QPSCALER_H_

#include "Scaler.h"
#include "OoqpVector.h"
#include "DoubleMatrix.h"
#include "QpGenData.h"
#include "QpGenResiduals.h"

class Data;


/**  * @defgroup QpScaler
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
  virtual void doObjScaling() = 0;

  /** get maximum absolute row ratio and write maximum row entries into vectors */
  virtual double maxRowRatio(OoqpVector& maxvecA, OoqpVector& maxvecC, OoqpVector& minvecA, OoqpVector& minvecC);

  /** get maximum absolute column ratio and write maximum column entries into vectors */
  virtual double maxColRatio(OoqpVector& maxvec, OoqpVector& minvec);
public:

  QpScaler(Data* prob, bool bitshifting = true);
  virtual ~QpScaler();

  /** scale */
  virtual void scale() = 0;

  /** unscale given objective value */
  virtual double getOrigObj(double objval);

};

//@}


#endif /* PIPS_IPM_CORE_QPSCALERS_QPSCALER_H_ */
