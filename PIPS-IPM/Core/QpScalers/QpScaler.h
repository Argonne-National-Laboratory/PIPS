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
  enum MatrixType {MATRIXTYPE_A, MATRIXTYPE_C, MATRIXTYPE_Q};

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

  //todo method for scaling Q

  virtual void applyScaling();
  virtual void doObjScaling() = 0;
  /*
  virtual void computeScalingVecs();

  virtual void computeRowscaleVec(GenMatrix* matrix, MatrixType type) = 0;
  virtual void computeColscaleVec(GenMatrix* matrixA, GenMatrix* matrixC, GenMatrix* matrixQ) = 0;


  virtual void scaleMatrix(GenMatrix* matrix);
  virtual void scaleVector(OoqpVector* vec, const OoqpVector* scalevec);
  virtual void scaleVector(OoqpVector* vec, double scalefactor);
*/
  /** get maximum absolute row ratio and write maximum row entries into vectors */
  virtual double maxRowRatio(OoqpVector& maxvecA, OoqpVector& maxvecC, OoqpVector& minvecA, OoqpVector& minvecC);

  /** get maximum absolute column ratio and write maximum column entries into vectors */
  virtual double maxColRatio(OoqpVector& maxvec, OoqpVector& minvec);
public:

  QpScaler(Data* prob, bool bitshifting = true);
  virtual ~QpScaler();

  /** scale */
  virtual void scale() = 0;;

  /** unscale given objective value */
  virtual double getOrigObj(double objval);

};

//@}


#endif /* PIPS_IPM_CORE_QPSCALERS_QPSCALER_H_ */
