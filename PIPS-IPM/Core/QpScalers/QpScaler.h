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

//class Data;


/**  * @defgroup QpScaler
 *
 * QP scaler
 * @{
 */

/**
 * Abstract base class for QP scalers.
 */
class QpScaler : Scaler
{
  enum MatrixType {MATRIXTYPE_A, MATRIXTYPE_C};
protected:

  // scaling vector
  OoqpVector* vec_rowscaleA;
  OoqpVector* vec_rowscaleC;
  OoqpVector* vec_colscale;

  // problem data
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

  //todo empty method for scaling Q

  virtual void computeScalingVecs();
  virtual void applyScaling(QpGenData * prob);

  virtual void computeRowscaleVec(GenMatrix* matrix, MatrixType type);
  virtual void computeColscaleVec(GenMatrix* matrixA, GenMatrix* matrixC);
  virtual void computeObjScaleFactor();

  virtual void scaleMatrix(GenMatrix* matrix);
  virtual void scaleVector(OoqpVector* vec, const OoqpVector* scalevec);
  virtual void scaleVector(OoqpVector* vec, double scalefactor);
public:

  QpScaler(QpGenData * prob);
  virtual ~QpScaler();

  /** scale */
  virtual void scale( ProblemFormulation * formulation,
        QpGenData * prob, QpGenVars * vars, QpGenResiduals * resid) = 0;

  /** unscale */
  virtual void unscale( ProblemFormulation * formulation,
        QpGenData * prob, QpGenVars * vars, QpGenResiduals * resid) = 0;

};

//@}


#endif /* PIPS_IPM_CORE_QPSCALERS_QPSCALER_H_ */
