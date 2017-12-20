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
 * Base class for QP scalers.
 */
class QpScaler : Scaler
{
  enum MatrixType {MATRIXTYPE_A, MATRIXTYPE_C};
protected:

  OoqpVector* vec_rowscaleA;
  OoqpVector* vec_colscaleA;
  OoqpVector* vec_rowscaleC;
  OoqpVector* vec_colscaleC;
  double factor_objscale;

  //todo empty method for scaling Q

  virtual void computeRowscaleVec(GenMatrix* matrix, MatrixType type);
  virtual void computeColscaleVec(GenMatrix* matrix, MatrixType type);

  virtual void scaleMatrix(GenMatrix* matrix);
  virtual void scaleVector(OoqpVector* vec, const OoqpVector* scalevec);
  virtual void scaleVector(OoqpVector* vec, double scalefactor);
public:

  QpScaler();
  virtual ~QpScaler();

  /** scale */
  virtual void scale( ProblemFormulation * formulation,
        QpGenData * prob, QpGenVars * vars, QpGenResiduals * resid);

  /** unscale */
  virtual void unscale( ProblemFormulation * formulation,
        QpGenData * prob, QpGenVars * vars, QpGenResiduals * resid);

};

//@}


#endif /* PIPS_IPM_CORE_QPSCALERS_QPSCALER_H_ */
