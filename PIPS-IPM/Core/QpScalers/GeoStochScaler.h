/*
 * GeoStochScaler.h
 *
 *  Created on: 16.07.2018
 *      Author: Svenja Uslu
 */

#ifndef PIPS_IPM_CORE_QPSCALERS_GEOSTOCHSCALER_H_
#define PIPS_IPM_CORE_QPSCALERS_GEOSTOCHSCALER_H_

#include "QpScaler.h"
#include "StochVector.h"

class Data;


/**  * @defgroup QpScaler
 *
 * QP scaler
 * @{
 */

/**
 * Derived class for QP scalers.
 */
class GeoStochScaler : public QpScaler
{
   bool equilibrate;
   int maxIters;
   double minImpr;
   double goodEnough;

protected:
  virtual void doObjScaling();
  void applyGeoMean(StochVector& maxvec, StochVector& minvec);

public:

  GeoStochScaler(Data * prob, bool bitshifting = true);
  virtual ~GeoStochScaler();

  /** scale */
  virtual void scale();
};

//@}

#endif /* PIPS_IPM_CORE_QPSCALERS_GEOSTOCHSCALER_H_ */
