/*
 * GeoStochScaler.h
 *
 *  Created on: 16.07.2018
 *      Author: Svenja Uslu
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_GEOSTOCHSCALER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_GEOSTOCHSCALER_H_

#include "../QpPreprocess/StochScaler.h"
#include "StochVector.h"

class Data;


/**  * @defgroup QpPreprocess
 *
 * Geometric scaler
 * @{
 */

/**
 * Derived class for Geometric scaler.
 */
class GeoStochScaler : public StochScaler
{
   bool equilibrate;
   int maxIters;
   double minImpr;
   double goodEnough;

protected:
  virtual void doObjScaling();
  void applyGeoMean(StochVector& maxvec, StochVector& minvec);
  void postEquiScale();
  void setScalingVecsToOne();

public:

  GeoStochScaler(Data * prob, bool equiScaling, bool bitshifting = false);
  virtual ~GeoStochScaler();

  /** scale */
  virtual void scale();
};

//@}

#endif /* PIPS_IPM_CORE_QPPREPROCESS_GEOSTOCHSCALER_H_ */
