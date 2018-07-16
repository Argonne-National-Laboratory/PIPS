/*
 * GeoStochScaler.h
 *
 *  Created on: 16.07.2018
 *      Author: bzfuslus
 */

#ifndef PIPS_IPM_CORE_QPSCALERS_GEOSTOCHSCALER_H_
#define PIPS_IPM_CORE_QPSCALERS_GEOSTOCHSCALER_H_

#include "QpScaler.h"

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
protected:
  virtual void doObjScaling();

public:

  GeoStochScaler(Data * prob, bool bitshifting = true);
  virtual ~GeoStochScaler();

  /** scale */
  virtual void scale();
};

//@}

#endif /* PIPS_IPM_CORE_QPSCALERS_GEOSTOCHSCALER_H_ */
