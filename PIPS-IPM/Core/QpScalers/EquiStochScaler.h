/*
 * EquiStochScaler.h
 *
 *  Created on: 20.12.2017
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPSCALERS_EQUISTOCHSCALER_H_
#define PIPS_IPM_CORE_QPSCALERS_EQUISTOCHSCALER_H_

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
class EquiStochScaler : public QpScaler
{
protected:
  virtual void doObjScaling();

public:

  EquiStochScaler(Data * prob, bool bitshifting = true);
  virtual ~EquiStochScaler();

  /** scale */
  virtual void scale();
};

//@}


#endif /* PIPS_IPM_CORE_QPSCALERS_EQUISTOCHSCALER_H_ */
