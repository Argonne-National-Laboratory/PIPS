/*
 * StochScaler.h
 *
 *  Created on: 18.07.2019
 *      Author: Nils-Christian Kempke
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHSCALER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHSCALER_H_

#include "QpScaler.h"

class Data;

class StochScaler : public QpScaler
{
public:
	StochScaler(Data* prob, bool bitshifting);
	virtual ~StochScaler();

	virtual void scale() = 0;

	virtual Variables* getUnscaledVariables(const Variables& vars) const;
  	virtual Residuals* getUnscaledResiduals(const Residuals& resids) const;

protected:
	virtual void doObjScaling() = 0;
};













































#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHSCALER_H_ */
