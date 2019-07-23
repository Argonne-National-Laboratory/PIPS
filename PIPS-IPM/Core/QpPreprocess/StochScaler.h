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

	virtual void scale() override = 0;

	virtual Variables* getVariablesUnscaled(const Variables& vars) const override;
  	virtual Residuals* getResidualsUnscaled(const Residuals& resids) const override;

protected:
	virtual void doObjScaling() override = 0;
};













































#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHSCALER_H_ */
