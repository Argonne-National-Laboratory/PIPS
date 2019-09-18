/*
 * StochScaler.C
 *
 *  Created on: 18.07.2019
 *      Author: Nils-Christian
 */

#include "StochScaler.h"
#include "pipsdef.h"
#include "sVars.h"
#include "sResiduals.h"


StochScaler::StochScaler(Data* prob, bool bitshifting)
	: QpScaler(prob, bitshifting)
{

}

StochScaler::~StochScaler()
{
}

Variables* StochScaler::getVariablesUnscaled(const Variables& vars) const
{
   Variables* s_vars = new sVars(dynamic_cast<const sVars&>(vars)); 
   assert(s_vars);
   assert( dynamic_cast<sVars*>(s_vars)->x);

   unscaleVars(*s_vars);

   return s_vars;
};

Residuals* StochScaler::getResidualsUnscaled(const Residuals& resids) const
{
   Residuals* s_resids = new sResiduals(dynamic_cast<const sResiduals&>(resids));
   assert(s_resids);

   unscaleResids(*s_resids);

   return s_resids;
};
