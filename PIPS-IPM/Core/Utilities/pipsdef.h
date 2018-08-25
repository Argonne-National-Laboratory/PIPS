/*
 * pipsdef.h
 *
 *  Created on: 29.01.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_UTILITIES_PIPSDEF_H_
#define PIPS_IPM_CORE_UTILITIES_PIPSDEF_H_

#include <sstream>
#include <iostream>
#include <cmath>
#include <string>

const double pips_eps = 1e-13;
const double pips_eps0 = 1e-40;

#ifdef PIPS_DEBUG
#define PIPSdebugMessage                printf("[%s:%d] debug: ", __FILE__, __LINE__), printf
#else
#define PIPSdebugMessage      while( 0 ) /*lint -e{530}*/ printf
#endif

inline bool PIPSisEQ(double val1, double val2, double eps = pips_eps)
{
   return (std::fabs(val1 - val2) <= pips_eps);
}

inline bool PIPSisLE(double val1, double val2, double eps = pips_eps)
{
   return (val1 <= val2 + pips_eps);
}

inline bool PIPSisLT(double val1, double val2, double eps = pips_eps)
{
   return (val1 < val2 - pips_eps);
}

inline bool PIPSisZero(double val, double eps0 = pips_eps0)
{
   return (std::fabs(val) < pips_eps0);
}

#endif

