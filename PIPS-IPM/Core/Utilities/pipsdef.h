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

const double pips_eps = 10e-14;

#ifdef PIPS_DEBUG
#define PIPSdebugMessage                printf("[%s:%d] debug: ", __FILE__, __LINE__), printf
#else
#define PIPSdebugMessage      while( 0 ) /*lint -e{530}*/ printf
#endif

inline bool PIPSisEQ(double val1, double val2)
{
   return (std::fabs(val1 - val2) <= pips_eps);
}

inline bool PIPSisLE(double val1, double val2)
{
   return (val1 <= val2 + pips_eps);
}

inline bool PIPSisLT(double val1, double val2)
{
   return (val1 < val2 - pips_eps);
}

inline bool PIPSisEQ_withTolerance(double val1, double val2, double tolerance)
{
   return (std::fabs(val1 - val2) <= tolerance);
}

#endif

