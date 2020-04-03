/*
 * Options.h
 *
 *  Created on: 03.04.2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_
#define PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_

#include <string>

#include "Options.h"
#include "pipsport.h"
/**
 * Abstract base class for options class.
 * Implements something close to a singleton pattern.
 * The getInstanceMethod must be specified in the BaseClasses.
 */


class StochOptions : public Options
{
public:
   static const StochOptions& getInstance()
   {
      static StochOptions options;
      return options;
   }
private:
   void setDefaults() override;
   StochOptions();

   virtual ~StochOptions() {};
};

#endif /* PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_ */
