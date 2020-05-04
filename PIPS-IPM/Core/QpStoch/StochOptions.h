/*
 * StochOptions.h
 *
 *  Created on: 03.04.2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_
#define PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_

#include "Options.h"
#include "pipsport.h"

#include <cassert>

/**
 * The getInstanceMethod must be specified in this BaseClasse.
 * Defines default options for StochPIPS.
 */

namespace pips_options
{
   void setOptions(std::string opt_file);
   int getIntParameter(std::string identifier);
   double getDoubleParameter(std::string identifier);
   bool getBoolParameter(std::string identifier);

class StochOptions : public Options
{

private:
   friend void pips_options::setOptions(std::string opt_file);
   friend int pips_options::getIntParameter(std::string identifier);
   friend double pips_options::getDoubleParameter(std::string identifier);
   friend bool pips_options::getBoolParameter(std::string identifier);

   static StochOptions& getInstance()
   {
      static StochOptions opt;
      return opt;
   }

   void setDefaults() override;
   StochOptions();

   virtual ~StochOptions() {};
};

   inline void setOptions(std::string opt_file)
   {
      return StochOptions::getInstance().fillOptionsFromFile(opt_file);
   }

   inline int getIntParameter(std::string identifier)
   {
      return StochOptions::getInstance().getIntParam(identifier);
   }

   inline bool getBoolParameter(std::string identifier)
   {
      return StochOptions::getInstance().getBoolParam(identifier);
   }

   inline double getDoubleParameter(std::string identifier)
   {
      return StochOptions::getInstance().getDoubleParam(identifier);
   }
}

#endif /* PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_ */
