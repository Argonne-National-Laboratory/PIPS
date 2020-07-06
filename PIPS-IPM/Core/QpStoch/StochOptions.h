/*
 * StochOptions.h
 *
 *  Created on: 03.04.2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_
#define PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_

#include "QpGenOptions.h"
#include "pipsport.h"

#include <cassert>

/**
 * The getInstanceMethod must be specified in this BaseClass.
 * Defines default options for StochPIPS.
 */

namespace pips_options
{
   void setOptions(const std::string& opt_file);
   int getIntParameter(const std::string& identifier);
   double getDoubleParameter(const std::string& identifier);
   bool getBoolParameter(const std::string& identifier);

   class StochOptions : public qpgen_options::QpGenOptions
   {

   private:
      friend void pips_options::setOptions(const std::string& opt_file);
      friend int pips_options::getIntParameter(const std::string& identifier);
      friend double pips_options::getDoubleParameter(const std::string& identifier);
      friend bool pips_options::getBoolParameter(const std::string& identifier);

      static StochOptions& getInstance()
      {
         static StochOptions opt;
         return opt;
      }

      void setDefaults() override;
      StochOptions();

      virtual ~StochOptions() {};
   };

   inline void setOptions(const std::string& opt_file)
   {
      return StochOptions::getInstance().fillOptionsFromFile(opt_file);
   }

   inline int getIntParameter(const std::string& identifier)
   {
      return StochOptions::getInstance().getIntParam(identifier);
   }

   inline bool getBoolParameter(const std::string& identifier)
   {
      return StochOptions::getInstance().getBoolParam(identifier);
   }

   inline double getDoubleParameter(const std::string& identifier)
   {
      return StochOptions::getInstance().getDoubleParam(identifier);
   }
}

#endif /* PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_ */
