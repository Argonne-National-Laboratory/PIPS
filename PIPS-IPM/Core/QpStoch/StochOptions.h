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
   void setIntParameter(const std::string& identifier, int value);
   void setDoubleParameter(const std::string& identifier, double value);
   void setBoolParameter(const std::string& identifier, bool value);

   int getIntParameter(const std::string& identifier);
   double getDoubleParameter(const std::string& identifier);
   bool getBoolParameter(const std::string& identifier);

   class StochOptions : public qpgen_options::QpGenOptions
   {

   private:
      friend void setIntParameter(const std::string& identifier, int value);
      friend void setDoubleParameter(const std::string& identifier, double value);
      friend void setBoolParameter(const std::string& identifier, bool value);
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
      void setPresolveDefaults();
      StochOptions();

      virtual ~StochOptions() {};
   };

   inline void setOptions(const std::string& opt_file)
   {
      return StochOptions::getInstance().fillOptionsFromFile(opt_file);
   }

   inline void setIntParameter(const std::string& identifier, int value)
   {
      StochOptions::getInstance().setIntParam(identifier, value);
   }

   inline void setDoubleParameter(const std::string& identifier, double value)
   {
      StochOptions::getInstance().setIntParam(identifier, value);
   }

   inline void setBoolParameter(const std::string& identifier, bool value)
   {
      StochOptions::getInstance().setIntParam(identifier, value);
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
