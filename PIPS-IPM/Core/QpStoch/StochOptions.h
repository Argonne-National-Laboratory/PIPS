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

static int PIPSgetIntParameter(const std::string& identifier);
static double PIPSgetDoubleParameter(const std::string& identifier);
static bool PIPSgetBoolParameter(const std::string& identifier);
static void PIPSfillParametersFromFile(const std::string& filename);

class StochOptions : public Options
{

private:
   static StochOptions& getInstance()
   {
      static StochOptions opt;
      return opt;
   }

   friend int PIPSgetIntParameter(const std::string& identifier);
   friend double PIPSgetDoubleParameter(const std::string& identifier);
   friend bool PIPSgetBoolParameter(const std::string& identifier);
   friend void PIPSfillParametersFromFile(const std::string& filename);

   void setDefaults() override;
   StochOptions();

   virtual ~StochOptions() {};
};

inline int PIPSgetIntParameter(const std::string& identifier)
{
   const StochOptions& opt = StochOptions::getInstance();

   assert(opt.isIdentifierUnique(identifier));

   return opt.getIntParam(identifier);
}

inline double PIPSgetDoubleParameter(const std::string& identifier)
{
   const StochOptions& opt = StochOptions::getInstance();

   assert(opt.isIdentifierUnique(identifier));

   return opt.getDoubleParam(identifier);
}

inline bool PIPSgetBoolParameter(const std::string& identifier)
{
   const StochOptions& opt = StochOptions::getInstance();

   assert(opt.isIdentifierUnique(identifier));

   return opt.getBoolParam(identifier);
}

inline void PIPSfillParametersFromFile(const std::string& filename)
{
   StochOptions::getInstance().fillParamsFromFile(filename);
}

#endif /* PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_ */
