/*
 * Options.h
 *
 *  Created on: 03.04.2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_ABSTRACT_OPTIONS_H_
#define PIPS_IPM_CORE_ABSTRACT_OPTIONS_H_

#include <string>
#include <map>

#include "pipsport.h"
#include "Singleton.h"

/**
 * Abstract base class for options class.
 */

class Options : public Singleton
{
public:
   int getIntParam(const std::string& identifier) const;
   double getDoubleParam(const std::string& identifier) const;
   bool getBoolParam(const std::string& identifier) const;

   void setIntParam(const std::string& identifier, int value);
   void setDoubleParam(const std::string& identifier, double value);
   void setBoolParam(const std::string& identifier, bool value);

   void fillParamsFromFile(const std::string& filename);

private:
   std::map<std::string, double> double_options;
   std::map<std::string, int> int_options;
   std::map<std::string, bool> bool_options;

   virtual void setDefaults() = 0;
protected:
   Options() {};
   virtual ~Options() {};
};

#endif /* PIPS_IPM_CORE_ABSTRACT_OPTIONS_H_ */
