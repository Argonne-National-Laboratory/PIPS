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
 * abstract base class for options class.
 */

enum OptionType
{
   INT = 0,
   DOUBLE = 1,
   BOOL = 2
};

class Options : public Singleton
{
private:
   virtual void setDefaults() = 0;

protected:
   // TODO : there is no hash_map in C++03 I think
   std::map<std::string, double> double_options;
   std::map<std::string, int> int_options;
   std::map<std::string, bool> bool_options;

   Options() {};
   virtual ~Options() {};

   bool isIdentifierUnique( const std::string& identifier ) const;

   void fillOptionsFromFile(const std::string& filename);

   int getIntParam(const std::string& identifier) const;
   double getDoubleParam(const std::string& identifier) const;
   bool getBoolParam(const std::string& identifier) const;
};

#endif /* PIPS_IPM_CORE_ABSTRACT_OPTIONS_H_ */
