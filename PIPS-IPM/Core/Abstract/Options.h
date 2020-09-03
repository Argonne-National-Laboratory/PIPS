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
 * base class for options class.
 */

namespace base_options
{
   int getIntParameter(const std::string& identifier);
   double getDoubleParameter(const std::string& identifier);
   bool getBoolParameter(const std::string& identifier);

   class Options : public Singleton
   {
   private:
      virtual void setDefaults(){};

   protected:
      // not thread safe when modified..
      // TODO : there is no hash_map in C++03 I think
      static std::map<std::string, double> double_options;
      static std::map<std::string, int> int_options;
      static std::map<std::string, bool> bool_options;

      friend int base_options::getIntParameter(const std::string& identifier);
      friend double base_options::getDoubleParameter(const std::string& identifier);
      friend bool base_options::getBoolParameter(const std::string& identifier);

      Options();
      virtual ~Options() {};

      static Options& getInstance()
      {
         static Options opt;
         return opt;
      }

      bool isIdentifierUnique( const std::string& identifier ) const;
      bool identifierExists( const std::string& identifier ) const;
      void fillOptionsFromFile(const std::string& filename);

      int getIntParam(const std::string& identifier) const;
      double getDoubleParam(const std::string& identifier) const;
      bool getBoolParam(const std::string& identifier) const;

      void setIntParam(const std::string& param, int value);
      void setBoolParam(const std::string& param, int value);
      void setDoubleParam(const std::string& param, int value);
   };

   inline int getIntParameter(const std::string& identifier)
   {
      return Options::getInstance().getIntParam(identifier);
   }

   inline bool getBoolParameter(const std::string& identifier)
   {
      return Options::getInstance().getBoolParam(identifier);
   }

   inline double getDoubleParameter(const std::string& identifier)
   {
      return Options::getInstance().getDoubleParam(identifier);
   }
}

#endif /* PIPS_IPM_CORE_ABSTRACT_OPTIONS_H_ */
