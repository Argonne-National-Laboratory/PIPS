#include "Options.h"

int Options::getIntParam(const std::string& identifier) const
{
   const std::map<std::string, int>::const_iterator& it = int_options.find(identifier);

   if( it != int_options.end() )
      return it->second;
   else
      throw new std::string("No element \"" + identifier + "\" of type int in options");
}

double Options::getDoubleParam(const std::string& identifier) const
{
   const std::map<std::string, double>::const_iterator it = double_options.find(identifier);

   if( it != double_options.end() )
      return it->second;
   else
      throw new std::string("No element \"" + identifier + "\" of type double in options");
}

bool Options::getBoolParam(const std::string& identifier) const
{
   const std::map<std::string, bool>::const_iterator it = bool_options.find(identifier);

   if( it != bool_options.end() )
      return it->second;
   else
      throw new std::string("No element \"" + identifier + "\" of type double in options");
}

void Options::setDoubleParam(const std::string& identifier, double value)
{
   double_options[identifier] = value;
}

void Options::setIntParam(const std::string& identifier, int value)
{
   int_options[identifier] = value;
}

void Options::setBoolParam(const std::string& identifier, bool value)
{
   bool_options[identifier] = value;
}

void Options::fillParamsFromFile(const std::string& filename)
{
   // TODO...
}
