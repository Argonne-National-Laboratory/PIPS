#include "Options.h"
#include "pipsdef.h"

#include <fstream>
#include <sstream>

#ifdef PRE_CPP11
#include <cstdlib>
#endif

namespace base_options
{
   /// initialize static options storages
   std::map<std::string, double> Options::double_options;
   std::map<std::string, int> Options::int_options;
   std::map<std::string, bool> Options::bool_options;

   Options::Options()
   {
      /// INTERIOR-POINT ALGORITHM
      bool_options["IP_ACCURACY_REDUCED"] = false;
      bool_options["IP_PRINT_TIMESTAMP"] = false;
      bool_options["IP_STEPLENGTH_CONSERVATIVE"] = false;
   }

   int Options::getIntParam(const std::string& identifier) const
   {
      const std::map<std::string, int>::const_iterator& it = int_options.find(identifier);

      if( it != int_options.end() )
         return it->second;
      else
      {
         std::cout << "No element \"" << identifier << "\" of type int in options - using 0 instead" << std::endl;
         return 0;
      }
   }

   double Options::getDoubleParam(const std::string& identifier) const
   {
      const std::map<std::string, double>::const_iterator it = double_options.find(identifier);

      if( it != double_options.end() )
         return it->second;
      else
      {
         std::cout << "No element \"" << identifier << "\" of type double in options - using 0 instead" << std::endl;
         return 0;
      }
   }

   bool Options::getBoolParam(const std::string& identifier) const
   {
      const std::map<std::string, bool>::const_iterator it = bool_options.find(identifier);

      if( it != bool_options.end() )
         return it->second;
      else
      {
         std::cout << "No element \"" << identifier << "\" of type bool in options - using false instead" << std::endl;
         return false;
      }
   }

   void Options::setIntParam(const std::string& param, int value)
   {
      int_options[param] = value;
   }

   void Options::setBoolParam(const std::string& param, int value)
   {
      bool_options[param] = value;
   }

   void Options::setDoubleParam(const std::string& param, int value)
   {
      double_options[param] = value;
   }

   void Options::fillOptionsFromFile(const std::string& filename)
   {
      std::ifstream params;
      params.open(filename, std::ios::in);
      const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);

      if( !params.good() )
      {
         if( my_rank == 0 )
         {
            std::cout << "Failed to open provided options file \"" << filename << "\"" << std::endl;
            std::cout << "Using default configuration" << std::endl;
         }
         return;
      }

      std::string line;
      while( std::getline(params, line) )
      {
         /* skip empty line */
         if( line.compare("") == 0 )
            continue;

         std::istringstream iss(line);

         std::string identifier;
         std::string value;
         std::string type;
         if( !(iss >> identifier >> value >> type) )
         {
            /* some error while reading that line occured */
            if( my_rank == 0 )
               std::cout << "Error while reading line \"" << line << "\" from options file - skipping that line" << std::endl;
            continue;
         }

         if( identifier.length() == 0 || identifier[0] == '#' )
            continue;

         try
         {
            if( !identifierExists(identifier) )
            {
               if( my_rank == 0 )
                  std::cout << "Warning - unknown identifier - skipping it: " << identifier << std::endl;
            }
            else if( type.compare("int") == 0 || type.compare("integer") == 0 )
            {
   #ifdef PRE_CPP11
               int_options[identifier] = atoi(value.c_str());
   #else
               int_options[identifier] = std::stoi(value);
   #endif
            }
            else if( type.compare("double") == 0 )
            {
   #ifdef PRE_CPP11
               double_options[identifier] = atof(value.c_str());
   #else
               double_options[identifier] = std::stod(value);
   #endif
            }
            else if( type.compare("bool") == 0 || type.compare("boolean") == 0 )
            {
               if( value.compare("true") == 0 || value.compare("TRUE") == 0 || value.compare("True") == 0 )
                  bool_options[identifier] = true;
               else if( value.compare("false") == 0 || value.compare("FALSE") == 0 || value.compare("False") == 0 )
                  bool_options[identifier] = false;
               else
                  if( my_rank == 0 )
                     std::cout << "Unknown value \"" << value << "\" for bool parameter \"" << identifier << "\" in options file - skipping that line" << std::endl;
            }
            else
            {
               if( my_rank == 0 )
                  std::cout << "Unknown type \"" << type << "\" in options file (supported: bool/boolean, double, int/integer)- skipping that line" << std::endl;
            }
         }
         catch( std::exception& e )
         {
            if( my_rank == 0 )
                std::cout << "Error while reading line \"" << line << "\" from options file - could not parse value of " <<
                   identifier << " " << value << " - skipping that line" << std::endl;
         }
      }
   }

   bool Options::identifierExists( const std::string& identifier ) const
   {
      return int_options.find(identifier) != int_options.end()
            || bool_options.find(identifier) != bool_options.end()
            || double_options.find(identifier) != double_options.end();
   }

   bool Options::isIdentifierUnique( const std::string& identifier ) const
   {
      int n_found = 0;

      if( int_options.find(identifier) != int_options.end() )
         n_found++;
      if( bool_options.find(identifier) != bool_options.end() )
         n_found++;
      if( double_options.find(identifier) != double_options.end() )
         n_found++;

      return n_found <= 1;
   }
}
