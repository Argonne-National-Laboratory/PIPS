/*
 * PreprocessFactory.h
 *
 *  Created on: Dec 30, 2017
 *      Author: Daniel Rehfeldt
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_PREPROCESSFACTORY_H_
#define PIPS_IPM_CORE_QPPREPROCESS_PREPROCESSFACTORY_H_

#include "EquiStochScaler.h"
#include "StochPresolver.h"

enum ScalerType {SCALER_NONE, SCALER_EQUI_STOCH, SCALER_GEO_STOCH};
enum PresolverType {PRESOLVER_NONE, PRESOLVER_STOCH};

class PreprocessFactory
{
public:
      static Scaler* makeScaler(Data* data, ScalerType type)
      {
         switch( type )
            {
            case SCALER_EQUI_STOCH:
               return new EquiStochScaler(data);
            default:
               return 0;
            }
      };

      static Presolver* makePresolver(const Data* data, PresolverType type)
      {
         switch( type )
            {
            case PRESOLVER_STOCH:
               return new StochPresolver(data);
            default:
               return 0;
            }
      };

      static PreprocessFactory& getInstance()
      {
         static PreprocessFactory* instance = new PreprocessFactory;
         return *instance;
      }
private:
      // we need to make some given functions private to finish the definition of the singleton
      PreprocessFactory(){}

      PreprocessFactory(const PreprocessFactory &old); // disallow copy constructor
      const PreprocessFactory& operator=(const PreprocessFactory &old); //disallow assignment operator

      ~PreprocessFactory(){}
};



#endif /* PIPS_IPM_CORE_QPPREPROCESS_PREPROCESSFACTORY_H_ */
