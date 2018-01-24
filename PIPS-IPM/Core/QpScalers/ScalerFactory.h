/*
 * ScalerFactory.h
 *
 *  Created on: Dec 30, 2017
 *      Author: Daniel Rehfeldt
 */

#ifndef PIPS_IPM_CORE_QPSCALERS_SCALERFACTORY_H_
#define PIPS_IPM_CORE_QPSCALERS_SCALERFACTORY_H_

#include "EquiStochScaler.h"

enum ScalerType {SCALER_NONE, SCALER_EQUI_STOCH, SCALER_GEO_STOCH};

class ScalerFactory
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

      static ScalerFactory& getInstance()
      {
         static ScalerFactory* instance = new ScalerFactory;
         return *instance;
      }
private:
      // we need to make some given functions private to finish the definition of the singleton
      ScalerFactory(){}

      ScalerFactory(const ScalerFactory &old); // disallow copy constructor
      const ScalerFactory& operator=(const ScalerFactory &old); //disallow assignment operator

      ~ScalerFactory(){}
};



#endif /* PIPS_IPM_CORE_QPSCALERS_SCALERFACTORY_H_ */
