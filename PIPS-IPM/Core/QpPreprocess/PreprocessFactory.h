/*
 * PreprocessFactory.h
 *
 *  Created on: Dec 30, 2017
 *      Author: Daniel Rehfeldt
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_PREPROCESSFACTORY_H_
#define PIPS_IPM_CORE_QPPREPROCESS_PREPROCESSFACTORY_H_

#include "EquiStochScaler.h"
#include "GeoStochScaler.h"
#include "StochPresolver.h"
#include "IotrRefCount.h"
#include "SmartPointer.h"
#include "StochPostsolver.h"
#include "sData.h"

enum ScalerType {SCALER_NONE, SCALER_EQUI_STOCH, SCALER_GEO_STOCH, SCALER_GEO_EQUI_STOCH};
enum PresolverType {PRESOLVER_NONE, PRESOLVER_STOCH};

class PreprocessFactory : public IotrRefCount
{
public:
  
      static Scaler* makeScaler(Data* data, ScalerType type)
      {
         switch( type )
            {
            case SCALER_EQUI_STOCH:
               return new EquiStochScaler(data, false);
            case SCALER_GEO_STOCH:
               return new GeoStochScaler(data, false, false);
            case SCALER_GEO_EQUI_STOCH:
               return new GeoStochScaler(data, true, false);
            default:
               return 0;
            }
      };

      static Presolver* makePresolver(const Data* data, PresolverType type, Postsolver* postsolver = NULL)
      {
         switch( type )
            {
            case PRESOLVER_STOCH:
               return new StochPresolver(data, postsolver);
            default:
               return 0;
            }
      };

      // todo : not sure about the factory design here
      static Postsolver* makePostsolver(const Data* original_problem)
      {
         return new StochPostsolver(dynamic_cast<const sData&>(*original_problem));
      }
};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_PREPROCESSFACTORY_H_ */
