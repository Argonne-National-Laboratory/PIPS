/*
 * StochVectorUtilities.h
 *
 *  Created on: 13.11.2019
 *      Author: bzfkempk
 */
#ifndef PIPS_IPM_CORE_STOCHLINEARALGEBRA_STOCHVECTORUTILITIES_H_
#define PIPS_IPM_CORE_STOCHLINEARALGEBRA_STOCHVECTORUTILITIES_H_

#include "OoqpVector.h"
#include "StochVector.h"


/* Utility functions to return for a given StochVectorBase<T> a certain SimpleVectorBase<T> 
 *	either from one of the children (node = 0, ..., nChildren - 1 )
 * or from the parent node = -1
 *
 * for the col vec function SimpleVectorBase<T> vec is returned 
 *
 * for row vec function either the associated SimpleVectorBase<T> vec (linking = false) or vecl (linking = true) is returned
 *
 * asserts existence of these vectors as well as the specified child
 */
template <typename T>
SimpleVectorBase<T>& getSimpleVecFromStochVec(const StochVectorBase<T>& stochvec, int node, bool linking)
{
   assert(-1 <= node && node < static_cast<int>(stochvec.children.size()) );

   if(node == -1)
   {
      if( linking )
      {
         assert(stochvec.vecl);
         return dynamic_cast<SimpleVectorBase<T>&>(*(stochvec.vecl));
      }
      else
      {
         assert(stochvec.vec);
         return dynamic_cast<SimpleVectorBase<T>&>(*(stochvec.vec));
      }
   }
   else
   {
      if( linking )
      {
         assert(stochvec.vecl);
         return dynamic_cast<SimpleVectorBase<T>&>(*(stochvec.vecl));
      }
      else
      {
         assert(stochvec.children[node]->vec);
         return dynamic_cast<SimpleVectorBase<T>&>(*(stochvec.children[node]->vec));
      }
   }
}

template <typename T>
SimpleVectorBase<T>& getSimpleVecFromRowStochVec(const OoqpVectorBase<T>& ooqpvec, int node, bool linking)
{ return getSimpleVecFromStochVec(dynamic_cast<const StochVectorBase<T>&>(ooqpvec), node, linking); };

template <typename T>
SimpleVectorBase<T>& getSimpleVecFromColStochVec(const OoqpVectorBase<T>& ooqpvec, int node)
   { return getSimpleVecFromStochVec(dynamic_cast<const StochVectorBase<T>&>(ooqpvec), node, false); };



/// clone the structure of a StochVectorBase<T> into one of type U
template<typename T, typename U>
StochVectorBase<U>* cloneStochVector(const StochVectorBase<T>& svec)
{
  StochVectorBase<U>* clone;
  
  if( svec.isKindOf(kStochDummy) )
    return new StochDummyVectorBase<U>();

  if( svec.vecl )
    clone = new StochVectorBase<U>( svec.vec->length(), svec.vecl->length(), svec.mpiComm, -1);
  else
    clone = new StochVectorBase<U>( svec.vec->length(), svec.mpiComm);

  for(size_t it = 0; it < svec.children.size(); it++) {
    clone->AddChild( cloneStochVector<T,U>(*svec.children[it]) );
  }
  return clone;
}

template<typename T, typename U>
StochDummyVectorBase<U>* cloneStochVector(const StochDummyVectorBase<T>& dummyvec)
{
  return new StochDummyVectorBase<U>();
}

template<typename T, typename U>
StochVectorBase<U>* cloneStochVector(const OoqpVectorBase<T>& ooqpvec)
{
  if( ooqpvec.isKindOf(kStochDummy) )
    return cloneStochVector<T, U>( dynamic_cast<const StochDummyVectorBase<T>&>(ooqpvec) );
  else
    return cloneStochVector<T, U>( dynamic_cast<const StochVectorBase<T>&>(ooqpvec) );
}

#endif