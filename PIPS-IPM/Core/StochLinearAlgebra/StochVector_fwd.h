/*
* StochVector_fwd.h
*
*  Created on: 08.11.2019
*      Author: bzfkempk
*/
 
#ifndef STOCHVECTOR_FWD_H
#define STOCHVECTOR_FWD_H
 
template<typename T>
class StochVectorBase;
 
template<typename T>
class StochDummyVectorBase;
 
 
typedef StochVectorBase<double> StochVector;
typedef StochDummyVectorBase<double> StochDummyVector;

#endif