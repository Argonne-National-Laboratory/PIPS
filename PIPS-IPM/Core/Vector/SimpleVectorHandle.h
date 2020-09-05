/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SIMPLEVECTORHANDLE
#define SIMPLEVECTORHANDLE

#include "OoqpVectorHandle.h"
#include "SimpleVector_fwd.h"

#include "pipsport.h"

#ifndef PRE_CPP11
	template<typename T>
	using SimpleVectorBaseHandle = SmartPointer<SimpleVectorBase<T> >;
#else
//	#define SimpleVectorBaseHandle<int> SmartPointer<StochVectorBase<int> >
//	#define SimpleVectorBaseHandle<double> SmartPointer<StochVectorBase<double> >
#endif 

typedef SmartPointer<SimpleVector> SimpleVectorHandle;

#endif
