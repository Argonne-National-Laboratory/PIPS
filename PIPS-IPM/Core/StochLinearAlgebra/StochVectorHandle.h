#ifndef STOCHVECTORHANDLE
#define STOCHVECTORHANDLE

#include "IotrRefCount.h"
#include "SmartPointer.h"
#include "StochVector_fwd.h"
#include "pipsport.h"

#ifndef PRE_CPP11
	template<typename T>
	using StochVectorBaseHandle = SmartPointer<StochVectorBase<T> >;
#else
//	#define StochVectorBaseHandle<int> SmartPointer<StochVectorBase<int> >
//	#define StochVectorBaseHandle<double> SmartPointer<StochVectorBase<double> >
#endif 

typedef SmartPointer<StochVector> StochVectorHandle;

#endif
