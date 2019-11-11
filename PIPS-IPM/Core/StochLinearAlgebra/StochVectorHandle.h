#ifndef STOCHVECTORHANDLE
#define STOCHVECTORHANDLE

#include "IotrRefCount.h"
#include "SmartPointer.h"
#include "StochVector_fwd.h"

template<typename T>
using StochVectorBaseHandle = SmartPointer<StochVectorBase<T> >;

using StochVectorHandle = SmartPointer<StochVector>;

#endif
