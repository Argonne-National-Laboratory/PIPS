/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SIMPLEVECTORHANDLE
#define SIMPLEVECTORHANDLE

#include "OoqpVectorHandle.h"
#include "SmartPointer.h"
#include "SimpleVector_fwd.h"

// template<typename T>
// using SimpleVectorBaseHandle = SmartPointer<SimpleVectorBase<T> >;

typedef SmartPointer<SimpleVector> SimpleVectorHandle;

#endif
