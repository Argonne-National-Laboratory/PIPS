/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef OoqpVectorHandle_H
#define OoqpVectorHandle_H

#include "IotrRefCount.h"
#include "SmartPointer.h"
#include "OoqpVector_fwd.h"
#include "pipsport.h"

#ifndef PRE_CPP11
	template<typename T>
	using OoqpVectorBaseHandle = SmartPointer<OoqpVectorBase<T> >;
#else
//	typedef SmartPointer<OoqpVectorBase<int> > OoqpVectorBaseHandle<int>;
//	typedef SmartPointer<OoqpVectorBase<int> > OoqpVectorBaseHandle<int>
//	typedef SmartPointer<OoqpVectorBase<double> > OoqpVectorBaseHandle<double>
#endif 

typedef SmartPointer<OoqpVector> OoqpVectorHandle;

#endif
