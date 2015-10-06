/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef PETSCINTERFACE_VECTOR_H
#define PETSCINTERFACE_VECTOR_H

#include "OoqpVector.h"
#include "petscksp.h"

class PetscInterfaceVector{
public:
  PetscInterfaceVector( OoqpVector *vec_in);
  virtual ~PetscInterfaceVector();
};

#endif





