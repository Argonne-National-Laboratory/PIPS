/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QpGenScaDriver.h"
#include "QpGenSca.h"
#include "MehrotraSolver.h"

int main( int argc, char *argv[] )
{
  MehrotraSolver  * solver = 0;
  QpGenSca * qpgen  = 0;

  // Passing nil arguments to keep old compilers that can't handle
  // explicit template instatiation happy.
  int result = qpgensca_solve( argc, argv, solver, qpgen );
  
  return result;

}
