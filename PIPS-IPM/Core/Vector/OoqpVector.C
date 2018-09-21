/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "OoqpVector.h"

OoqpVector::OoqpVector( int n_ )
{
  n = n_;
}

OoqpVector::~OoqpVector()
{
}

void
OoqpVector::writefToStreamStats( ostream& out, std::string prestring)
{
   double min;
   double max;
   int dummy;

   this->min(min, dummy);
   this->max(max, dummy);

   std::cout << prestring << " length=" << n << " min="<< min <<  "max=" << max << " infnorm=" << this->infnorm() << std::endl;
}
