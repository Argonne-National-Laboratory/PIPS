/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

template<typename T>
OoqpVectorBase<T>::OoqpVectorBase( int n_ )
{
  n = n_;
}

template<typename T>
OoqpVectorBase<T>::~OoqpVectorBase()
{
}

template<typename T>
void OoqpVectorBase<T>::writefToStreamStats( std::ostream& out, std::string prestring)
{
   T min;
   T max;
   int dummy;

   this->min(min, dummy);
   this->max(max, dummy);

   std::cout << prestring << " length=" << n << " min=" << min <<  " max=" <<
      max << " infnorm=" << this->infnorm() << std::endl;
}
