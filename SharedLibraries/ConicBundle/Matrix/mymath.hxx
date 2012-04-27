/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/mymath.hxx

    This file is part of ConciBundle, a C/C++ library for convex optimization.

    ConicBundle is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ConicBundle is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

***************************************************************************** */



#ifndef CH_MATRIX_CLASSES_MYMATH_HXX
#define CH_MATRIX_CLASSES_MYMATH_HXX

/**  @file mymath.hxx
    @brief Header defining simple functions like max, min, abs.

    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg

*/

#include <cmath>
#include <cstdlib>

namespace CH_Matrix_Classes {

  /**@defgroup mymath simple functions like max, min, abs 
   */
  //@{
  
/// absolute value of a double
inline double abs(double d)
{
  return std::fabs(d);
}
  
/// absolute value of an int
inline int abs(int d)
{
  return std::abs(d);
}

/// absolute value of a long
inline long abs(long d)
{
  return std::labs(d);
}
  
/// maximum value of two double variables
inline double max(double a,double b)
{
 return (a>=b)?a:b;
}

/// minimum value of two double variables
inline double min(double a,double b)
{
 return (a<=b)?a:b;
}

/// maximum value of two int variables
inline int max(int a,int b)
{
 return (a>=b)?a:b;
}

/// minimum value of two int variables
inline int min(int a,int b)
{
 return (a<=b)?a:b;
}

/// maximum value of two long vriables
inline long max(long a,long b)
{
 return (a>=b)?a:b;
}

/// minimum value of two long variables 
inline long min(long a,long b)
{
 return (a<=b)?a:b;
}

/// swaps two double variables
inline void swap(double& a,double& b)
{
 double h=a; a=b; b=h;
}

/// swaps two long variables
inline void swap(long& a,long& b)
{
 long h=a; a=b; b=h;
}

/// swpas two int variables
inline void swap(int& a,int& b)
{
 int h=a; a=b; b=h;
}

/// swpas two bool variables
inline void swap(bool& a,bool& b)
{
 bool h=a; a=b; b=h;
}

/// return a*a for int a
inline int sqr(int a)
{
 return a*a;
}

/// return a*a for long a
inline long sqr(long a)
{
 return a*a;
}

/// return a*a for double a
inline double sqr(double a)
{
 return a*a;
}

/// return sqrt for int a
inline double sqrt(int a)
{
  return std::sqrt(double(a));
}

/// return sqrt for long a
inline double sqrt(long a)
{
  return std::sqrt(double(a));
}

/// return the signum of an int a (1 for a>0,-1 for a<0,0 for a==0)
inline int sign(int a)
{
 if (a>0) return 1;
 if (a<0) return -1;
 return 0;
}

/// return the signum of a long a (1 for a>0,-1 for a<0,0 for a==0)
inline long sign(long a)
{
 if (a>0) return 1;
 if (a<0) return -1;
 return 0;
}

/// return the signum of a double a (1. for a>0.,-1. for a<0.,0. for a==0.)
inline double sign(double a)
{
 return (a>=0)?1.:-1.;
}

/// return the signum of a double a with tolerance (1. for a>tol,-1. for a<-tol,0. otherwise)
inline double sign(double a,double tol)
{
 if (a>tol) return 1.;
 if (a<-tol) return -1.;
 return 0.;
}

//---- functions for f2c translations

/// return a if a and b have the same sign, return -a otherwise
inline double d_sign(double a,double b)
{
 return (((b>=0)&&(a>=0))||((b<0)&&(a<0)))?a:-a;
}

//@}

}

#endif

