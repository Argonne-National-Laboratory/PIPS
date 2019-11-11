/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef VECTORUTILITIES
#define VECTORUTILITIES

#include <iostream>
#include <fstream>
#include <cstring>

template<typename T>
void set_to_zero( T v[], int n, int stride )
{
  int i;
  int lenv = stride*n;
  for( i = 0; i < lenv; i += stride ) v[i] = 0;
}

template<typename T>
void writef_to_stream( T v[], int n, int stride,
		       std::ostream& out, const char format[] )
{
  int i;
  int lenv = n * stride;

  for( i = 0; i < lenv; i += stride ) {
    int j = 0;
    char c;
    while( (c = format[j]) != 0 ) {
      if( c != '%' ) {
	out << c;
      } else {
	// Brain-dead variable substitution, but good enough for this
	// simple case
	if( 0 == strncmp( "{value}", &format[j + 1], 7 ) ) {
	  out << v[i];
	  j += 7;
	} else if ( 0 == strncmp( "{index}", &format[j + 1], 7 ) ) {
	  out << i;
	  j += 7;
	} else {
	  out << c;
	}
      }
      j++;
    }
    out << std::endl;
  }
}

template<typename T>
void set_to_constant( T v[], int n, int stride, T c )
{
  int i;
  int lenv = n * stride;
  for( i = 0; i < lenv; i += stride ) v[i] = c;
}

template<typename T>
void add_constant( T v[], int n, int stride, T c )
{
  int i;
  int lenv = n * stride;
  for( i = 0; i < lenv; i += stride ) v[i] += c;
}

template<typename T>
T stepbound( T v[], int n, int incv,
				  T s[], int incs, T max )
{
  T bound = max;
  T *pv = v;
  T *ps = s, *lasts = s + n * incs;
  for( ; ps < lasts; ps += incs, pv += incv ) {
	T ss = *ps;
	if( ss < 0 ) {
	  T cbnd = - *pv / ss;
	  if( cbnd < bound ) bound = cbnd;
	}
  }
  return bound;
}

template<typename T>
T stepbound( T v[], int n, int incv,
				  T s[], int incs,
				  T b[], int incb, T u[], int incu,
				  T max )
{
  T bound = max;
  T *pv = v;
  T *ps = s, *lasts = s + n * incs;
  int i;
  for( i = 0; ps < lasts; ps += incs, pv += incv, i++ ) {
	T ss = *ps;
	T cbnd;
	if( ss > 0 ) {
	  cbnd = ( u[i*incu] - *pv ) / ss;
	} else if ( ss < 0 ) {
	  cbnd = ( b[i*incb] - *pv ) / ss;
	} else {
	  continue; // Next element
	}
	if( cbnd < bound ) bound = cbnd;
  }
  return bound;
}

template<typename T>
void axdzpy( int n, T alpha,
	     T x[], int incx, T z[], int incz,
	     T y[], int incy )
{
  T *px = x, *py = y, *pz = z, *lastx = x + incx * n;
  for( ; px < lastx; px += incx, py += incy, pz += incz ) {
    *py += alpha * (*px / *pz);
  }
}

template<typename T>
T find_blocking( T w[],     int n, int incw,
		    T wstep[],        int incwstep,
		    T u[],            int incu,
		    T ustep[],        int incustep,
		    T maxStep,
		    T *w_elt,          T *wstep_elt,
		    T *u_elt,          T *ustep_elt,
		    int& first_or_second )
{
  T bound = maxStep;

  int i = n - 1, lastBlocking = -1;

  // Search backward so that we find the blocking constraint of lowest
  // index. We do this to make things consistent with MPI's MPI_MINLOC,
  // which returns the processor with smallest rank where a min occurs.
  //
  // Still, going backward is ugly!
  T *pw     = w     + (n - 1) * incw;
  T *pwstep = wstep + (n - 1) * incwstep;
  T *pu     = u     + (n - 1) * incu;
  T *pustep = ustep + (n - 1) * incustep;

  while( i >= 0 ) {
    T temp = *pwstep;
    if( *pw > 0 && temp < 0 ) {
      temp = -*pw/temp;
      if( temp <= bound ) {
         bound = temp;
         lastBlocking = i;
         first_or_second = 1;
      }
    }
    temp = *pustep;
    if( *pu > 0 && temp < 0 ) {
      temp = -*pu/temp;
      if( temp <= bound ) {
         bound = temp;
         lastBlocking = i;
         first_or_second = 2;
      }
    }

    i--;
    if( i >= 0 ) {
      // It is safe to decrement the pointers
      pw     -= incw;
      pwstep -= incwstep;
      pu     -= incu;
      pustep -= incustep;
    }
  }

  if( lastBlocking > -1 ) {
    // fill out the elements
    *w_elt     = w[lastBlocking];
    *wstep_elt = wstep[lastBlocking];
    *u_elt     = u[lastBlocking];
    *ustep_elt = ustep[lastBlocking];
  } 
  return bound;
}

template<typename T>
void find_blocking_pd( const T w[], const int n,
          const T wstep[],
		    const T u[],
		    const T ustep[],
		    T& maxStep_primal, T& maxStep_dual,
		    T& w_elt,          T& wstep_elt,
		    T& u_elt,          T& ustep_elt,
		    T& w_elt_d,        T& wstep_elt_d,
		    T& u_elt_d,        T& ustep_elt_d,
			bool& primalBlocking, bool& dualBlocking )
{
   int lastBlockingPrimal = -1, lastBlockingDual = -1;

   for( int i = 0; i < n; i++ )
   {
      T temp = wstep[i];
      if( w[i] > 0 && temp < 0 )
      {
         temp = -w[i] / temp;
         if( temp < maxStep_primal )
         {
            maxStep_primal = temp;
            lastBlockingPrimal = i;
            primalBlocking = true;
         }
      }
      temp = ustep[i];
      if( u[i] > 0 && temp < 0 )
      {
         temp = -u[i] / temp;
         if( temp < maxStep_dual )
         {
            maxStep_dual = temp;
            lastBlockingDual = i;
            dualBlocking = true;
         }
      }
   }

  if( lastBlockingPrimal > -1 ) {
    // fill out the elements
    w_elt     = w[lastBlockingPrimal];
    wstep_elt = wstep[lastBlockingPrimal];
    u_elt     = u[lastBlockingPrimal];
    ustep_elt = ustep[lastBlockingPrimal];
  }
  if( lastBlockingDual > -1 ) {
     // fill out the elements
     w_elt_d     = w[lastBlockingDual];
     wstep_elt_d = wstep[lastBlockingDual];
     u_elt_d     = u[lastBlockingDual];
     ustep_elt_d = ustep[lastBlockingDual];
   }
}



#endif
