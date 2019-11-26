/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "SimpleVectorHandle.h"
#include "VectorUtilities.h"
#include "SimpleVector.h"
#include "OoqpBlas.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <algorithm>

template <typename T>
long long SimpleVectorBase<T>::numberOfNonzeros() const
{
  long long i, count = 0;
  for( i = 0; i < this->n; i++ ) {
    if( v[i] != 0 ) count++;
  }
  return count;
}

template<typename T>
void SimpleVectorBase<T>::min( T& m, int& index ) const
{
  if( this->n == 0 ) {
    m = std::numeric_limits<T>::max();
    return;
  }
  index = 0;
  m     = v[0];
  for( int i = 0; i < this->n; i++ ) {
    if( v[i] < m ) {
      m     = v[i];
      index = i;
    }
  }
}

template<typename T>
void SimpleVectorBase<T>::absminVecUpdate( OoqpVectorBase<T>& absminvec) const
{
   const SimpleVectorBase<T>& absminvecSimple = dynamic_cast<const SimpleVectorBase<T>&>(absminvec);
   assert( absminvecSimple.length() == this->n );
   T* const absminvecArr = absminvecSimple.elements();

   for( int i = 0; i < this->n; i++ )
   {
      const T abs = fabs(v[i]);
      if( abs < absminvecArr[i] && abs > pips_eps )
         absminvecArr[i] = abs;
   }
}

template<typename T>
void SimpleVectorBase<T>::absmaxVecUpdate( OoqpVectorBase<T>& absmaxvec) const
{
   const SimpleVectorBase<T>& absmaxvecSimple = dynamic_cast<const SimpleVectorBase<T>&>(absmaxvec);
   assert( absmaxvecSimple.length() == this->n );
   T* const absmaxvecArr = absmaxvecSimple.elements();

   for( int i = 0; i < this->n; i++ )
   {
      const T abs = fabs(v[i]);
      if( abs > absmaxvecArr[i] )
         absmaxvecArr[i] = abs;
   }
}

template<typename T>
void SimpleVectorBase<T>::absmin(T& min) const
{
   if (this->n == 0) {
     min = std::numeric_limits<T>::max();
     return;
   }
   min = fabs(v[0]);
   for( int i = 0; i < this->n; i++ ) {
     if( fabs(v[i]) < min ) {
        min = fabs(v[i]);
     }
   }
}

/** Compute the min absolute value that is larger than zero_eps.
 * If there is no such value, return -1.0 */
 template<typename T>
void SimpleVectorBase<T>::absminNonZero(T& m, T zero_eps) const
{
   assert(zero_eps >= 0.0);

   m = -1.0;

   if( this->n == 0 )
      return;

   T min = std::numeric_limits<T>::max();

   for( int i = 0; i < this->n; i++ )
   {
      if( fabs(v[i]) < min && fabs(v[i]) > zero_eps )
         min = fabs(v[i]);
   }

   if( min < std::numeric_limits<T>::max() )
      m = min;
}

template<typename T>
void SimpleVectorBase<T>::max( T& m, int& index ) const
{
   if( this->n == 0 )
   {
      index = -1;
      m = -std::numeric_limits<T>::max();
      return;
   }
   index = 0;
   m = v[0];
   for( int i = 0; i < this->n; i++ )
   {
      if( v[i] > m )
      {
         m = v[i];
         index = i;
      }
   }
}

template<typename T>
bool SimpleVectorBase<T>::isKindOf( int kind ) const
{
  return (kind == kSimpleVector);
}

template<typename T>
void SimpleVectorBase<T>::copyIntoArray( T w[] ) const
{
  memcpy( w, this->v, this->n * sizeof( T ) );
}

template<typename T>
void SimpleVectorBase<T>::copyFromArray( const T w[] )
{
  memcpy( this->v, w, this->n * sizeof( T ) );
}

template<typename T>
void SimpleVectorBase<T>::copyFromArray( const char w[] )
{
  int i;
  for( i = 0; i < this->n; i++ ) {
    this->v[i] = w[i];
  }
}

template<typename T>
SimpleVectorBase<T>::SimpleVectorBase( int n_ ) : OoqpVectorBase<T>( n_ )
{
  assert(this->n >= 0);
  preserveVec = 0;
  v = new T[this->n];
  memset(v, 0, this->n * sizeof(T));
}

template<typename T>
SimpleVectorBase<T>::SimpleVectorBase( T * v_, int n_ )
  : OoqpVectorBase<T>( n_ )
{
  preserveVec = 1;
  v = v_;
}

template<typename T>
SimpleVectorBase<T>::~SimpleVectorBase()
{
  if( !preserveVec ) {
    delete [] v;
  }
}

template<typename T>
OoqpVectorBase<T>* SimpleVectorBase<T>::cloneFull() const
{
   SimpleVectorBase<T>* clone = new SimpleVectorBase<T>(this->n);
   clone->copyFromArray(v);

   return clone;
}

template<typename T>
bool SimpleVectorBase<T>::isZero() const
{
	bool is_zero = true;

	for(int i = 0; i < this->n; ++i)
		is_zero = (is_zero && v[i] == 0.0);

	return is_zero;
}

template<typename T>
void SimpleVectorBase<T>::setToZero()
{
  int i;
  for( i = 0; i < this->n; i++ ) v[i] = 0.0;
}

template<typename T>
void SimpleVectorBase<T>::setToConstant( T c)
{
  int i;
  for( i = 0; i < this->n; i++ ) v[i] = c;
}

// specialiced for double only
template<>
void SimpleVectorBase<double>::randomize( double alpha, double beta, double *ix )
{
  assert( beta > alpha);

  double drand(double *);
  double scale = beta - alpha;
  double shift = alpha/scale;

  int i;
  for( i = 0; i < this->n; i++ ) {
    v[i] = scale * (drand(ix) + shift);
  }
}

template<typename T>
void SimpleVectorBase<T>::randomize( T alpha, T beta, T *ix )
{
  assert( 0 && "not implemented here" );
  return;
}

template<typename T>
void SimpleVectorBase<T>::copyFrom( const OoqpVectorBase<T>& vec )
{
  assert( vec.length() == this->n );

  vec.copyIntoArray( this->v );
}

template<typename T>
void SimpleVectorBase<T>::copyFromAbs(const OoqpVectorBase<T>& vec )
{
   const SimpleVectorBase<T>& vecSimple = dynamic_cast<const SimpleVectorBase<T>&>(vec);
   assert( vec.length() == this->n );
   T* const vecArr = vecSimple.elements();

   for( int i = 0; i < this->n; i++ )
     v[i] = fabs( vecArr[i] );
}

template<typename T>
T SimpleVectorBase<T>::infnorm() const
{
  T temp, norm = 0;
  int i;
  for( i = 0; i < this->n; i++ ) {
    temp = fabs( v[i] );
    // Subtle reversal of the logic to handle NaNs
    if( ! ( temp <=  norm ) ) norm = temp;
  }

  return norm;
}

template<typename T>
T SimpleVectorBase<T>::onenorm() const
{
  T temp, norm = 0;
  int i;
  for( i = 0; i < this->n; i++ ) {
    temp = fabs( v[i] );
    norm += temp;
  }
  return norm;
}

template<typename T>
double SimpleVectorBase<T>::twonorm() const
{
  T temp = dotProductWith(*this);
  return sqrt(temp);
}

template<typename T>
void SimpleVectorBase<T>::componentMult( const OoqpVectorBase<T>& vec )
{
  assert( this->n == vec.length() );
  const SimpleVectorBase<T> & sv = dynamic_cast<const SimpleVectorBase<T> &>(vec);
  T * y = sv.v;
  int i;
  for( i = 0; i < this->n; i++ ) v[i] *= y[i];
}

template<typename T>
void SimpleVectorBase<T>::scalarMult( T num)
{
  int i;
  for( i = 0; i < this->n; i++) v[i] *= num;
}

// Print first 10 entries of solution vector to stderr.
// Useful for debugging purposes...
template<typename T>
void SimpleVectorBase<T>::printSolutionToStdErr( OoqpVectorBase<T> &vec)
{
  int i;
  for( i = 0; i < 10; i++ )
  {
     std::cerr << v[i] << std::endl;
   }
  std::cerr << "******" << std::endl;
}

template<typename T>
void SimpleVectorBase<T>::componentDiv ( const OoqpVectorBase<T>& vec )
{
  assert( this->n == vec.length() );
  T * pv = v, *lv = v + this->n;

  const SimpleVectorBase<T> & sv = dynamic_cast<const SimpleVectorBase<T> &>(vec);
  T * y = sv.v;

  for( ; pv < lv; pv++, y++ ) *pv /= *y;
}

template<typename T>
void SimpleVectorBase<T>::writeToStream(std::ostream& out) const
{
  this->writefToStream( out, "%{value}" );
}

template<typename T>
void SimpleVectorBase<T>::writeToStreamAll(std::ostream& out) const
{
   for( int i = 0; i < this->n; i++ )
      out << v[i] << "\n";
}

template<typename T>
void SimpleVectorBase<T>::writeToStreamAllStringStream(std::stringstream& sout) const
{
   for( int i = 0; i < this->n; i++ )
      sout << v[i] << "\n";
}

template<typename T>
void SimpleVectorBase<T>::writefToStream( std::ostream& out, const char format[] ) const
{
  SmartPointer<SimpleVectorBase<T>> empty( new SimpleVectorBase<T>(0) );
  this->writefSomeToStream( out, format, *empty );
}

template<typename T>
void SimpleVectorBase<T>::writefSomeToStream( std::ostream& out,
					   const char format[],
					   const OoqpVectorBase<T>& select ) const
{
  const SimpleVectorBase<T> & sselect = dynamic_cast<const SimpleVectorBase<T> &>(select);
  T * s = 0;
  if( select.length() > 0 ) {
    s = sselect.v;
  }
  int i;

  for( i = 0; i < this->n; i++ ) {
    if( !s || s[i] != 0.0 ) {
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
}

template<typename T>
void SimpleVectorBase<T>::writeMPSformatOnlyRhs(std::ostream& out, const std::string rowName, const OoqpVectorBase<T>* irhs) const
{
   if( irhs )
      assert( this->n == irhs->n );

   for( int i = 0; i < this->n; i++ )
   {
      if( !irhs || ( irhs && dynamic_cast<const SimpleVectorBase<T> *>(irhs)->elements()[i] != 0.0 ) )
         out << rowName << i << " " << v[i] << "\n";
   }
}

template<typename T>
void SimpleVectorBase<T>::writeMPSformatBoundsWithVar(std::ostream& out, const std::string varStub, const OoqpVectorBase<T>* ix, bool upperBound) const
{
   assert( this->n == dynamic_cast<const SimpleVectorBase<T>*>(ix)->n );
   std::string boundType = (upperBound) ? " UP" : " LO";
   std::string infiniteBound = (upperBound) ? " PL" : " MI";

   for( int i = 0; i < this->n; i++ )
   {
      if( dynamic_cast<const SimpleVectorBase<T>*>(ix)->elements()[i] != 0.0  )
         out << boundType << " BND " << varStub << i << " " << v[i] << "\n";
      else
         out << infiniteBound << " BND " << varStub << i << " " << "\n";
   }
}

template<>
void SimpleVectorBase<double>::scale( double alpha )
{
  int one = 1;
  dscal_( &this->n, &alpha, v, &one );
}

// generic implementation without boost 
template<typename T>
void SimpleVectorBase<T>::scale( T alpha )
{
  assert(0 && "not implemented here");
   // std::transform( this->v, this->v + this->n, this->v, [alpha](T a)->T { return alpha * a; } );
}

template<>
void SimpleVectorBase<double>::axpy( double alpha, const OoqpVectorBase<double>& vec )
{
  assert( this->n == vec.length() );
  const SimpleVectorBase<double> & sv = dynamic_cast<const SimpleVectorBase<double> &>(vec);

  int one = 1;
  daxpy_( &this->n, &alpha, sv.v, &one, v, &one );
}

template<typename T>
void SimpleVectorBase<T>::axpy( T alpha, const OoqpVectorBase<T>& vec )
{
  assert(0 && "not implemented here");

  // assert( this->n == vec.length() );
  // const SimpleVectorBase<T> & sv = dynamic_cast<const SimpleVectorBase<T> &>(vec);
  // std::transform( this->v, this->v + this->n, sv.v, this->v,
      // [alpha](T a, T b)->T { return a + alpha * b; });
}

template<typename T>
void SimpleVectorBase<T>::addConstant( T c )
{
  int i;
  for( i = 0; i < this->n; i++ ) v[i] += c;
}

template<typename T>
void SimpleVectorBase<T>::gondzioProjection( T rmin, T rmax )
{
  int i;
  for( i = 0; i < this->n; i++ ) {
    if( v[i] < rmin ) {
      v[i] = rmin - v[i];
    } else if ( v[i] > rmax ) {
      v[i] = rmax - v[i];
    } else {
      v[i] = 0.0;
    }

    if( v[i] < -rmax ) v[i] = -rmax;
  }
}

template<typename T>
void SimpleVectorBase<T>::axzpy( T alpha, const OoqpVectorBase<T>& xvec,
			      const OoqpVectorBase<T>& zvec )
{
  assert( this->n == xvec.length() &&
	  this->n == zvec.length() );

  const SimpleVectorBase<T> & sxvec = dynamic_cast<const SimpleVectorBase<T> &>(xvec);
  const SimpleVectorBase<T> & szvec = dynamic_cast<const SimpleVectorBase<T> &>(zvec);

  T * x = sxvec.v;
  T * z = szvec.v;
  T * lx = x + this->n;
  T * w = v;

  if( alpha == 1.0 ) {
    while( x < lx ) {
      *w += *x * *z;
      w++; x++; z++;
    }
  } else if ( alpha == -1 ) {
    while( x < lx ) {
      *w -= *x * *z;
      w++; x++; z++;
    }
  } else {
    while( x < lx ) {
      *w += alpha * *x * *z;
      w++; x++; z++;
    }
  }
}

template<typename T>
void SimpleVectorBase<T>::axdzpy( T alpha, const OoqpVectorBase<T>& xvec,
   const OoqpVectorBase<T>& zvec )
{
  const SimpleVectorBase<T> & sxvec = dynamic_cast<const SimpleVectorBase<T> &>(xvec);
  T * x = sxvec.v;
  const SimpleVectorBase<T> & szvec = dynamic_cast<const SimpleVectorBase<T> &>(zvec);
  T * z = szvec.v;

  assert( this->n == xvec.length() &&
	  this->n == zvec.length() );

  int i;
  for( i = 0; i < this->n; i++ ) {
    //if(x[i] > 0 && z[i] > 0)
      v[i] += alpha * x[i] / z[i];
  }
}

template<typename T>
void SimpleVectorBase<T>::axdzpy( T alpha, const OoqpVectorBase<T>& xvec,
			       const OoqpVectorBase<T>& zvec, const OoqpVectorBase<T>& select )
{
  assert( this->n == xvec.length() &&
	  this->n == zvec.length() );

  const SimpleVectorBase<T> & sxvec = dynamic_cast<const SimpleVectorBase<T> &>(xvec);
  T * x = sxvec.v;
  const SimpleVectorBase<T> & szvec = dynamic_cast<const SimpleVectorBase<T> &>(zvec);
  T * z = szvec.v;
  const SimpleVectorBase<T> & sselect = dynamic_cast<const SimpleVectorBase<T> &>(select);
  T * s   = sselect.v;
  int i;
  if( alpha == 1.0 ) {
    for( i = 0; i < this->n; i++ ) {
      if( 0.0 != s[i] ) v[i] += x[i] / z[i];
    }
  } else if ( alpha == -1.0 ) {
    for( i = 0; i < this->n; i++ ) {
      if( 0.0 != s[i] ) v[i] -= x[i] / z[i];
    }
  } else {
    for( i = 0; i < this->n; i++ ) {
      if( 0.0 != s[i] ) v[i] += alpha * x[i] / z[i];
    }
  }
}

template<typename T>
T SimpleVectorBase<T>::dotProductWith( const OoqpVectorBase<T>& vec ) const
{
  assert( this->n == vec.length() );
  const SimpleVectorBase<T> & svec = dynamic_cast<const SimpleVectorBase<T> &>(vec);
  T * vvec = svec.v;

  T dot1 = 0.0;
  T dot2 = 0.0;

  const int size = 8196;
  int kmax       = this->n / size;
  int i          = 0;

  int k, imax;
  for( k = 0; k < kmax; k++ ) {
    imax = (k + 1) * size;
    for( ; i < imax; i++ ) {
      dot1 += v[i] * vvec[i];
    }
    dot2 += dot1;
    dot1  = 0;
  }
  for( ; i < this->n; i++ ) {
    dot1 += v[i] * vvec[i];
  }

  return dot2 + dot1;
}

template<typename T>
T SimpleVectorBase<T>::dotProductSelf( T scaleFactor ) const
{
   assert(scaleFactor >= 0.0);

   T dot = 0.0;

   if( scaleFactor == 1.0 )
   {
      for( int i = 0; i < this->n; i++ )
         if( !PIPSisZero(v[i]) )
            dot += v[i] * v[i];
   }
   else
   {
      for( int i = 0; i < this->n; i++ )
      {
         const T valScaled = v[i] * scaleFactor;
         if( !PIPSisZero(valScaled) )
            dot += valScaled * valScaled;
      }
   }
   return dot;
}

template<typename T>
T SimpleVectorBase<T>::shiftedDotProductWith( T alpha, const OoqpVectorBase<T>& mystep,
					 const OoqpVectorBase<T>& yvec,
					 T beta, const OoqpVectorBase<T>& ystep ) const
{
  assert( this->n == mystep.length() &&
	  this->n == yvec.length() &&
	  this->n == ystep.length() );

  const SimpleVectorBase<T> & syvec = dynamic_cast<const SimpleVectorBase<T> &>(yvec);
  T * y = syvec.v;

  const SimpleVectorBase<T> & smystep = dynamic_cast<const SimpleVectorBase<T> &>(mystep);
  T * p = smystep.v;

  const SimpleVectorBase<T> & systep = dynamic_cast<const SimpleVectorBase<T> &>(ystep);
  T * q = systep.v;

  T dot1 = 0.0;
  T dot2 = 0.0;

  const int size = 8196;
  int kmax       = this->n / size;
  int i          = 0;

  int k, imax;
  for( k = 0; k < kmax; k++ ) {
    imax = (k + 1) * 8196;
    for( ; i < imax; i++ ) {
      dot1 += (v[i] + alpha * p[i]) * (y[i] + beta * q[i] );
    }
    dot2 += dot1;
    dot1  = 0;
  }
  for( ; i < this->n; i++ ) {
    dot1 += (v[i] + alpha * p[i]) * (y[i] + beta * q[i] );
  }

  return dot2 + dot1;
}

template<typename T>
void SimpleVectorBase<T>::negate()
{
  int i;
  for( i = 0; i < this->n; i++ ) v[i] = -v[i];
}

template<typename T>
void SimpleVectorBase<T>::invert()
{
  for( int i = 0; i < this->n; i++ )
  {
    assert(v[i] != 0.0);
    v[i] = 1 / v[i];
  }
}

template<typename T>
void SimpleVectorBase<T>::invertSave( T zeroReplacementVal )
{
  for( int i = 0; i < this->n; i++ )
  {
     if( v[i] != 0.0 )
        v[i] = 1 / v[i];
     else
        v[i] = zeroReplacementVal;
  }
}

template<typename T>
void SimpleVectorBase<T>::applySqrt()
{
   for( int i = 0; i < this->n; i++)
   {
      assert( v[i] >= 0.0 );
      v[i] = std::sqrt(v[i]);
   }
}

template<typename T>
void SimpleVectorBase<T>::roundToPow2()
{
  for( int i = 0; i < this->n; i++ )
  {
     int exp;
#if 0
     const double mantissa = std::frexp(v[i], &exp);

     if( mantissa >= 0.75 )
        v[i] = std::ldexp(0.5, exp + 1);
     else
        v[i] = std::ldexp(0.5, exp);
#else
     (void) std::frexp(v[i], &exp);
     v[i] = std::ldexp(0.5, exp);
#endif
  }
}

template<typename T>
bool SimpleVectorBase<T>::allPositive() const
{
  int i;
  for( i = 0; i < this->n; i++ ) {
    if( v[i] <= 0 ) return false;
  }
  return true;
}

template<typename T>
T SimpleVectorBase<T>::stepbound( const OoqpVectorBase<T> & pvec, T maxStep ) const
{
  assert( this->n == pvec.length() );

  const SimpleVectorBase<T> & spvec = dynamic_cast<const SimpleVectorBase<T> &>(pvec);
  T * p = spvec.v;
  T * w = v;
  T bound = maxStep;

  int i;
  for( i = 0; i < this->n; i++ ) {
    T temp = p[i];
    if( w[i] >= 0 && temp < 0 ) {
      temp = -w[i]/temp;
      if( temp < bound ) {
	        bound = temp;
      }
    }
  }
  return bound;
}

template<typename T>
T SimpleVectorBase<T>::findBlocking( const OoqpVectorBase<T> & wstep_vec,
				      const OoqpVectorBase<T> & u_vec,
				      const OoqpVectorBase<T> & ustep_vec,
                  T maxStep,
				      T *w_elt,
				      T *wstep_elt,
				      T *u_elt,
				      T *ustep_elt,
				      int& first_or_second) const
{
  T * w     = v;
  const SimpleVectorBase<T> & swstep = dynamic_cast<const SimpleVectorBase<T> &>(wstep_vec);
  T * wstep = swstep.v;

  const SimpleVectorBase<T> & su_vec = dynamic_cast<const SimpleVectorBase<T> &>(u_vec);
  T * u     = su_vec.v;

  const SimpleVectorBase<T> & sustep_vec = dynamic_cast<const SimpleVectorBase<T> &>(ustep_vec);
  T * ustep = sustep_vec.v;

  return ::find_blocking( w, this->n, 1, wstep, 1, u, 1, ustep, 1, maxStep,
			  w_elt, wstep_elt, u_elt, ustep_elt,
			  first_or_second );
}

template<typename T>
void SimpleVectorBase<T>::findBlocking_pd(const OoqpVectorBase<T> & wstep_vec,
						const OoqpVectorBase<T> & u_vec, const OoqpVectorBase<T> & ustep_vec,
						T& maxStepPri, T& maxStepDual,
						T& w_elt_p, T& wstep_elt_p, T& u_elt_p, T& ustep_elt_p,
						T& w_elt_d, T& wstep_elt_d, T& u_elt_d, T& ustep_elt_d,
						bool& primalBlocking, bool& dualBlocking) const {
	const T * w = v;
	const SimpleVectorBase<T> & swstep = dynamic_cast<const SimpleVectorBase<T> &>(wstep_vec);
	const T * wstep = swstep.v;

	const SimpleVectorBase<T> & su_vec = dynamic_cast<const SimpleVectorBase<T> &>(u_vec);
	const T * u = su_vec.v;

	const SimpleVectorBase<T> & sustep_vec = dynamic_cast<const SimpleVectorBase<T> &>(ustep_vec);
	const T * ustep = sustep_vec.v;

	::find_blocking_pd(w, this->n, wstep, u, ustep, maxStepPri,
			maxStepDual, w_elt_p, wstep_elt_p, u_elt_p, ustep_elt_p, w_elt_d,
			wstep_elt_d, u_elt_d, ustep_elt_d,
			primalBlocking, dualBlocking);
}

template<typename T>
bool SimpleVectorBase<T>::matchesNonZeroPattern( const OoqpVectorBase<T>& select ) const
{
  const SimpleVectorBase<T> & sselect = dynamic_cast<const SimpleVectorBase<T> &>(select);
  T * map = sselect.v;

  T * lmap = map + this->n;
  assert( this->n == select.length() );

  T *w = v;
  while( map < lmap ) {
    if( *map == 0.0 && *w != 0.0  ) return false;
    map++;
    w++;
  }

  return true;
}

template<typename T>
void SimpleVectorBase<T>::selectNonZeros( const OoqpVectorBase<T>& select )
{
  const SimpleVectorBase<T> & sselect = dynamic_cast<const SimpleVectorBase<T> &>(select);
  T * map = sselect.v;

  assert( this->n == select.length() );
  int i;
  for( i = 0; i < this->n; i++ ) {
    if( 0.0 == map[i] ) v[i] = 0.0;
  }
}

template<typename T>
void SimpleVectorBase<T>::addSomeConstants( T c, const OoqpVectorBase<T>& select )
{
  const SimpleVectorBase<T> & sselect = dynamic_cast<const SimpleVectorBase<T> &>(select);
  T * map = sselect.v;

  int i;
  assert( this->n == select.length() );
  for( i = 0; i < this->n; i++ ) {
    if( map[i] ) v[i] += c;
  }
}

template<typename T>
bool SimpleVectorBase<T>::somePositive( const OoqpVectorBase<T>& select ) const
{
  const SimpleVectorBase<T> & sselect = dynamic_cast<const SimpleVectorBase<T> &>(select);
  T * map = sselect.v;

  assert( this->n == select.length() );

  int i;
  for( i = 0; i < this->n; i++ ) {
    if( 0.0 != map[i] && v[i] <= 0 ) {
      std::cout << "Element " << i << " is nonpositive: " << v[i] << std::endl;
      return false;
    }
  }
  return true;
}

template<typename T>
void SimpleVectorBase<T>::divideSome( const OoqpVectorBase<T>& div, const OoqpVectorBase<T>& select )
{
  if( this->n == 0 ) return;

  const SimpleVectorBase<T> & sselect = dynamic_cast<const SimpleVectorBase<T> &>(select);
  T * map = sselect.v;

  const SimpleVectorBase<T> & sdiv = dynamic_cast<const SimpleVectorBase<T> &>(div);
  T * q   = sdiv.v;
  assert( this->n == div.length() && this->n == select.length() );

  // todo Daniel Rehfeldt: this expects non-aliasing; think there is no need for pointer arithmetic
  // as the run-time is dominated by the division.
#if 0
  T * lmap = map + this->n;
  T * w = v;
  while( map < lmap ) {
    if( 0 != *map ) {
      *w  /= *q;
    }
    map++;
    w++;
    q++;
  }
#else
  for( int i = 0; i < this->n; i++ )
  {
     if( 0.0 != map[i] )
     {
        assert(q[i] != 0.0);
        v[i] /= q[i];
     }
  }
#endif

}

template<typename T>
void SimpleVectorBase<T>::removeEntries(const OoqpVectorBase<T>& select)
{
   const SimpleVectorBase<T>& selectSimple = dynamic_cast<const SimpleVectorBase<T>&>(select);
   const T* const selectArr = selectSimple.v;

   assert(this->n == selectSimple.length());

   int nNew = 0;

   for( int i = 0; i < this->n; i++ )
      if( selectArr[i] != 0.0 )
         v[nNew++] = v[i];

   this->n = nNew;
}

template<typename T>
void SimpleVectorBase<T>::permuteEntries(const std::vector<unsigned int>& permvec)
{
   if( this->n == 0 )
      return;

   assert(this->n > 0);
   assert(permvec.size() == size_t(this->n));

   T* buffer = new T[this->n];

   for( size_t i = 0; i < permvec.size(); i++ )
   {
      assert(permvec[i] < unsigned(this->n));
      buffer[i] = v[permvec[i]];
   }

   std::swap(v, buffer);

   delete[] buffer;
}

template class SimpleVectorBase<int>;
// template class SimpleVectorBase<bool>;
template class SimpleVectorBase<double>;
