/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */
/* Modified by Cosmin Petra and Miles Lubin */


#ifndef OOQPVECTOR_H
#define OOQPVECTOR_H

#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>

#include "OoqpVector_fwd.h"

#include "IotrRefCount.h"
#include "OoqpVectorHandle.h"
#include "pipsdef.h"



/** An abstract class representing the implementation of a OoqpVectorTemplate.
 *
 *  Do not create instances of OoqpVectorBase. Create instance of subclasses
 *  of OoqpVectorBase instead.
 *
 *  @ingroup AbstractLinearAlgebra
 */
template<typename T>
class OoqpVectorBase : public IotrRefCount {
public:
  int n;
  /** Return the length of this vector. */
  int length() const { return n; }

  OoqpVectorBase( int n_ = 0 );
  virtual ~OoqpVectorBase();

  /** Set all elements of this OoqpVector to zero. */
  virtual void setToZero() = 0;
  /** Check if all elemens in the vector are equal to zero. */
  virtual bool isZero() const = 0;
  /** Set all elements of this OoqpVector to the constant value c */
  virtual void setToConstant( T c ) = 0;
  /** Fill this OoqpVector with random elements
   *	  @param alpha
   *      @param beta the elements will be in the interval [alpha, beta]
   *      @param ix an aribitray number used to seed the random number
   *             generator
   */
  virtual void randomize( T alpha, T beta, T *ix ) = 0;
  /** Copy the elements of v into this OoqpVector object. */
  virtual void copyFrom( const OoqpVectorBase<T>& v ) = 0;

  /** Return the infinity norm of this OoqpVector object. */
  virtual double twonorm() const = 0;
  /** Return the infinity norm of this OoqpVector object. */
  virtual T infnorm() const = 0;
  /** Return the one norm of this OoqpVector object. */
  virtual T onenorm() const = 0;

  /** Multiply the components of this OoqpVector by the components of v. */
  virtual void componentMult( const OoqpVectorBase<T>& v ) = 0;
  /** Divide the components of this OoqpVector by the components of v. */
  virtual void componentDiv ( const OoqpVectorBase<T>& v ) = 0;
  /** Write the components of this OoqpVector, one element per line, to
   *  the stream out.
   */
  /** Multiply the components of this OoqpVector by num. */
  virtual void scalarMult( T num ) = 0;

  virtual void writeToStreamAll( std::ostream& out ) const {assert(0 && "not implemented here");};
  virtual void writeToStreamAllStringStream( std::stringstream& sout ) const {assert(0 && "not implemented here");};
  virtual void writeToStreamAllChild( std::stringstream& sout ) const {assert(0 && "not implemented here");};


  virtual void writeToStream( std::ostream& out) const = 0;
  /** Write the components of this OoqpVector to a stream, subject to
   *  a format.
   *  @param out a C++-style output stream
   *  @param format a string used to format the output. The substring
   *         %{index} will be substituted by the index of the current element
   *         and the string %{value} will be substituted with the element's
   *         value.
   */
  virtual void writefToStream( std::ostream& out, const char format[] ) const = 0;

  virtual void writeMPSformatOnlyRhs(std::ostream& out, const std::string rowName, const OoqpVectorBase<T>* irhs) const {assert(0 && "not implemented here");};
  virtual void writeMPSformatRhs(std::ostream& out, int rowType, const OoqpVectorBase<T>* irhs) const {assert(0 && "not implemented here");};
  virtual void writeMPSformatBounds(std::ostream& out, const OoqpVectorBase<T>* ix, bool upperBound) const {assert(0 && "not implemented here");};
  virtual void writeMPSformatBoundsWithVar(std::ostream& out, const std::string varStub, const OoqpVectorBase<T>* ix, bool upperBound) const {assert(0 && "not implemented here");};
  void writefToStreamStats( std::ostream& out, const std::string prestring);

  /** Scale each element of this OoqpVector by the constant alpha */
  virtual void scale( T alpha ) = 0;

  /** this += alpha * x */
  virtual void axpy  ( T alpha, const OoqpVectorBase<T>& x ) = 0;
  /** this += alpha * x * z */
  virtual void axzpy ( T alpha, const OoqpVectorBase<T>& x, const OoqpVectorBase<T>& z ) = 0;
  /** this += alpha * x / z */
  virtual void axdzpy( T alpha, const OoqpVectorBase<T>& x, const OoqpVectorBase<T>& z ) = 0;

  /** Add c to the elements of this OoqpVector object */
  virtual void addConstant( T c ) = 0;

  /** Perform the projection needed by Gondzio's multiple corrector method.
   *
   * @see SimpleVector::gondzioProjection
   */
  virtual void gondzioProjection( T rmin, T rmax ) = 0;

  /** Return the minimum value in this vector, and the index at
   *  which it occurs. */
  virtual void min( T& m, int& index ) const = 0;

  /** Return the maximum value in this vector, and the index at
   *  which it occurs. */
  virtual void max( T& m, int& index ) const = 0;

  /** Return the absolute minimum value of this vector */
  virtual void absmin(T& m) const = 0;

  virtual void absminVecUpdate(OoqpVectorBase<T>& absminvec) const = 0;

  virtual void absmaxVecUpdate(OoqpVectorBase<T>& absmaxvec) const = 0;

  /** Return the absolute minimum value of this vector bigger than eps_zero, or -1 if none could be found */
  virtual void absminNonZero(T& m, T zero_eps) const = 0;

  /** Return the dot product of this OoqpVector with v */
  virtual T dotProductWith( const OoqpVectorBase<T>& v ) const = 0;

  /** Return the scaled dot product of this (scaled) vector with itself  */
  virtual T dotProductSelf(T scaleFactor) const = 0;

  /** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
   */
  virtual T shiftedDotProductWith( T alpha, const OoqpVectorBase<T>& mystep,
					const OoqpVectorBase<T>& yvec,
					T beta,  const OoqpVectorBase<T>& ystep ) const = 0;
  /** Negate all the elements of this OoqpVector object. */
  virtual void negate() = 0;

  /** Invert (1/x) the elements of this OoqpVector. */
  virtual void invert() = 0;

  /** Invert (1/x) the elements of this OoqpVector, but don't divide by zero and replace by zeroReplacementVal instead */
  virtual void invertSave( T zeroReplacementVal = 0.0 ) = 0;

  /** Take the square root of each element of this OoqpVector. */
  virtual void applySqrt() = 0;

  /** Rounds vector entries to nearest power of two values */
  virtual void roundToPow2() = 0;

  /** True if all elements of this OoqpVector are positive. */
  virtual bool allPositive() const = 0;

  /** Return the number of non-zero elements in this OoqpVector. */
  virtual int numberOfNonzeros() const = 0;

  /** True if this OoqpVector has the same non-zero pattern as select. */
  virtual bool matchesNonZeroPattern( const OoqpVectorBase<T>& select ) const = 0;

  /** Set each element of this OoqpVector to zero if the corresponding
   *  element in select is zero.
   */
  virtual void selectNonZeros( const OoqpVectorBase<T>& select ) = 0;
  
  /** Add the constant c to some of the elements of this OoqpVector
   *  @param c The constant to be added
   *  @param select a mask OoqpVector. The constant c is added to an element
   *         of this OoqpVector only if the corresponding element of select is
   *         non-zero.
   */
  virtual void addSomeConstants( T c, const OoqpVectorBase<T>& select ) = 0;

  /** Write some elements of this OoqpVector to a stream, subject to a format.
   *  @param out a C++-style output stream
   *  @param format a string used to format the output. The substring
   *         %{index} will be substituted by the index of the current element
   *         and the string %{value} will be substituted with the element's
   *         value.
   *  @param select a mask OoqpVector. An element if this OoqpVector is
   *         written only if the corresponding element of selec is non-zero.
   */

  virtual void writefSomeToStream( std::ostream& out,
				   const char format[],
				   const OoqpVectorBase<T>& select ) const = 0;
  /** this += alpha * x / z
   *  @param select only perform the division on elements of x and z if the
   *         corresponding element of select is non-zero. Generally we avoid
   *         performing the division if we know that it will result in
   *         division by zero. The OoqpVector select may be x, z or a third
   *         OoqpVector.
   */
  virtual void axdzpy( T alpha, const OoqpVectorBase<T>& x,
		       const OoqpVectorBase<T>& z, const OoqpVectorBase<T>& select ) = 0;

  /** True if selected elements of this OoqpVector are positive
   *  @param select Each element of this OoqpVector must be positive
   *                if the corresponding element of select is non-zero.
   */
  virtual bool somePositive( const OoqpVectorBase<T>& select ) const = 0;
  /** Divide selected elements of this OoqpVector by the corresponding
   *  element in div.
   *  @param div If element i of this OoqpVector is selected, then it will
   *             be divided by element i of div.
   *  @param select Perform division on elements of this OoqpVector only if
   *             the corresponding element in select is non-zero
   */
  virtual void divideSome( const OoqpVectorBase<T>& div, const OoqpVectorBase<T>& select ) = 0;

  /** True if this OoqpVector identifies itself as having the type kind.
   *
   *  Classes overriding this method must call the inherited version, so that
   *  the class hierarchy is supported.
   */
  virtual bool isKindOf( int kind ) const = 0;

  /** Return the largest value of alpha in the interval [0, maxStep] for
   *  which: this + alpha * v >= 0. Set firstBlocking to be the
   *  "blocking" index i - the one that limits the step length to alpha.
   */
  virtual T stepbound(const OoqpVectorBase<T> & v, T maxStep ) const = 0;

  /** Return the largest value of alpha in the interval [0,1] for
   *  which: this + alpha * wstep_vec >= 0 and u_vec + alpha *
   *  ustep_vec >= 0. Also return the components this[i],
   *  wstep_vec[i], u_vec[i], ustep_vec[i], where i is the index of
   *  the "blocking" component - the one that limits the step length
   *  to alpha. Also return first_or_second=1 if the blocking
   *  component is in "this", and first_or_second=2 if the blocking
   *  component is in u_vec.
   */
  virtual T findBlocking(const OoqpVectorBase<T> & wstep_vec,
			      const OoqpVectorBase<T> & u_vec,
			      const OoqpVectorBase<T> & ustep_vec,
			      T maxStep,
			      T *w_elt,
			      T *wstep_elt,
			      T *u_elt,
			      T *ustep_elt,
			      int& first_or_second) const = 0;

  virtual void findBlocking_pd(const OoqpVectorBase<T> & wstep_vec,
    			      const OoqpVectorBase<T> & u_vec,
    			      const OoqpVectorBase<T> & ustep_vec,
    			      T& maxStepPri, T& maxStepDual,
    			      T& w_elt_p,
    			      T& wstep_elt_p,
    			      T& u_elt_p,
    			      T& ustep_elt_p,
    				  T& w_elt_d, T& wstep_elt_d, T& u_elt_d, T& ustep_elt_d,
					  bool& primalBlocking, bool& dualBlocking) const = 0;

  /** Copy the elements of this OoqpVector into the C-style array v. */
  virtual void copyIntoArray( T v[] ) const = 0;
  /** Copy the elements of the C-style array v into this OoqpVector. */
  virtual void copyFromArray( const T v[] ) = 0;
  /** Copy the elements of the C-style char array v into this OoqpVector. */
  virtual void copyFromArray( const char v[] ) = 0;

  /** remove entries i for which select[i] == 0 */
  virtual void removeEntries( const OoqpVectorBase<T>& select ) { assert(0 && "not implemented here"); };

  /** Copy the absolute values of elements of v_in into this OoqpVector object. */
  virtual void copyFromAbs(const OoqpVectorBase<T>& v) = 0;

  virtual OoqpVectorBase<T>* clone() const { assert(0 && "not implemented here"); return NULL; };
  virtual OoqpVectorBase<T>* cloneFull() const { assert(0 && "not implemented here"); return NULL; };
};


enum { kSimpleVector = 0, kPetscVector, kStochVector, kStochDummy,
	kScaVector, kEmtlVector};

#include "OoqpVector.C"

#endif
