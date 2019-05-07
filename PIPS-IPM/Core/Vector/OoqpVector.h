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
using namespace std;
#include "IotrRefCount.h"
#include "OoqpVectorHandle.h"
#include "pipsdef.h"

/** An abstract class representing the implementation of a OoqpVector.
 *
 *  Do not create instances of OoqpVector. Create instance of subclasses
 *  of OoqpVector instead.
 *
 *  @ingroup AbstractLinearAlgebra
 */
class OoqpVector : public IotrRefCount {
public:
  int n;
  /** Return the length of this vector. */
  long long length() const { return n; }

  OoqpVector( int n_ = 0 );
  virtual ~OoqpVector();

  /** Set all elements of this OoqpVector to zero. */
  virtual void setToZero() = 0;
  /** Check if all elemens in the vector are equal to zero. */
  virtual bool isZero() const = 0;
  /** Set all elements of this OoqpVector to the constant value c */
  virtual void setToConstant( double c ) = 0;
  /** Fill this OoqpVector with random elements
   *	  @param alpha
   *      @param beta the elements will be in the interval [alpha, beta]
   *      @param ix an aribitray number used to seed the random number
   *             generator
   */
  virtual void randomize( double alpha, double beta, double *ix ) = 0;
  /** Copy the elements of v into this OoqpVector object. */
  virtual void copyFrom( OoqpVector& v ) = 0;

  /** Return the infinity norm of this OoqpVector object. */
  virtual double twonorm() const = 0;
  /** Return the infinity norm of this OoqpVector object. */
  virtual double infnorm() const = 0;
  /** Return the one norm of this OoqpVector object. */
  virtual double onenorm() const = 0;

  /** Multiply the components of this OoqpVector by the components of v. */
  virtual void componentMult( OoqpVector& v ) = 0;
  /** Divide the components of this OoqpVector by the components of v. */
  virtual void componentDiv ( OoqpVector& v ) = 0;
  /** Write the components of this OoqpVector, one element per line, to
   *  the stream out.
   */
  /** Multiply the components of this OoqpVector by num. */
  virtual void scalarMult( double num) = 0;

  virtual void writeToStreamAll(ostream& out) const {assert(0 && "not implemented here");};
  virtual void writeToStreamAllStringStream( stringstream& sout ) const {assert(0 && "not implemented here");};
  virtual void writeToStreamAllChild( stringstream& sout ) const {assert(0 && "not implemented here");};


  virtual void writeToStream(ostream& out) const = 0;
  /** Write the components of this OoqpVector to a stream, subject to
   *  a format.
   *  @param out a C++-style output stream
   *  @param format a string used to format the output. The substring
   *         %{index} will be substituted by the index of the current element
   *         and the string %{value} will be substituted with the element's
   *         value.
   */
  virtual void writefToStream( ostream& out,
			       const char format[] ) const = 0;

  virtual void writeMPSformatOnlyRhs(ostream& out, string rowName, OoqpVector* irhs) const {assert(0 && "not implemented here");};
  virtual void writeMPSformatRhs(ostream& out, int rowType, OoqpVector* irhs) const {assert(0 && "not implemented here");};
  virtual void writeMPSformatBounds(ostream& out, OoqpVector* ix, bool upperBound) const {assert(0 && "not implemented here");};
  virtual void writeMPSformatBoundsWithVar(ostream& out, string varStub, OoqpVector* ix, bool upperBound) const {assert(0 && "not implemented here");};
  void writefToStreamStats( ostream& out, std::string prestring);

  /** Scale each element of this OoqpVector by the constant alpha */
  virtual void scale( double alpha ) = 0;

  /** this += alpha * x */
  virtual void axpy  ( double alpha, OoqpVector& x ) = 0;
  /** this += alpha * x * z */
  virtual void axzpy ( double alpha, OoqpVector& x, OoqpVector& z ) = 0;
  /** this += alpha * x / z */
  virtual void axdzpy( double alpha, OoqpVector& x, OoqpVector& z ) = 0;

  /** Add c to the elements of this OoqpVector object */
  virtual void addConstant( double c ) = 0;

  /** Perform the projection needed by Gondzio's multiple corrector method.
   *
   * @see SimpleVector::gondzioProjection
   */
  virtual void gondzioProjection( double rmin, double rmax ) = 0;

  /** Return the minimum value in this vector, and the index at
   *  which it occurs. */
  virtual void min( double& m, int& index ) const = 0;

  /** Return the maximum value in this vector, and the index at
   *  which it occurs. */
  virtual void max( double& m, int& index ) const = 0;

  /** Return the absolute minimum value of this vector */
  virtual void absmin(double& m) const = 0;

  virtual void absminVecUpdate(OoqpVector& absminvec) const = 0;

  virtual void absmaxVecUpdate(OoqpVector& absmaxvec) const = 0;

  /** Return the absolute minimum value of this vector bigger than eps_zero, or -1 if none could be found */
  virtual void absminNonZero(double& m, double zero_eps) const = 0;

  /** Return the dot product of this OoqpVector with v */
  virtual double dotProductWith( const OoqpVector& v ) const = 0;

  /** Return the scaled dot product of this (scaled) vector with itself  */
  virtual double dotProductSelf(double scaleFactor) const = 0;

  /** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
   */
  virtual double shiftedDotProductWith( double alpha, OoqpVector& mystep,
					OoqpVector& yvec,
					double beta,  OoqpVector& ystep ) = 0;
  /** Negate all the elements of this OoqpVector object. */
  virtual void negate() = 0;

  /** Invert (1/x) the elements of this OoqpVector. */
  virtual void invert() = 0;

  /** Invert (1/x) the elements of this OoqpVector, but don't divide by zero and replace by zeroReplacementVal instead */
  virtual void invertSave( double zeroReplacementVal = 0.0 ) = 0;

  /** Take the square root of each element of this OoqpVector. */
  virtual void applySqrt() = 0;

  /** Rounds vector entries to nearest power of two values */
  virtual void roundToPow2() = 0;

  /** True if all elements of this OoqpVector are positive. */
  virtual int allPositive() = 0;

  /** Return the number of non-zero elements in this OoqpVector. */
  virtual long long numberOfNonzeros() = 0;

  /** True if this OoqpVector has the same non-zero pattern as select. */
  virtual int matchesNonZeroPattern( OoqpVector& select ) = 0;

  /** Set each element of this OoqpVector to zero if the corresponding
   *  element in select is zero.
   */
  virtual void selectNonZeros( OoqpVector& select ) = 0;
  /** Add the constant c to some of the elements of this OoqpVector
   *  @param c The constant to be added
   *  @param select a mask OoqpVector. The constant c is added to an element
   *         of this OoqpVector only if the corresponding element of select is
   *         non-zero.
   */
  virtual void addSomeConstants( double c, OoqpVector& select ) = 0;

  /** Write some elements of this OoqpVector to a stream, subject to a format.
   *  @param out a C++-style output stream
   *  @param format a string used to format the output. The substring
   *         %{index} will be substituted by the index of the current element
   *         and the string %{value} will be substituted with the element's
   *         value.
   *  @param select a mask OoqpVector. An element if this OoqpVector is
   *         written only if the corresponding element of selec is non-zero.
   */

  virtual void writefSomeToStream( ostream& out,
				   const char format[],
				   OoqpVector& select ) const = 0;
  /** this += alpha * x / z
   *  @param select only perform the division on elements of x and z if the
   *         corresponding element of select is non-zero. Generally we avoid
   *         performing the division if we know that it will result in
   *         division by zero. The OoqpVector select may be x, z or a third
   *         OoqpVector.
   */
  virtual void axdzpy( double alpha, OoqpVector& x,
		       OoqpVector& z, OoqpVector& select ) = 0;

  /** True if selected elements of this OoqpVector are positive
   *  @param select Each element of this OoqpVector must be positive
   *                if the corresponding element of select is non-zero.
   */
  virtual int somePositive( OoqpVector& select ) = 0;
  /** Divide selected elements of this OoqpVector by the corresponding
   *  element in div.
   *  @param div If element i of this OoqpVector is selected, then it will
   *             be divided by element i of div.
   *  @param select Perform division on elements of this OoqpVector only if
   *             the corresponding element in select is non-zero
   */
  virtual void divideSome( OoqpVector& div, OoqpVector& select ) = 0;

  /** True if this OoqpVector identifies itself as having the type kind.
   *
   *  Classes overriding this method must call the inherited version, so that
   *  the class hierarchy is supported.
   */
  virtual int isKindOf( int kind ) const = 0;

  /** Return the largest value of alpha in the interval [0, maxStep] for
   *  which: this + alpha * v >= 0. Set firstBlocking to be the
   *  "blocking" index i - the one that limits the step length to alpha.
   */
  virtual double stepbound(OoqpVector & v, double maxStep ) = 0;

  /** Return the largest value of alpha in the interval [0,1] for
   *  which: this + alpha * wstep_vec >= 0 and u_vec + alpha *
   *  ustep_vec >= 0. Also return the components this[i],
   *  wstep_vec[i], u_vec[i], ustep_vec[i], where i is the index of
   *  the "blocking" component - the one that limits the step length
   *  to alpha. Also return first_or_second=1 if the blocking
   *  component is in "this", and first_or_second=2 if the blocking
   *  component is in u_vec.
   */
  virtual double findBlocking(OoqpVector & wstep_vec,
			      OoqpVector & u_vec,
			      OoqpVector & ustep_vec,
			      double maxStep,
			      double *w_elt,
			      double *wstep_elt,
			      double *u_elt,
			      double *ustep_elt,
			      int& first_or_second) = 0;

  virtual void findBlocking_pd(const OoqpVector & wstep_vec,
    			      const OoqpVector & u_vec,
    			      const OoqpVector & ustep_vec,
    			      double& maxStepPri, double& maxStepDual,
    			      double& w_elt_p,
    			      double& wstep_elt_p,
    			      double& u_elt_p,
    			      double& ustep_elt_p,
    				  double& w_elt_d, double& wstep_elt_d, double& u_elt_d, double& ustep_elt_d,
					  bool& primalBlocking, bool& dualBlocking) const = 0;

  /** Copy the elements of this OoqpVector into the C-style array v. */
  virtual void copyIntoArray( double v[] ) const = 0;
  /** Copy the elements of the C-style array v into this OoqpVector. */
  virtual void copyFromArray( double v[] ) = 0;
  /** Copy the elements of the C-style char array v into this OoqpVector. */
  virtual void copyFromArray( char v[] ) = 0;

  /** remove entries i for which select[i] == 0 */
  virtual void removeEntries( const OoqpVector& select ) { assert(0 && "not implemented here"); };

  /** Copy the absolute values of elements of v_in into this OoqpVector object. */
  virtual void copyFromAbs(const OoqpVector& v) = 0;

  virtual OoqpVector* clone() const { assert(0 && "not implemented here"); return NULL; };
  virtual OoqpVector* cloneFull() const { assert(0 && "not implemented here"); return NULL; };
};


enum { kSimpleVector = 0, kPetscVector, kStochVector, kStochDummy,
	kScaVector, kEmtlVector};

#endif
