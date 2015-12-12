/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */
/* Modified by Cosmin Petra and Miles Lubin */
/* 2015. Modified by Nai-Yuan Chiang for NLP*/


#ifndef OOQPVECTOR_H
#define OOQPVECTOR_H

#include <cstring>
#include <iostream>
#include <fstream>
#include "IotrRefCount.h"
#include "OoqpVectorHandle.h"

//using namespace std;

/** Norm types */
enum PIPS_NormType {
  PIPS_NORM_INF = -1,
  PIPS_NORM_ONE = 0,
  PIPS_NORM_TWO 
};


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

  double Norm(const int Norm_Type, OoqpVector *vec2 = NULL);
  double Norm(const int Norm_Type, const double norm1, const double norm2);

  /** Set all elements of this OoqpVector to zero. */
  virtual void setToZero() = 0;
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
  virtual double infnorm() = 0;
  /** Return the one norm of this OoqpVector object. */
  virtual double onenorm() = 0;

  /** Multiply the components of this OoqpVector by the components of v. */
  virtual void componentMult( OoqpVector& v ) = 0;
  /** Divide the components of this OoqpVector by the components of v. */
  virtual void componentDiv ( OoqpVector& v ) = 0;
  /** Write the components of this OoqpVector, one element per line, to 
   *  the stream out. 
   */
  /** Multiply the components of this OoqpVector by num. */
  virtual void scalarMult( double num) = 0;

  virtual void writeToStream(std::ostream& out) const = 0;
  /** Write the components of this OoqpVector to a stream, subject to 
   *  a format.
   *  @param out a C++-style output stream
   *  @param format a string used to format the output. The substring
   *         %{index} will be substituted by the index of the current element
   *         and the string %{value} will be substituted with the element's
   *         value.
   */
  virtual void writefToStream(std::ostream& out,
			       const char format[] ) const = 0;

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
  virtual void min( double& m, int& index ) = 0;
  /** Return the dot product of this OoqpVector with v */
  virtual double dotProductWith( OoqpVector& v ) = 0;

  /** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
   */
  virtual double shiftedDotProductWith( double alpha, OoqpVector& mystep,
					OoqpVector& yvec,
					double beta,  OoqpVector& ystep ) = 0;
  /** Negate all the elements of this OoqpVector object. */
  virtual void negate() = 0;

  /** Invert (1/x) the elements of this OoqpVector. */
  virtual void invert() = 0;

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

  virtual void writefSomeToStream( std::ostream& out,
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
  virtual int isKindOf( int kind ) = 0;

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

  /** Copy the elements of this OoqpVector into the C-style array v. */
  virtual void copyIntoArray( double v[] ) = 0;
  /** Copy the elements of the C-style array v into this OoqpVector. */
  virtual void copyFromArray( double v[] ) = 0;
  /** Copy the elements of the C-style char array v into this OoqpVector. */
  virtual void copyFromArray( char v[] ) = 0;


  /*add by Naiyuan 2013*/
  virtual void findBlockingPD(OoqpVector & wstep_vec, 
				      OoqpVector & u_vec, 
				      OoqpVector & ustep_vec, 
				      double *w_elt, 
				      double *wstep_elt,
				      double *u_elt, 
				      double *ustep_elt,
				      int& first_or_second, double * alphaPri, double * alphaDual) = 0;  
  virtual void print()=0;
  // compute the sum of the log function (of each component)
  virtual double sumLog(OoqpVector* select)=0;
  // 2-norm
  virtual double twonorm()=0;
  // sum of the power of elts, can be used to evaluate norm
  virtual double sumPowElt(const double pow_in)=0;
  // sum of elts
  virtual double sumElt()=0;

  virtual void MinComponentOrConstant( OoqpVector* vec_in, double minVal )=0;
  
  virtual void MinComponentBetween( OoqpVector* vec_in, OoqpVector *select=NULL)=0; 

  virtual void correctLargeVal( const double testVal, const double corVal, const int absFlag)=0;

  virtual void setToConstantFromTo( double c, int start, int end ) = 0; 

  virtual void SetComponentFromMaxXorY( OoqpVector* x_in, OoqpVector *y_in, OoqpVector *select=NULL)=0;

  virtual void SetComponentFromMaxXorY( OoqpVector* x_in, OoqpVector *y_in, 
			int Start, int End,int xStart, int xEnd, int yStart, int yEnd)=0;

  virtual void SetComponentFromMaxXorConstant( OoqpVector* x_in, const double y_in)=0;

  virtual void SetComponentFromMaxXorConstant( OoqpVector* x_in, const double y_in, 
		  int Start, int End,int xStart, int xEnd)=0;

  virtual void copyFromFromTo( OoqpVector* vec_in, int VStart, int VEnd, int VinStart, int VinEnd)=0;

  virtual void copyFromArrayFromTo( double *w, int VStart, int VEnd, int Wstart, int Wend)=0;


  /** Return the absolute value */
  virtual void absVal(OoqpVector *vec_in) = 0;

  virtual void copyIntoArrayWithIndex_AggVarCon( double *ResultArray, const int *vecmap,const int _length, bool isVar=true){assert("not done"&&0);};

};


enum { kSimpleVector = 0, kPetscVector, kStochVector, kStochDummy, 
	kScaVector, kEmtlVector};

#endif

