/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef SIMPLEVECTOR
#define SIMPLEVECTOR

#include "OoqpVector.h"
#include "SimpleVectorHandle.h"

/**
 * Simple sequential vectors with element access.
 * @ingroup SparseLinearAlgebra
 * @ingroup DenseLinearAlgebra
 */
class SimpleVector : public OoqpVector {
protected:
  int preserveVec;
  double * v;
public:
  SimpleVector( int nx );
  SimpleVector( double * v, int nx );
  //@{
  /**
   * Access the individual elements of this vector.
   */
  double & operator[]( int i ) { 
#ifdef RANGECHECKS
    assert( i > 0 && i < n );
#endif
    return v[i];
 }
  const double & operator[]( int i ) const
  {
#ifdef RANGECHECKS
    assert( i > 0 && i < n );
#endif
    return v[i];
  } 
  //@}
  virtual ~SimpleVector();

  virtual void copyIntoArray( double v[] ) ;
  virtual void copyFromArray( double v[] );
  virtual void copyFromArray( char   v[] );
  virtual void setToZero();
  virtual void setToConstant( double c );
  virtual void randomize( double alpha, double beta, double *ix );
  virtual void copyFrom( OoqpVector& v );
  virtual double twonorm();
  virtual double infnorm();
  virtual double onenorm();
  virtual void min( double& m, int& index );
  virtual void max( double& m, int& index );

  virtual void componentMult( OoqpVector& v );
  virtual void scalarMult( double num);
  virtual void printSolutionToStdErr( OoqpVector& v );
  virtual void componentDiv ( OoqpVector& v );
  virtual void writeToStream(ostream& out) const;
  virtual void writefToStream( ostream& out,
			       const char format[] ) const;

  virtual void scale( double alpha );

  virtual void axpy  ( double alpha, OoqpVector& x );
  virtual void axzpy ( double alpha, OoqpVector& x, OoqpVector& z );
  virtual void axdzpy( double alpha, OoqpVector& x, OoqpVector& z );

  virtual void addConstant( double c );

/** perform the projection operation required by Gondzio algorithm:
   * replace each component of the vector v by vp_i - v_i, where vp_i
   * is the projection of v_i onto the box [rmin, rmax]. Then if the
   * resulting value is less than -rmax, replace it by -rmax.
   * */
  virtual void gondzioProjection( double rmin, double rmax );
  virtual double dotProductWith( OoqpVector& v );
  virtual double shiftedDotProductWith( double alpha, OoqpVector& mystep,
					OoqpVector& yvec,
					double beta,  OoqpVector& ystep );
  virtual void negate();
  virtual void invert();
  virtual int allPositive();
  virtual long long numberOfNonzeros();

  virtual int matchesNonZeroPattern( OoqpVector& select );
  virtual void selectNonZeros( OoqpVector& select );
  virtual void addSomeConstants( double c, OoqpVector& select );
  virtual void writefSomeToStream( ostream& out,
				   const char format[],
				   OoqpVector& select ) const;
  virtual void axdzpy( double alpha, OoqpVector& x,
		       OoqpVector& z, OoqpVector& select );

  virtual int isKindOf( int kind );

  virtual int somePositive( OoqpVector& select );
  virtual void divideSome( OoqpVector& div, OoqpVector& select );

  virtual double stepbound(OoqpVector & v, double maxStep  );
  virtual double findBlocking(OoqpVector & wstep_vec, 
			      OoqpVector & u_vec, 
			      OoqpVector & ustep_vec, 
			      double maxStep,
			      double *w_elt, 
			      double *wstep_elt,
			      double *u_elt, 
			      double *ustep_elt,
			      int& first_or_second);
  /** Returns a pointer to the elemens of this vector. */
  double * elements() { return v; };



  /* following routines are added by Naiyuan 2013 */

  virtual void print();

  virtual void findBlockingPD(OoqpVector & wstep_vec, 
				      OoqpVector & u_vec, 
				      OoqpVector & ustep_vec, 
				      double *w_elt, 
				      double *wstep_elt,
				      double *u_elt, 
				      double *ustep_elt,
				      int& first_or_second, double * alphaPri, double * alphaDual);


  // log function
  virtual double sumLog(OoqpVector* select);
  // sum of the power of elts, can be used to evaluate norm
  virtual double sumPowElt(const double pow_in);
  // sum of elts
  virtual double sumElt();


  virtual void MinComponentOrConstant( OoqpVector* vec_in, double minVal  );

  virtual void correctLargeVal( const double testVal, const double corVal, const int absFlag);  
  
  virtual void MinComponentBetween( OoqpVector* vec_in, OoqpVector *select_in=NULL);

  virtual void setToConstantFromTo( double c, int start, int end );  

  virtual void SetComponentFromMaxXorY( OoqpVector* x_in, OoqpVector *y_in, OoqpVector *select=NULL);

  virtual void SetComponentFromMaxXorY( OoqpVector* x_in, OoqpVector *y_in, 
			  int Start, int End,int xStart, int xEnd, int yStart, int yEnd);

  virtual void SetComponentFromMaxXorConstant( OoqpVector* x_in, const double y_in);

  virtual void SetComponentFromMaxXorConstant( OoqpVector* x_in, const double y_in, 
			int Start, int End,int xStart, int xEnd);


  virtual void copyFromFromTo( OoqpVector* vec_in, int VStart, int VEnd, int VinStart, int VinEnd);

  virtual void copyFromArrayFromTo( double *w, int VStart, int VEnd, int Wstart, int Wend);

  virtual void copyIntoArrayFromTo( double *w, int VStart, int VEnd, int Wstart, int Wend);

  virtual void absVal(OoqpVector *vec_in);

  virtual void copyIntoArrayWithIndex_AggVarCon( double *ResultArray, const int *vecmap, const int _length, bool isVar=true);  
  
};

#endif
