/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SIMPLEVECTOR
#define SIMPLEVECTOR


#include "OoqpVector.h"
#include "SimpleVectorHandle.h"
#include "pipsdef.h"
#include "vector"

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
  SimpleVector( int nx = 0);
  SimpleVector( double * v, int nx );
  //@{
  /**
   * Access the individual elements of this vector.
   */
  double & operator[]( int i ) {
#ifdef RANGECHECKS
    assert( i >= 0 && i < n );
#endif
    return v[i];
 }
  const double & operator[]( int i ) const
  {
#ifdef RANGECHECKS
    assert( i >= 0 && i < n );
#endif
    return v[i];
  }
  //@}

  /* copy vector entries as well */
  virtual SimpleVector* cloneFull() const;

  virtual ~SimpleVector();

  virtual void copyIntoArray( double v[] ) const;
  virtual void copyFromArray( double v[] );
  virtual void copyFromArray( char   v[] );
  virtual bool isZero() const;
  virtual void setToZero() { setToConstant(0.0); };
  virtual void setToConstant( double c );
  virtual void randomize( double alpha, double beta, double *ix );
  virtual void copyFrom( OoqpVector& v );
  virtual void copyFromAbs(const OoqpVector& v);
  virtual double twonorm() const;
  virtual double infnorm() const;
  virtual double onenorm() const ;
  virtual void min( double& m, int& index ) const;
  virtual void max( double& m, int& index ) const;
  virtual void absminVecUpdate(OoqpVector& absminvec) const;
  virtual void absmaxVecUpdate(OoqpVector& absmaxvec) const;
  virtual void absmin( double& m) const;
  virtual void absminNonZero(double& m, double zero_eps) const;

  virtual void componentMult( OoqpVector& v );
  virtual bool componentEqual( const OoqpVector& vec, double tol) const;
  virtual void scalarMult( double num);
  virtual void printSolutionToStdErr( OoqpVector& v );
  virtual void componentDiv ( OoqpVector& v );
  virtual void writeToStream(ostream& out) const;
  virtual void writeToStreamAll(ostream& out) const;
  virtual void writeToStreamAllStringStream(stringstream& sout) const;
  virtual void writefToStream( ostream& out,
			       const char format[] ) const;
  virtual void writeMPSformatOnlyRhs(ostream& out, string rowName, OoqpVector* irhs) const;
  virtual void writeMPSformatBoundsWithVar(ostream& out, string varStub, OoqpVector* ix, bool upperBound) const;

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
  virtual double dotProductWith( const OoqpVector& v ) const;
  virtual double dotProductSelf( double scaleFactor ) const;

  virtual double shiftedDotProductWith( double alpha, OoqpVector& mystep,
					OoqpVector& yvec,
					double beta,  OoqpVector& ystep );
  virtual void negate();
  virtual void invert();
  virtual void invertSave( double zeroReplacementVal = 0.0 );
  virtual void applySqrt();
  virtual void roundToPow2();
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

  virtual int isKindOf( int kind ) const;

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

  virtual void findBlocking_pd(const OoqpVector & wstep_vec,
                      const OoqpVector & u_vec, const OoqpVector & ustep_vec,
                      double& maxStepPri, double& maxStepDual,
                      double& w_elt_p, double& wstep_elt_p, double& u_elt_p, double& ustep_elt_p,
                      double& w_elt_d, double& wstep_elt_d, double& u_elt_d, double& ustep_elt_d,
                      bool& primalBlocking, bool& dualBlocking) const;

  virtual void removeEntries(const OoqpVector& select);

  void permuteEntries(const std::vector<unsigned int>& permvec);

  /** Returns a pointer to the elements of this vector. */
  double * elements() const { return v; };
};

#endif
