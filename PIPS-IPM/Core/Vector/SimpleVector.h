/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SIMPLEVECTOR
#define SIMPLEVECTOR

#include "OoqpVector.h"
#include "SimpleVector_fwd.h"
#include "SimpleVectorHandle.h"
#include "pipsdef.h"

#include <vector>

/**
 * Simple sequential vectors with element access.
 * @ingroup SparseLinearAlgebra
 * @ingroup DenseLinearAlgebra
 */
template <typename T>
class SimpleVectorBase : public OoqpVectorBase<T> {
protected:
  int preserveVec;
  T * v;
public:
  SimpleVectorBase( int nx = 0 );
  SimpleVectorBase( T * v, int nx );
  //@{
  /**
   * Access the individual elements of this vector.
   */
  T & operator[]( int i ) {
#ifdef RANGECHECKS
    assert( i >= 0 && i < n );
#endif
    return v[i];
 }
  const T & operator[]( int i ) const
  {
#ifdef RANGECHECKS
    assert( i >= 0 && i < n );
#endif
    return v[i];
  }
  //@}

  OoqpVectorBase<T>* clone() const override;
  /* copy vector entries as well */
  OoqpVectorBase<T>* cloneFull() const override;

  virtual ~SimpleVectorBase();

  void copyIntoArray( T v[] ) const override;
  void copyFromArray( const T v[] ) override;
  void copyFromArray( const char v[] ) override;
  bool isZero() const override;
  void setToZero() override;
  void setToConstant( T c ) override;
  void randomize( T alpha, T beta, T *ix ) override;
  void copyFrom( const OoqpVectorBase<T>& v ) override;
  void copyFromAbs(const OoqpVectorBase<T>& v) override;
  double twonorm() const override;
  T infnorm() const override;
  T onenorm() const override;
  void min( T& m, int& index ) const override;
  void max( T& m, int& index ) const override;
  void absminVecUpdate(OoqpVectorBase<T>& absminvec) const override;
  void absmaxVecUpdate(OoqpVectorBase<T>& absmaxvec) const override;
  void absmin( T& m) const override;
  void absminNonZero(T& m, T zero_eps) const override;

  void componentMult( const OoqpVectorBase<T>& v ) override;
  void scalarMult( T num) override;
  virtual void printSolutionToStdErr( OoqpVectorBase<T>& v );
  void componentDiv( const OoqpVectorBase<T>& v ) override;
  bool componentEqual( const OoqpVectorBase<T>& vec, T tol) const override;
  void writeToStream( std::ostream& out ) const override;
  void writeToStreamAll( std::ostream& out ) const override;
  void writeToStreamAllStringStream( std::stringstream& sout ) const override;
  void writefToStream( std::ostream& out, const char format[] ) const override;
  void writeMPSformatOnlyRhs(std::ostream& out, const std::string rowName, const OoqpVectorBase<T>* irhs) const override;
  void writeMPSformatBoundsWithVar(std::ostream& out, const std::string varStub, const OoqpVectorBase<T>* ix, bool upperBound) const override;

  void scale( T alpha ) override;
  
  void axpy  ( T alpha, const OoqpVectorBase<T>& x ) override;
  void axzpy ( T alpha, const OoqpVectorBase<T>& x, const OoqpVectorBase<T>& z ) override;
  void axdzpy( T alpha, const OoqpVectorBase<T>& x, const OoqpVectorBase<T>& z ) override;

   void addConstant( T c ) override;

/** perform the projection operation required by Gondzio algorithm:
   * replace each component of the vector v by vp_i - v_i, where vp_i
   * is the projection of v_i onto the box [rmin, rmax]. Then if the
   * resulting value is less than -rmax, replace it by -rmax.
   * */
  void gondzioProjection( T rmin, T rmax ) override;
  T dotProductWith( const OoqpVectorBase<T>& v ) const override;
  T dotProductSelf( T scaleFactor ) const override;

  T shiftedDotProductWith( T alpha, const OoqpVectorBase<T>& mystep,
					const OoqpVectorBase<T>& yvec,
					T beta,  const OoqpVectorBase<T>& ystep ) const override;
   void negate() override;
   void invert() override;
   void invertSave( T zeroReplacementVal = 0.0 ) override;
   void applySqrt() override;
   void roundToPow2() override;
   bool allPositive() const override;
   long long numberOfNonzeros() const override;

  bool matchesNonZeroPattern( const OoqpVectorBase<T>& select ) const override;
  void selectNonZeros( const OoqpVectorBase<T>& select ) override;
  void addSomeConstants( T c, const OoqpVectorBase<T>& select ) override;
  void writefSomeToStream( std::ostream& out, const char format[],
				   const OoqpVectorBase<T>& select ) const override;
  void axdzpy( T alpha, const OoqpVectorBase<T>& x,
		       const OoqpVectorBase<T>& z, const OoqpVectorBase<T>& select ) override;

  bool isKindOf( int kind ) const override;

  bool somePositive( const OoqpVectorBase<T>& select ) const override;
  void divideSome( const OoqpVectorBase<T>& div, const OoqpVectorBase<T>& select ) override;

  T stepbound(const OoqpVectorBase<T>& v, T maxStep ) const override;
  T findBlocking(const OoqpVectorBase<T> & wstep_vec,
			      const OoqpVectorBase<T> & u_vec,
			      const OoqpVectorBase<T> & ustep_vec,
			      T maxStep,
			      T *w_elt,
			      T *wstep_elt,
			      T *u_elt,
			      T *ustep_elt,
			      int& first_or_second) const override;

  void findBlocking_pd(const OoqpVectorBase<T> & wstep_vec,
                      const OoqpVectorBase<T> & u_vec, const OoqpVectorBase<T> & ustep_vec,
                      T& maxStepPri, T& maxStepDual,
                      T& w_elt_p, T& wstep_elt_p, T& u_elt_p, T& ustep_elt_p,
                      T& w_elt_d, T& wstep_elt_d, T& u_elt_d, T& ustep_elt_d,
                      bool& primalBlocking, bool& dualBlocking) const override;

  void removeEntries(const OoqpVectorBase<int>& select) override;

  void permuteEntries(const std::vector<unsigned int>& permvec);

  /** Returns a pointer to the elements of this vector. */
  T * elements() const { return v; };
};

#endif
