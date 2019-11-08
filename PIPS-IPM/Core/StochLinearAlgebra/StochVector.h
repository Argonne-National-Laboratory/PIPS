#ifndef STOCHVECTOR_H
#define STOCHVECTOR_H

#include "OoqpVector.h"
#include "SimpleVector.h"
#include "StochVector_fwd.h"
#include "StochVectorHandle.h"
#include "mpi.h"

#include <vector>

class StochTree;

template <typename T>
class StochVectorBase : public OoqpVectorBase<T> {
private:

  virtual void writeToStreamAllChild( std::stringstream& sout ) const;

public:
  StochVectorBase( int n, MPI_Comm mpiComm, int isDistributed = -1);
  StochVectorBase( int n, int nl, MPI_Comm mpiComm, int isDistributed);
  virtual ~StochVectorBase();

  virtual void AddChild(StochVectorBase<T>* child);
  virtual void AddChild(OoqpVectorBase<T>* child);

  /** The data for this node. */
  OoqpVectorBase<T>*               vec;

  /** The linking constraint data for this node. */
  OoqpVectorBase<T>*               vecl;

  /** Children of this node */
  std::vector<StochVectorBase<T>*> children;

  /** Link to the parent of this node. Needed when we multiply a matrix
      with this vector
  */
  StochVectorBase<T>*              parent;


  /* MPI communicator */
  MPI_Comm mpiComm;
  /* flag used to indicate if the children are distributed or not. */
  int iAmDistrib;

  /** Creates and returns a vector of the type used to store data in this node,
      i.e., same type as 'vec'.
      NO data is copied, for this use one of the 'copy...' functions.
  */
  virtual OoqpVectorBase<T>* dataClone() const;
  virtual OoqpVectorBase<T>* dataCloneLinkCons() const;
  virtual StochVectorBase<T>* clone() const;
  /* copy vector entries as well */
  virtual StochVectorBase<T>* cloneFull() const;

  virtual void jointCopyFrom(const StochVectorBase<T>& v1, const StochVectorBase<T>& v2, const StochVectorBase<T>& v3);
  virtual void jointCopyFromLinkCons(const StochVectorBase<T>& vx, const StochVectorBase<T>& vy, const StochVectorBase<T>& vz);
  virtual void jointCopyTo(const StochVectorBase<T>& v1, const StochVectorBase<T>& v2, const StochVectorBase<T>& v3);
  virtual void jointCopyToLinkCons(const StochVectorBase<T>& vx, const StochVectorBase<T>& vy, const StochVectorBase<T>& vz);

  bool isKindOf( int kind ) const override;
  void setToZero() override;
  void setToConstant( T c ) override;
  bool isZero() const;

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
   void absmin( T& m ) const override;
   void absminNonZero(T& m, T zero_eps) const override;
   T stepbound(const OoqpVectorBase<T> & v, T maxStep ) const override;
   T findBlocking(const OoqpVectorBase<T> & wstep_vec,
			      const OoqpVectorBase<T> & u_vec,
			      const OoqpVectorBase<T> & ustep_vec,
			      T maxStep,
			      T *w_elt,
			      T *wstep_elt,
			      T *u_elt,
			      T *ustep_elt,
			      int& first_or_second) const override;
   void findBlocking_pd(const OoqpVectorBase<T>& wstep_vec,
               const OoqpVectorBase<T>& u_vec,
  			      const OoqpVectorBase<T>& ustep_vec,
  			      T& maxStepPri, T& maxStepDual,
  			      T& w_elt_p, T& wstep_elt_p, T& u_elt_p, T& ustep_elt_p,
  				   T& w_elt_d, T& wstep_elt_d, T& u_elt_d, T& ustep_elt_d,
				   bool& primalBlocking, bool& dualBlocking) const override;

   void componentMult( const OoqpVectorBase<T>& v ) override;
   void componentDiv ( const OoqpVectorBase<T>& v ) override;
   void scalarMult( T num) override;
   void writeToStream(std::ostream& out) const override;
   void writeToStreamAll(std::ostream& out) const override;
   void writefToStream( std::ostream& out,
			       const char format[] ) const override;
   void writeMPSformatOnlyRhs(std::ostream& out, const std::string rowName, const OoqpVectorBase<T>* irhs) const override {};
   void writeMPSformatRhs(std::ostream& out, int rowType, const OoqpVectorBase<T>* irhs) const override;
   void writeMPSformatBounds(std::ostream& out, const OoqpVectorBase<T>* ix, bool upperBound) const override;
   void writeMPSformatBoundsWithVar(std::ostream& out, const std::string varStub, const OoqpVectorBase<T>* ix, bool upperBound) const override {};

   void scale( T alpha ) override;

  /** this += alpha * x */
   void axpy  ( T alpha, const OoqpVectorBase<T>& x ) override;
  /** this += alpha * x * z */
   void axzpy ( T alpha, const OoqpVectorBase<T>& x, const OoqpVectorBase<T>& z ) override;
  /** this += alpha * x / z */
   void axdzpy( T alpha, const OoqpVectorBase<T>& x, const OoqpVectorBase<T>& z ) override;

   void addConstant( T c ) override;
   void gondzioProjection( T rmin, T rmax ) override;
   T dotProductWith( const OoqpVectorBase<T>& v ) const override;
   T dotProductSelf( T scaleFactor) const override;

  /** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
   */
   T shiftedDotProductWith( T alpha, const OoqpVectorBase<T>& mystep,
					const OoqpVectorBase<T>& yvec,
					T beta, const OoqpVectorBase<T>& ystep ) const override;
   void negate() override;
   void invert() override;
   void invertSave( T zeroReplacementVal = 0.0 ) override;
   void applySqrt() override;
   void roundToPow2() override;

   bool allPositive() const override;

   bool matchesNonZeroPattern( const OoqpVectorBase<T>& select ) const override;
   void selectNonZeros( const OoqpVectorBase<T>& select ) override;
   int numberOfNonzeros() const override;
   void addSomeConstants( T c, const OoqpVectorBase<T>& select ) override;
   void writefSomeToStream( std::ostream& out,
				   const char format[],
				   const OoqpVectorBase<T>& select ) const override;
   void axdzpy( T alpha, const OoqpVectorBase<T>& x,
		       const OoqpVectorBase<T>& z, const OoqpVectorBase<T>& select ) override;

   bool somePositive( const OoqpVectorBase<T>& select ) const override;
   void divideSome( const OoqpVectorBase<T>& div, const OoqpVectorBase<T>& select ) override;
   void copyIntoArray( T v[] ) const override;
   void copyFromArray( const T v[] ) override;
   void copyFromArray( const char v[] ) override;
   virtual void permuteVec0Entries(const std::vector<unsigned int>& permvec);
   virtual void permuteLinkingEntries(const std::vector<unsigned int>& permvec);
   virtual std::vector<T> gatherStochVector() const;

   /** remove entries i for which select[i] == 0 */
   void removeEntries( const OoqpVectorBase<T>& select ) override;

   int getSize() const { return this.n; };

   virtual bool isRootNodeInSync() const;

};

/** DUMMY VERSION
 *
 */
template <typename T>
class StochDummyVectorBase : public StochVectorBase<T> {
protected:

public:
  StochDummyVectorBase(  )
    : StochVectorBase<T>(0, MPI_COMM_NULL) {};

  virtual ~StochDummyVectorBase(){};

  void AddChild(StochVector* child) override {};
  void AddChild(OoqpVector* child) override {};

  /** Creates and returns a vector of the type used to store data in this node,
      i.e., same type as 'vec'.
      NO data is copied, for this use one of the 'copy...' functions.
  */
   OoqpVectorBase<T>* dataClone() const override { return new SimpleVectorBase<T>(0);}
   StochVectorBase<T>* clone() const override { return new StochDummyVectorBase<T>();}
   StochVectorBase<T>* cloneFull() const override { return new StochDummyVectorBase<T>();}

   void jointCopyFrom(StochVector& v1, StochVector& v2, StochVector& v3)override {};
   void jointCopyTo(StochVector& v1, StochVector& v2, StochVector& v3)override {};

   int isKindOf( int kind ) const override {return kind == kStochDummy;}
   void setToZero()override {};
   bool isZero() const override { return true; };

   void setToConstant( double c )override {};
   void randomize( double alpha, double beta, double *ix )override {};
   void copyFrom( OoqpVectorBase<T>& v )override {};
   void copyFromAbs(const OoqpVectorBase<T>& v) override {};
   double twonorm() const override {return 0.0;}
   double infnorm() const override {return 0.0;}
   double onenorm() const override {return 0.0;}
   void min( double& m, int& index ) const override {};
   void max( double& m, int& index ) const override {};
   void absminVecUpdate(OoqpVectorBase<T>& absminvec) const override {};
   void absmaxVecUpdate(OoqpVectorBase<T>& absmaxvec) const override {};
   void absmin( double& m) const override {};
   void absminNonZero(double& m, double zero_eps) const override {m=-1.0;};
   T stepbound(OoqpVectorBase<T> & v, double maxStep ) override {return maxStep;}
   double findBlocking(OoqpVectorBase<T> & wstep_vec,
			      OoqpVectorBase<T> & u_vec,
			      OoqpVectorBase<T> & ustep_vec,
			      double maxStep,
			      double *w_elt,
			      double *wstep_elt,
			      double *u_elt,
			      double *ustep_elt,
			      int& first_or_second)override {return maxStep;}

   void findBlocking_pd(const OoqpVectorBase<T>& wstep_vec,
               const OoqpVectorBase<T>& u_vec,
               const OoqpVectorBase<T>& ustep_vec,
               double& maxStepPri, double& maxStepDual,
               double& w_elt_p, double& wstep_elt_p, double& u_elt_p, double& ustep_elt_p,
               double& w_elt_d, double& wstep_elt_d, double& u_elt_d, double& ustep_elt_d,
               bool& primalBlocking, bool& dualBlocking) const override {};

   void componentMult( OoqpVectorBase<T>& v )override {};
   void componentDiv ( OoqpVectorBase<T>& v )override {};
   void scalarMult( double num)override {};
   void writeToStream(std::ostream& out) const override {};
   void writeToStreamAll(std::ostream& out) const override {};
   void writeToStreamAllChild( std::stringstream& sout ) const override {};
   void writefToStream( std::ostream& out,
			       const char format[] ) const override {};
   void writeMPSformatOnlyRhs(std::ostream& out, const std::string rowName, OoqpVectorBase<T>* irhs) const override {};
   void writeMPSformatRhs(std::ostream& out, int rowType, OoqpVectorBase<T>* irhs) const override {};
   void writeMPSformatBounds(std::ostream& out, OoqpVectorBase<T>* ix, bool upperBound) const override {};
   void writeMPSformatBoundsWithVar(std::ostream& out, const std::string varStub, OoqpVectorBase<T>* ix, bool upperBound) const override {};

   void scale( double alpha )override {};

  /** this += alpha * x */
   void axpy  ( double alpha, OoqpVectorBase<T>& x )override {};
  /** this += alpha * x * z */
   void axzpy ( double alpha, OoqpVectorBase<T>& x, OoqpVectorBase<T>& z )override {};
  /** this += alpha * x / z */
   void axdzpy( double alpha, OoqpVectorBase<T>& x, OoqpVectorBase<T>& z )override {};

   void addConstant( double c )override {};
   void gondzioProjection( double rmin, double rmax )override {};
   double dotProductWith( const OoqpVectorBase<T>& v ) const override {return 0.0;}
   double dotProductSelf(double scaleFactor = 1.0) const override {return 0.0;};

  /** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
   */
   double shiftedDotProductWith( double alpha, OoqpVectorBase<T>& mystep,
					OoqpVectorBase<T>& yvec,
					double beta,  OoqpVectorBase<T>& ystep )override {return 0.0;}
   void negate()override {};
   void invert()override {};
   void invertSave( double zeroReplacementVal = 0.0 )override {};
   void applySqrt()override {};
   void roundToPow2()override {};
   int allPositive()override {return 1;}

   int matchesNonZeroPattern( OoqpVectorBase<T>& select )override {return 1;}
   void selectNonZeros( OoqpVectorBase<T>& select )override {};
   long long numberOfNonzeros()override {return 0;}
   void addSomeConstants( double c, OoqpVectorBase<T>& select )override {};
   void writefSomeToStream( std::ostream& out,
				   const char format[],
				   OoqpVectorBase<T>& select ) const override {};
   void axdzpy( double alpha, OoqpVectorBase<T>& x,
		       OoqpVectorBase<T>& z, OoqpVectorBase<T>& select )override {};

   int somePositive( OoqpVectorBase<T>& select )override {return 1;}
   void divideSome( OoqpVectorBase<T>& div, OoqpVectorBase<T>& select )override {};
   void copyIntoArray( double v[] ) const override {};
   void copyFromArray( double v[] )override {};
   void copyFromArray( char v[] )override {};

   void removeEntries( const OoqpVectorBase<T>& select ) override {};
   void permuteVec0(const std::vector<unsigned int>& permvec) override {};
   void permuteLinkingEntries(const std::vector<unsigned int>& permvec) override {};
   std::vector<double> gatherStochVector() const override {return std::vector<double>(0);};

  int getSize() override { return 0; };

  bool isRootNodeInSync() const override { return true; };
};

using StochDummyVector = StochDummyVectorBase<double>;

#endif
