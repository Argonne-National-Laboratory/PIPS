#ifndef STOCHVECTOR_H
#define STOCHVECTOR_H

#include "StochVectorHandle.h"
#include "OoqpVector.h"
#include "SimpleVector.h"
#include "mpi.h"

#include <vector>

class StochTree;

class StochVector : public OoqpVector {
private:

  virtual void writeToStreamAllChild( stringstream& sout ) const;

protected:

public:
  StochVector( int n, MPI_Comm mpiComm, int isDistributed=-1);
  StochVector( int n, int nl, MPI_Comm mpiComm, int isDistributed);
  virtual ~StochVector();

  void AddChild(StochVector* child);
  void AddChild(OoqpVector* child);

  /** The data for this node. */
  OoqpVector*               vec;

  /** The linking constraint data for this node. */
  OoqpVector*               vecl;

  /** Children of this node */
  std::vector<StochVector*> children;

  /** Link to the parent of this node. Needed when we multiply a matrix
      with this vector
  */
  StochVector*              parent;


  /* MPI communicator */
  MPI_Comm mpiComm;
  /* flag used to indicate if the children are distributed or not. */
  int iAmDistrib;

  /** Creates and returns a vector of the type used to store data in this node,
      i.e., same type as 'vec'.
      NO data is copied, for this use one of the 'copy...' functions.
  */
  virtual OoqpVector* dataClone() const;
  virtual OoqpVector* dataCloneLinkCons() const;
  virtual StochVector* clone() const;
  /* copy vector entries as well */
  virtual StochVector* cloneFull() const;

  virtual void jointCopyFrom(StochVector& v1, StochVector& v2, StochVector& v3);
  virtual void jointCopyFromLinkCons(StochVector& vx, StochVector& vy, StochVector& vz);
  virtual void jointCopyTo(StochVector& v1, StochVector& v2, StochVector& v3);
  virtual void jointCopyToLinkCons(StochVector& vx, StochVector& vy, StochVector& vz);

  virtual int isKindOf( int kind ) const;
  virtual void setToZero();
  virtual void setToConstant( double c );
  virtual bool isZero() const;

  virtual void randomize( double alpha, double beta, double *ix );
  virtual void copyFrom( OoqpVector& v );
  virtual void copyFromAbs(const OoqpVector& v);
  virtual double twonorm() const;
  virtual double infnorm() const;
  virtual double onenorm() const;
  virtual void min( double& m, int& index ) const;
  virtual void max( double& m, int& index ) const;
  virtual void absminVecUpdate(OoqpVector& absminvec) const;
  virtual void absmaxVecUpdate(OoqpVector& absmaxvec) const;
  virtual void absmin( double& m ) const;
  virtual void absminNonZero(double& m, double zero_eps) const;
  virtual double stepbound(OoqpVector & v, double maxStep );
  virtual double findBlocking(OoqpVector & wstep_vec,
			      OoqpVector & u_vec,
			      OoqpVector & ustep_vec,
			      double maxStep,
			      double *w_elt,
			      double *wstep_elt,
			      double *u_elt,
			      double *ustep_elt,
			      int& first_or_second);
  virtual void findBlocking_pd(const OoqpVector& wstep_vec,
               const OoqpVector& u_vec,
  			      const OoqpVector& ustep_vec,
  			      double& maxStepPri, double& maxStepDual,
  			      double& w_elt_p, double& wstep_elt_p, double& u_elt_p, double& ustep_elt_p,
  				   double& w_elt_d, double& wstep_elt_d, double& u_elt_d, double& ustep_elt_d,
				   bool& primalBlocking, bool& dualBlocking) const;

  virtual void componentMult( OoqpVector& v );
  virtual void componentDiv ( OoqpVector& v );
  virtual bool componentEqual( const OoqpVector& v, double tol) const;
  virtual void scalarMult( double num);
  virtual void writeToStream(ostream& out) const;
  virtual void writeToStreamAll(ostream& out) const;
  virtual void writefToStream( ostream& out,
			       const char format[] ) const;
  virtual void writeMPSformatOnlyRhs(ostream& out, string rowName, OoqpVector* irhs) const {};
  virtual void writeMPSformatRhs(ostream& out, int rowType, OoqpVector* irhs) const;
  virtual void writeMPSformatBounds(ostream& out, OoqpVector* ix, bool upperBound) const;
  virtual void writeMPSformatBoundsWithVar(ostream& out, string varStub, OoqpVector* ix, bool upperBound) const {};

  virtual void scale( double alpha );

  /** this += alpha * x */
  virtual void axpy  ( double alpha, OoqpVector& x );
  /** this += alpha * x * z */
  virtual void axzpy ( double alpha, OoqpVector& x, OoqpVector& z );
  /** this += alpha * x / z */
  virtual void axdzpy( double alpha, OoqpVector& x, OoqpVector& z );

  virtual void addConstant( double c );
  virtual void gondzioProjection( double rmin, double rmax );
  virtual double dotProductWith( const OoqpVector& v ) const;
  virtual double dotProductSelf(double scaleFactor) const;

  /** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
   */
  virtual double shiftedDotProductWith( double alpha, OoqpVector& mystep,
					OoqpVector& yvec,
					double beta,  OoqpVector& ystep );
  virtual void negate();
  virtual void invert();
  virtual void invertSave( double zeroReplacementVal = 0.0 );
  virtual void applySqrt();
  virtual void roundToPow2();

  virtual int allPositive();

  virtual int matchesNonZeroPattern( OoqpVector& select );
  virtual void selectNonZeros( OoqpVector& select );
  virtual long long numberOfNonzeros();
  virtual void addSomeConstants( double c, OoqpVector& select );
  virtual void writefSomeToStream( ostream& out,
				   const char format[],
				   OoqpVector& select ) const;
  virtual void axdzpy( double alpha, OoqpVector& x,
		       OoqpVector& z, OoqpVector& select );

  virtual int somePositive( OoqpVector& select );
  virtual void divideSome( OoqpVector& div, OoqpVector& select );
  virtual void copyIntoArray( double v[] ) const;
  virtual void copyFromArray( double v[] );
  virtual void copyFromArray( char v[] );
  virtual void permuteVec0Entries(const std::vector<unsigned int>& permvec);
  virtual void permuteLinkingEntries(const std::vector<unsigned int>& permvec);
  virtual std::vector<double> gatherStochVector() const;

  /** remove entries i for which select[i] == 0 */
  virtual void removeEntries( const OoqpVector& select );

  int getSize() { return n; };

  virtual bool isRootNodeInSync() const;

};

/** DUMMY VERSION
 *
 */
class StochDummyVector : public StochVector {
protected:

public:
  StochDummyVector(  )
    : StochVector(0, MPI_COMM_NULL) {};

  virtual ~StochDummyVector(){};

  void AddChild(StochVector* child){};
  void AddChild(OoqpVector* child){};

  /** Creates and returns a vector of the type used to store data in this node,
      i.e., same type as 'vec'.
      NO data is copied, for this use one of the 'copy...' functions.
  */
  virtual OoqpVector* dataClone() const { return new SimpleVector(0);}
  virtual StochVector* clone() const { return new StochDummyVector();}
  virtual StochVector* cloneFull() const { return new StochDummyVector();}

  virtual void jointCopyFrom(StochVector& v1, StochVector& v2, StochVector& v3){};
  virtual void jointCopyTo(StochVector& v1, StochVector& v2, StochVector& v3){};

  virtual int isKindOf( int kind ) const {return kind == kStochDummy;}
  virtual void setToZero(){};
  virtual bool isZero() const { return true; };

  virtual void setToConstant( double c ){};
  virtual void randomize( double alpha, double beta, double *ix ){};
  virtual void copyFrom( OoqpVector& v ){};
  virtual void copyFromAbs(const OoqpVector& v) {};
  virtual double twonorm() const {return 0.0;}
  virtual double infnorm() const {return 0.0;}
  virtual double onenorm() const {return 0.0;}
  virtual void min( double& m, int& index ) const {};
  virtual void max( double& m, int& index ) const {};
  virtual void absminVecUpdate(OoqpVector& absminvec) const {};
  virtual void absmaxVecUpdate(OoqpVector& absmaxvec) const {};
  virtual void absmin( double& m) const {};
  virtual void absminNonZero(double& m, double zero_eps) const {m=-1.0;};
  virtual double stepbound(OoqpVector & v, double maxStep ){return maxStep;}
  virtual double findBlocking(OoqpVector & wstep_vec,
			      OoqpVector & u_vec,
			      OoqpVector & ustep_vec,
			      double maxStep,
			      double *w_elt,
			      double *wstep_elt,
			      double *u_elt,
			      double *ustep_elt,
			      int& first_or_second){return maxStep;}

  virtual void findBlocking_pd(const OoqpVector& wstep_vec,
               const OoqpVector& u_vec,
               const OoqpVector& ustep_vec,
               double& maxStepPri, double& maxStepDual,
               double& w_elt_p, double& wstep_elt_p, double& u_elt_p, double& ustep_elt_p,
               double& w_elt_d, double& wstep_elt_d, double& u_elt_d, double& ustep_elt_d,
               bool& primalBlocking, bool& dualBlocking) const {};

  virtual void componentMult( OoqpVector& v ){};
  virtual void componentDiv ( OoqpVector& v ){};
  virtual bool componentEqual( const OoqpVector& v, double tol) const { if(!v.isKindOf(kStochDummy)) std::cout << "lol one should never end up here"
        << std::endl; return v.isKindOf(kStochDummy); };
  virtual void scalarMult( double num){};
  virtual void writeToStream(ostream& out) const{};
  virtual void writeToStreamAll(ostream& out) const{};
  virtual void writeToStreamAllChild( stringstream& sout ) const{};
  virtual void writefToStream( ostream& out,
			       const char format[] ) const{};
  virtual void writeMPSformatOnlyRhs(ostream& out, string rowName, OoqpVector* irhs) const{};
  virtual void writeMPSformatRhs(ostream& out, int rowType, OoqpVector* irhs) const{};
  virtual void writeMPSformatBounds(ostream& out, OoqpVector* ix, bool upperBound) const {};
  virtual void writeMPSformatBoundsWithVar(ostream& out, string varStub, OoqpVector* ix, bool upperBound) const {};

  virtual void scale( double alpha ){};

  /** this += alpha * x */
  virtual void axpy  ( double alpha, OoqpVector& x ){};
  /** this += alpha * x * z */
  virtual void axzpy ( double alpha, OoqpVector& x, OoqpVector& z ){};
  /** this += alpha * x / z */
  virtual void axdzpy( double alpha, OoqpVector& x, OoqpVector& z ){};

  virtual void addConstant( double c ){};
  virtual void gondzioProjection( double rmin, double rmax ){};
  virtual double dotProductWith( const OoqpVector& v ) const {return 0.0;}
  virtual double dotProductSelf(double scaleFactor = 1.0) const {return 0.0;};

  /** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
   */
  virtual double shiftedDotProductWith( double alpha, OoqpVector& mystep,
					OoqpVector& yvec,
					double beta,  OoqpVector& ystep ){return 0.0;}
  virtual void negate(){};
  virtual void invert(){};
  virtual void invertSave( double zeroReplacementVal = 0.0 ){};
  virtual void applySqrt(){};
  virtual void roundToPow2(){};
  virtual int allPositive(){return 1;}

  virtual int matchesNonZeroPattern( OoqpVector& select ){return 1;}
  virtual void selectNonZeros( OoqpVector& select ){};
  virtual long long numberOfNonzeros(){return 0;}
  virtual void addSomeConstants( double c, OoqpVector& select ){};
  virtual void writefSomeToStream( ostream& out,
				   const char format[],
				   OoqpVector& select ) const{};
  virtual void axdzpy( double alpha, OoqpVector& x,
		       OoqpVector& z, OoqpVector& select ){};

  virtual int somePositive( OoqpVector& select ){return 1;}
  virtual void divideSome( OoqpVector& div, OoqpVector& select ){};
  virtual void copyIntoArray( double v[] ) const{};
  virtual void copyFromArray( double v[] ){};
  virtual void copyFromArray( char v[] ){};

  virtual void removeEntries( const OoqpVector& select ) {};
  virtual void permuteVec0(const std::vector<unsigned int>& permvec) {};
  virtual void permuteLinkingEntries(const std::vector<unsigned int>& permvec) {};
  virtual std::vector<double> gatherStochVector() const {return std::vector<double>(0);};

  int getSize() { return 0; };

  virtual bool isRootNodeInSync() const { return true; };
};


#endif
