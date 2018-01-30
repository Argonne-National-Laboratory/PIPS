#ifndef STOCHVECTOR_H
#define STOCHVECTOR_H

#include "StochVectorHandle.h"
#include "OoqpVector.h"
#include "SimpleVector.h"
#include "mpi.h"

#include <vector>

class StochTree;

class StochVector : public OoqpVector {
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

  virtual int isKindOf( int kind );
  virtual void setToZero();
  virtual void setToConstant( double c );
  virtual void randomize( double alpha, double beta, double *ix );
  virtual void copyFrom( OoqpVector& v );
  virtual double twonorm();
  virtual double infnorm();
  virtual double onenorm();
  virtual void min( double& m, int& index );
  virtual void max( double& m, int& index );
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
  virtual void findBlocking_pd(OoqpVector & wstep_vec,
  			      OoqpVector & u_vec,
  			      OoqpVector & ustep_vec,
  			      double maxStepPri, double maxStepDual,
  			      double *w_elt_p,
  			      double *wstep_elt_p,
  			      double *u_elt_p,
  			      double *ustep_elt_p,
  				  double *w_elt_d, double *wstep_elt_d, double *u_elt_d, double *ustep_elt_d,
  				  double& stepPrimal, double& stepDual,
				  bool& primalBlocking, bool& dualBlocking);

  virtual void componentMult( OoqpVector& v );
  virtual void componentDiv ( OoqpVector& v );
  virtual void scalarMult( double num);
  virtual void writeToStream(ostream& out) const;
  virtual void writefToStream( ostream& out,
			       const char format[] ) const;

  virtual void scale( double alpha );

  /** this += alpha * x */
  virtual void axpy  ( double alpha, OoqpVector& x );
  /** this += alpha * x * z */
  virtual void axzpy ( double alpha, OoqpVector& x, OoqpVector& z );
  /** this += alpha * x / z */
  virtual void axdzpy( double alpha, OoqpVector& x, OoqpVector& z );

  virtual void addConstant( double c );
  virtual void gondzioProjection( double rmin, double rmax );
  virtual double dotProductWith( OoqpVector& v );
  /** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
   */
  virtual double shiftedDotProductWith( double alpha, OoqpVector& mystep,
					OoqpVector& yvec,
					double beta,  OoqpVector& ystep );
  virtual void negate();
  virtual void invert();
  virtual void invertSave( double zeroReplacementVal = 0.0 );
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

  int getSize() { return n; };
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

  virtual int isKindOf( int kind ){return kind==kStochDummy;}
  virtual void setToZero(){};
  virtual void setToConstant( double c ){};
  virtual void randomize( double alpha, double beta, double *ix ){};
  virtual void copyFrom( OoqpVector& v ){};
  virtual double twonorm(){return 0.0;}
  virtual double infnorm(){return 0.0;}
  virtual double onenorm(){return 0.0;}
  virtual void min( double& m, int& index ){};
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
  virtual void findBlocking_pd(OoqpVector & wstep_vec,
    			      OoqpVector & u_vec,
    			      OoqpVector & ustep_vec,
    			      double maxStepPri, double maxStepDual,
    			      double *w_elt_p,
    			      double *wstep_elt_p,
    			      double *u_elt_p,
    			      double *ustep_elt_p,
    				  double *w_elt_d, double *wstep_elt_d, double *u_elt_d, double *ustep_elt_d,
    				  double& stepPrimal, double& stepDual,
					  bool& primalBlocking, bool& dualBlocking){};

  virtual void componentMult( OoqpVector& v ){};
  virtual void componentDiv ( OoqpVector& v ){};
  virtual void scalarMult( double num){};
  virtual void writeToStream(ostream& out) const{};
  virtual void writefToStream( ostream& out,
			       const char format[] ) const{};

  virtual void scale( double alpha ){};

  /** this += alpha * x */
  virtual void axpy  ( double alpha, OoqpVector& x ){};
  /** this += alpha * x * z */
  virtual void axzpy ( double alpha, OoqpVector& x, OoqpVector& z ){};
  /** this += alpha * x / z */
  virtual void axdzpy( double alpha, OoqpVector& x, OoqpVector& z ){};

  virtual void addConstant( double c ){};
  virtual void gondzioProjection( double rmin, double rmax ){};
  virtual double dotProductWith( OoqpVector& v ){return 0.0;}
  /** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
   */
  virtual double shiftedDotProductWith( double alpha, OoqpVector& mystep,
					OoqpVector& yvec,
					double beta,  OoqpVector& ystep ){return 0.0;}
  virtual void negate(){};
  virtual void invert(){};
  virtual void invertSave( double zeroReplacementVal = 0.0 ){};
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

  int getSize() { return 0; };
};


#endif

