/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef STOCHVECTOR_H
#define STOCHVECTOR_H

#include "mpi.h"
#include "StochVectorHandle.h"
#include "OoqpVector.h"
#include "SimpleVector.h"
#include <unordered_map>

#include <vector>

class StochTree;
class StochDummyVector;

template <class pTree, class sTree> class VectorCompressedDummy {

public:

 size_t size() const {
   return size_;
 }
 void clear() {
   size_ = 0;
   children.clear();
 }

 pTree& operator[](std::size_t idx) {
   if (children.find(idx) == children.end()) {
     return sTree::dummy;
   } else {
     return children[idx];
   }
 }

 const pTree& operator[](std::size_t idx) const {
   if (children.find(idx) == children.end()) {
     return sTree::dummy;
   } else {
     return children[idx];
   }
 }
 void push_back(const pTree child) {
   if (child != sTree::dummy) {
     children[size_] = child;
   }
   size_++;
 }

private:
  std::unordered_map<int, pTree> children;
  size_t size_ = 0;
  size_t tablesize_ = 0;
};


class StochVector : public OoqpVector {
protected:
	bool firstDoNumOfNonZero;
	int nnzNonZeros;
	
public:
  StochVector( int n, MPI_Comm mpiComm, int isDistributed=-1);
  StochVector( int n, int const treeIDX_in, MPI_Comm mpiComm, int isDistributed=-1);
  
  virtual ~StochVector();

  void AddChild(StochVector* child);
  void AddChild(OoqpVector* child);

  // corresponding tree index
  int treeIDX;

  /** The data for this node. */
  OoqpVector*               vec;

  /** Children of this node */
  VectorCompressedDummy<StochVector*, StochVector> children;

  /** Link to the parent of this node. Needed when we multiply a matrix 
      with this vector
  */
  StochVector*              parent;


  /* MPI communicator */
  MPI_Comm mpiComm;
  /* flag used to indicate if the children are distributed or not. */
  int iAmDistrib;
  static StochVector *dummy;

  /** Creates and returns a vector of the type used to store data in this node,
      i.e., same type as 'vec'.
      NO data is copied, for this use one of the 'copy...' functions.
  */
  virtual OoqpVector* dataClone() const;
  virtual StochVector* clone();

  virtual void jointCopyFrom(StochVector& v1, StochVector& v2, StochVector& v3);
  virtual void jointCopyTo(StochVector& v1, StochVector& v2, StochVector& v3);
  virtual void jointCopyFromXSYZ(StochVector& v1, StochVector& v2, StochVector& v3, StochVector& v4);
  virtual void jointCopyToXSYZ(StochVector& v1, StochVector& v2, StochVector& v3, StochVector& v4); 

  virtual int isKindOf( int kind );
  virtual void setToZero();
  virtual void setToConstant( double c );
  virtual void randomize( double alpha, double beta, double *ix );
  virtual void copyFrom( OoqpVector& v );
  virtual double twonorm();
  virtual double infnorm();
  virtual double onenorm();

  // FIXME_ME: index is wrong! Do not use it in parallel
  virtual void min( double& m, int& index );
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

  virtual void componentMult( OoqpVector& v );
  virtual void componentDiv ( OoqpVector& v );
  virtual void scalarMult( double num);
  virtual void writeToStream(std::ostream& out) const;
  virtual void writefToStream( std::ostream& out,
			       const char format[] );

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
  virtual int allPositive();

  virtual int matchesNonZeroPattern( OoqpVector& select );
  virtual void selectNonZeros( OoqpVector& select );
  virtual long long numberOfNonzeros();
  virtual void addSomeConstants( double c, OoqpVector& select );
  virtual void writefSomeToStream( std::ostream& out,
				   const char format[],
				   OoqpVector& select );
  virtual void axdzpy( double alpha, OoqpVector& x,
		       OoqpVector& z, OoqpVector& select );

  virtual int somePositive( OoqpVector& select );
  virtual void divideSome( OoqpVector& div, OoqpVector& select );
  virtual void copyIntoArray( double v[] ) ;
  virtual void copyFromArray( double v[] );
  virtual void copyFromArray( char v[] );

  int getSize() { return n; };




 virtual void print();

  virtual void findBlockingPD(OoqpVector & wstep_vec, 
				      OoqpVector & u_vec, 
				      OoqpVector & ustep_vec, 
				      double *w_elt, 
				      double *wstep_elt,
				      double *u_elt, 
				      double *ustep_elt,
				      int& first_or_second, double * alphaPri, double * alphaDual);

  /* following routines are added by Naiyuan 2013 */
  // log function
  virtual double sumLog(OoqpVector* select);
  // sum of the power of elts, can be used to evaluate norm
  virtual double sumPowElt(const double pow_in);
  // sum of elts
  virtual double sumElt();


  virtual void MinComponentOrConstant( OoqpVector* vec_in, double minVal  );

  virtual void correctLargeVal( const double testVal, const double corVal, const int absFlag);  
  
  virtual void MinComponentBetween( OoqpVector* vec_in,OoqpVector *select_in=NULL);

  virtual void setToConstantFromTo( double c, int start, int end );  

  virtual void SetComponentFromMaxXorY( OoqpVector* x_in, OoqpVector *y_in,OoqpVector *select_in=NULL);

  virtual void SetComponentFromMaxXorY( OoqpVector* x_in, OoqpVector *y_in, 
			  int Start, int End,int xStart, int xEnd, int yStart, int yEnd);

  virtual void SetComponentFromMaxXorConstant( OoqpVector* x_in, const double y_in);

  virtual void SetComponentFromMaxXorConstant( OoqpVector* x_in, const double y_in, 
			int Start, int End,int xStart, int xEnd);


  virtual void copyFromFromTo( OoqpVector* vec_in, int VStart, int VEnd, int VinStart, int VinEnd);

  virtual void copyFromArrayFromTo( double *w, int VStart, int VEnd, int Wstart, int Wend);
  
  virtual void copyIntoArrayFromTo( double *w, int VStart, int VEnd, int Wstart, int Wend);

  virtual void absVal(OoqpVector *vec_in);

  virtual void copyIntoArrayWithIndex_AggVarCon( double *ResultArray, const int *vecmap, const int length, bool isVar=true);
  
};

/** DUMMY VERSION 
 *
 */
class StochDummyVector : public StochVector {
protected:
  bool firstDoNumOfNonZero;
  int nnzNonZeros;

public:
  StochDummyVector(  )
    : StochVector(0, MPI_COMM_NULL) {};

  virtual ~StochDummyVector(){};

  static StochDummyVector *dummy;

  void AddChild(StochVector* child){};
  void AddChild(OoqpVector* child){};

  /** Creates and returns a vector of the type used to store data in this node,
      i.e., same type as 'vec'.
      NO data is copied, for this use one of the 'copy...' functions.
  */
  virtual OoqpVector* dataClone() const { return new SimpleVector(0);}
  virtual StochVector* clone() const { return new StochDummyVector();}

  virtual void jointCopyFrom(StochVector& v1, StochVector& v2, StochVector& v3){};
  virtual void jointCopyTo(StochVector& v1, StochVector& v2, StochVector& v3){};

  virtual void jointCopyFromXSYZ(StochVector& v1, StochVector& v2, StochVector& v3, StochVector& v4){};
  virtual void jointCopyToXSYZ(StochVector& v1, StochVector& v2, StochVector& v3, StochVector& v4){};

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

  virtual void componentMult( OoqpVector& v ){};
  virtual void componentDiv ( OoqpVector& v ){};
  virtual void scalarMult( double num){};
  virtual void writeToStream(std::ostream& out) const{};
  virtual void writefToStream( std::ostream& out,
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
  virtual int allPositive(){return 1;}

  virtual int matchesNonZeroPattern( OoqpVector& select ){return 1;}
  virtual void selectNonZeros( OoqpVector& select ){};
  virtual long long numberOfNonzeros(){return 0;}
  virtual void addSomeConstants( double c, OoqpVector& select ){};
  virtual void writefSomeToStream( std::ostream& out,
				   const char format[],
				   OoqpVector& select ) const{};
  virtual void axdzpy( double alpha, OoqpVector& x,
		       OoqpVector& z, OoqpVector& select ){};

  virtual int somePositive( OoqpVector& select ){return 1;}
  virtual void divideSome( OoqpVector& div, OoqpVector& select ){};
  virtual void copyIntoArray( double v[] ) {};
  virtual void copyFromArray( double v[] ){};
  virtual void copyFromArray( char v[] ){};

  int getSize() { return 0; };



  /* following routines are added by Naiyuan 2013 */

 virtual void print(){};

  virtual void findBlockingPD(OoqpVector & wstep_vec, 
				      OoqpVector & u_vec, 
				      OoqpVector & ustep_vec, 
				      double *w_elt, 
				      double *wstep_elt,
				      double *u_elt, 
				      double *ustep_elt,
				      int& first_or_second, double * alphaPri, double * alphaDual){};


  // log function
  virtual double sumLog(OoqpVector& select){return 0.0;};
  // sum of the power of elts, can be used to evaluate norm
  virtual double sumPowElt(const double pow_in){return 0.0;};
  // sum of elts
  virtual double sumElt(){return 0.0;};


  virtual void MinComponentOrConstant( OoqpVector* vec_in, double minVal  ){};

  virtual void correctLargeVal( const double testVal, const double corVal, const int absFlag){};  
  
  virtual void MinComponentBetween( OoqpVector* vec_in,OoqpVector *select_in){};

  virtual void setToConstantFromTo( double c, int start, int end ){};  

  virtual void SetComponentFromMaxXorY( OoqpVector* x_in, OoqpVector *y_in,OoqpVector *select_in){};

  virtual void SetComponentFromMaxXorY( OoqpVector* x_in, OoqpVector *y_in, 
			  int Start, int End,int xStart, int xEnd, int yStart, int yEnd){};

  virtual void SetComponentFromMaxXorConstant( OoqpVector* x_in, const double y_in){};

  virtual void SetComponentFromMaxXorConstant( OoqpVector* x_in, const double y_in, 
			int Start, int End,int xStart, int xEnd){};


  virtual void copyFromFromTo( OoqpVector* vec_in, int VStart, int VEnd, int VinStart, int VinEnd){};

  virtual void copyFromArrayFromTo( double *w, int VStart, int VEnd, int Wstart, int Wend){};

  virtual void copyIntoArrayFromTo( double *w, int VStart, int VEnd, int Wstart, int Wend){};

  virtual void absVal(OoqpVector *vec_in){};
  
  virtual void copyIntoArrayWithIndex_AggVarCon( double *ResultArray, const int *vecmap, const int length, bool isVar=true){};
};


#endif

