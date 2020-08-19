/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef QPGENDATA
#define QPGENDATA

#include "Data.h"
#include "OoqpVectorHandle.h"
#include "OoqpVector.h"
#include "DoubleMatrixHandle.h"

class MpsReader;
class LinearAlgebraPackage;
class QpGenVars;

#ifdef TESTING
class QpGenDataTester;
#endif

/**
 * Data for the general QP formulation.
 *
 * @ingroup QpGen
 */

class QpGenData : public Data {
#ifdef TESTING
  friend QpGenDataTester;
#endif
private:

  /** as part of setting up a random test problem, generate a random
   *  set of upper, lower, and two-sided bounds */
  void randomlyChooseBoundedVariables( OoqpVector& x, OoqpVector& dualx,
				       OoqpVector& blx, OoqpVector& ixlow,
				       OoqpVector& bux, OoqpVector& ixupp,
				       double * ix,
				       double percentLowerOnly,
				       double percentUpperOnly,
				       double percentBound );
protected:
  LinearAlgebraPackage * la;
  QpGenData();
public:
  SymMatrixHandle Q;
  GenMatrixHandle A;
  GenMatrixHandle C;
  OoqpVectorHandle    g; // objective
  OoqpVectorHandle    bA; // rhs equality
  OoqpVectorHandle    bux; // upper bounds x
  OoqpVectorHandle    ixupp; // index for upper bounds
  OoqpVectorHandle    blx; // lower bounds x
  OoqpVectorHandle    ixlow; // index for lower bounds
  OoqpVectorHandle    bu; // upper bounds C
  OoqpVectorHandle    icupp; // index upper bounds
  OoqpVectorHandle    bl; // lower bounds C
  OoqpVectorHandle    iclow; // index lower bounds
  OoqpVectorHandle    sc; // scale

  long long nx, my, mz;

  /** constructor that makes data objects of the specified dimensions */
  QpGenData(LinearAlgebraPackage * la,
	    long long nx_, long long my_, long long mz_,
	    long long nnzQ, long long nnzA, long long nnzC);

  /** constructor that sets up pointers to the data objects that are
      passed as arguments */
  QpGenData( LinearAlgebraPackage * la,
	     OoqpVector * c, SymMatrix * Q,
	     OoqpVector * xlow, OoqpVector * ixlow,
	     OoqpVector * xupp, OoqpVector * ixupp,
	     GenMatrix * A, OoqpVector * bA,
	     GenMatrix * C,
	     OoqpVector * clow, OoqpVector * iclow,
	     OoqpVector * cupp, OoqpVector * ciupp );

  /** insert the Hessian Q into the matrix M for the fundamental
      linear system, where M is stored as a GenMatrix */
  virtual void putQIntoAt( GenMatrix& M, int row, int col );

  /** insert the constraint matrix A into the matrix M for the
      fundamental linear system, where M is stored as a GenMatrix */
  virtual void putAIntoAt( GenMatrix& M, int row, int col );

  /** insert the constraint matrix C into the matrix M for the
      fundamental linear system, where M is stored as a GenMatrix */
  virtual void putCIntoAt( GenMatrix& M, int row, int col );

  /** insert the Hessian Q into the matrix M for the fundamental
      linear system, where M is stored as a SymMatrix */
  virtual void putQIntoAt( SymMatrix& M, int row, int col );

  /** insert the constraint matrix A into the matrix M for the
      fundamental linear system, where M is stored as a SymMatrix */
  virtual void putAIntoAt( SymMatrix& M, int row, int col );

  /** insert the constraint matrix C into the matrix M for the
      fundamental linear system, where M is stored as a SymMatrix */
  virtual void putCIntoAt( SymMatrix& M, int row, int col );

  /** y = beta * y + alpha * Q * x */
  virtual void Qmult( double beta,  OoqpVector& y,
		      double alpha, const OoqpVector& x ) const;

  /** y = beta * y + alpha * A * x */
  virtual void Amult( double beta,  OoqpVector& y,
		      double alpha, const OoqpVector& x) const;

  /** y = beta * y + alpha * C * x   */
  virtual void Cmult( double beta,  OoqpVector& y,
		      double alpha, const OoqpVector& x ) const;

  /** y = beta * y + alpha * A\T * x */
  virtual void ATransmult( double beta,  OoqpVector& y,
			   double alpha, const OoqpVector& x ) const;

  /** y = beta * y + alpha * C\T * x */
  virtual void CTransmult( double beta,  OoqpVector& y,
			   double alpha, const OoqpVector& x ) const;

  //  virtual void addSymProdCRowToAt(double alpha, int i,
  //				  SymMatrix& M, int rowcol );

  virtual void getg(  OoqpVector& cout ) const;
  virtual void getbA( OoqpVector& bout ) const;

  /** extract the diagonal of Q and put it in the OoqpVector dQ */
  virtual void getDiagonalOfQ( OoqpVector& dQ );

  virtual OoqpVector&  xupperBound() { return *bux; };
  virtual OoqpVector& ixupperBound() { return *ixupp; };
  virtual OoqpVector&  xlowerBound() { return *blx; };
  virtual OoqpVector& ixlowerBound() { return *ixlow; };
  virtual OoqpVector&  supperBound() { return *bu; };
  virtual OoqpVector& isupperBound() { return *icupp; };
  virtual OoqpVector&  slowerBound() { return *bl; };
  virtual OoqpVector& islowerBound() { return *iclow; };
  virtual OoqpVector& scale() { return *sc; };

  virtual void createScaleFromQ();
  virtual void scaleQ();
  virtual void scaleA();
  virtual void scaleC();
  virtual void scaleg();
  virtual void scalexupp();
  virtual void scalexlow();

  virtual void flipg();
  virtual void flipQ();

  double datanorm() const override;
  virtual void datainput() {};
  virtual void datainput( MpsReader * reader, int& iErr );
  /** Create a random problem
   *  @param (x,y,z,s) the solution to the random problem
   */
  virtual void datarandom( OoqpVector  & x, OoqpVector  & y,
			    OoqpVector & z, OoqpVector & s );
  virtual void print();

  virtual double objectiveValue( const QpGenVars * vars ) const;

  virtual ~QpGenData();
};

#endif
