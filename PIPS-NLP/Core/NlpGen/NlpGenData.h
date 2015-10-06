/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef NLPGENDATA
#define NLPGENDATA

#include "Data.h"
#include "OoqpVectorHandle.h"
#include "OoqpVector.h"
#include "DoubleMatrixHandle.h"
#include <map>


class MpsReader;
class LinearAlgebraPackage;


class NlpGenVars;
class NlpInfo;
class NlpGenResiduals;
class Residuals;


#ifdef TESTING
class NlpGenDataTester;
#endif


class NlpGenData : public Data {
#ifdef TESTING
	  friend NlpGenDataTester;
#endif

 protected:
  //
  // The input class (base class)
  // Derived classes will implement functionality for dealing with AMPL, SML, or PYOMO
  // When coding a new problem by hand, a class that inherits this will be needed to be implemented

  LinearAlgebraPackage * la;
  NlpGenData();
  
 public:
  SymMatrixHandle  H; //Hessian
  GenMatrixHandle  Jeq;
  GenMatrixHandle  Jineq;

  
  OoqpVectorHandle grad; 

  OoqpVectorHandle trialBarrGrad_x;    //barrier gradient for trial step, for x
  OoqpVectorHandle trialBarrGrad_s;    //barrier gradient for trial step, for s


  OoqpVectorHandle dampind_xL_v; 	//damping index vector for x lb
  OoqpVectorHandle dampind_xU_w; 	//damping index vector for x ub
  OoqpVectorHandle dampind_sL_t; 	//damping index vector for ineq cons lb
  OoqpVectorHandle dampind_sU_u; 	//damping index vector for ineq cons ub
  
  //
  // Assuming that the NLP problem is formulated as QPs are in OOQP
  //

  // c(x) = b
  OoqpVectorHandle    bA;
  OoqpVectorHandle    CeqBody;
  OoqpVectorHandle    trialCeqBody;

  // bl <= c(x) <= bu 
  OoqpVectorHandle    CIneqBody;
  OoqpVectorHandle    trialCIneqBody;  
  OoqpVectorHandle    bu; 
  OoqpVectorHandle    icupp;
  OoqpVectorHandle    bl;
  OoqpVectorHandle    iclow;
  OoqpVectorHandle    sc;


  //blx <= x <=  bux
  OoqpVectorHandle    bux; //upper bounds
  OoqpVectorHandle    ixupp; //indexes  showing which variable are bounded
  OoqpVectorHandle    blx; // lower bounds 
  OoqpVectorHandle    ixlow;


  long long nx;
  long long nxupp, nxlow;
  long long my;
  long long mz;
  long long mcupp, mclow;

  double currMu;

  double PriObj;
  double BarrObj;
  double MeritObj;

  double kktTh;
  
  
  NlpInfo* inputNlp;

  long long nxOri,myOri,mzOri;
  long long nxUOri,nxLOri,nsUOri,nsLOri;

  double linsysRes, linsysRes_Full;

  // number of krylov iterations
  int KryIter;

  int* schurVarConID;
  int schurSize;


  int* var_Part_idx_in; 
  int* con_Part_idx_in; 

  /** constructor that makes data objects of the specified dimensions */
  NlpGenData(LinearAlgebraPackage * la_,
		     long long  nx_, long long my_, long long mz_,
		     long long nnzQ_, long long nnzA_, long long nnzC_,
		     long long nxL_in,long long nxU_in,long long nsL_in,long long nsU_in);


  NlpGenData(LinearAlgebraPackage * la_,
		     int nx_, int my_, int mz_,
		     int nnzQ_, int nnzA_, int nnzC_);



  /** constructor that sets up pointers to the data objects that are
      passed as arguments */
  NlpGenData( LinearAlgebraPackage * la_in,
		      OoqpVector * grad_in, SymMatrix * H_in,
		      OoqpVector * xlow_in, OoqpVector * ixlow_in, long long nxlow_,
		      OoqpVector * xupp_in, OoqpVector * ixupp_in, long long nxupp_,
		      GenMatrix  * A_in, OoqpVector * bA_in,
		      GenMatrix  * C_in,
		      OoqpVector * clow_in, OoqpVector * iclow_in, long long mclow_,
		      OoqpVector * cupp_in, OoqpVector * icupp_in, long long mcupp_,
		      OoqpVector * CeqBody_in, OoqpVector * CIneqBody_in,
		      OoqpVector * trialBarrGrad_x_in,OoqpVector * trialBarrGrad_s_in,
	     	  OoqpVector * trialCeqBody, OoqpVector *trialCIneqBody,
	     	  OoqpVector * dampind_xL_v,OoqpVector * dampind_xU_w,
	     	  OoqpVector * dampind_sL_t, OoqpVector *dampind_sU_u);


  /** constructor that sets up pointers to the data objects that are
      passed as arguments */
  NlpGenData( LinearAlgebraPackage * la_in,
		      OoqpVector * grad_in, SymMatrix * H_in,
		      OoqpVector * xlow_in, OoqpVector * ixlow_in,
		      OoqpVector * xupp_in, OoqpVector * ixupp_in,
		      GenMatrix  * A_in, OoqpVector * bA_in,
		      GenMatrix  * C_in,
		      OoqpVector * clow_in, OoqpVector * iclow_in,
		      OoqpVector * cupp_in, OoqpVector * icupp_in);


  void setDampingVarMap();

  //
  // NlpGenData functionality below from QpGenData 
  //
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
		      double alpha, OoqpVector& x );

  /** y = beta * y + alpha * A * x */
  virtual void Amult( double beta,  OoqpVector& y,
		      double alpha, OoqpVector& x);

  /** y = beta * y + alpha * C * x   */
  virtual void Cmult( double beta,  OoqpVector& y,
		      double alpha, OoqpVector& x );

  /** y = beta * y + alpha * A\T * x */
  virtual void ATransmult( double beta,  OoqpVector& y,
			   double alpha, OoqpVector& x );

  /** y = beta * y + alpha * C\T * x */
  virtual void CTransmult( double beta,  OoqpVector& y,
			   double alpha, OoqpVector& x );

  //  virtual void addSymProdCRowToAt(double alpha, int i, 
  //				  SymMatrix& M, int rowcol );

  virtual void getg(  OoqpVector& cout );
  


  virtual void getbA( OoqpVector& bout );
  virtual void getInEqCons( OoqpVector& InEqCon_out );

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
  virtual OoqpVector& scale(){ return *sc; };

  virtual void createScaleFromQ();
  virtual void scaleQ();
  virtual void scaleA();
  virtual void scaleC();
  virtual void scaleg();
  virtual void scalexupp();
  virtual void scalexlow();

  virtual void flipg();
  virtual void flipQ();

  virtual double datanorm();
  virtual void datainput() {};
  virtual void datainput( MpsReader * reader, int& iErr );
  /** Create a random problem 
   *  @param (x,y,z,s) the solution to the random problem
   */
//  virtual void datarandom( OoqpVector  & x, OoqpVector  & y,
//			    OoqpVector & z, OoqpVector & s );
  virtual void print();



  virtual ~NlpGenData();









   virtual double objectiveValue( Variables* vars );
   virtual double BarrObjValue( NlpGenVars * vars, double PriObj_in, double dampingFact=0 );
   virtual double evalMeritFunc( NlpGenVars * vars, double penalty, double BarrObj_in);
   virtual double evalMeritFunc( double priErr, double penalty, double BarrObj_in);

   virtual double evalMeritFunc( double BarrObj_in, Variables *iterate_in, Residuals *resid_in);
   virtual double evalScaledConstraintNorm( Variables *iterate_in, Residuals *resid_in, const int isTrialStep=0 );
   virtual double getConTimesD( 	Variables* vars_in, Variables* steps_in, Residuals *resid_in);   




   virtual void evalData(Variables* vars);






	virtual void evalConstraintBody(Variables* vars_in, const int IfTrialStep=0);

   

	virtual double getBarrGradTimesD(Variables* vars_in, Variables* steps_in, 
						const int IfTrialStep=0, const double dampingFact=0);

	virtual void moveBounds(OoqpVector *priWrk_X, OoqpVector *priWrk_S, const double tol);

	virtual void getInitX(OoqpVector *initVecX);



virtual void setQIntoAt( SymMatrix& M, int row, int col, bool firstCall, std::map<int,int> &ValIdxMap );
virtual void setAIntoAt( SymMatrix& M, int row, int col, bool firstCall, std::map<int,int> &ValIdxMap );
virtual void setCIntoAt( SymMatrix& M, int row, int col, bool firstCall, std::map<int,int> &ValIdxMap );



};
#endif
