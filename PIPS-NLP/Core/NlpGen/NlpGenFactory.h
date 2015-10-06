/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef SPSEQNLPGENFACTORY
#define SPSEQNLPGENFACTORY

#include "ProblemFormulation.h"

class Data;
class Residuals;
class LinearSystem;
class Variables;
class LinearAlgebraPackage;
class OoqpVector;

/**
 * @defgroup NLPGen
 * 
 *  <pre>
 *  minimize    f(x)         ; 
 *  subject to                      C_e(x)  = b    ;
 *                         clow <=  C_i(x) <= cupp ;
 *                         xlow <=    x <= xupp ;
 *  </pre> 
 *
 *  The general linear inequality constraints must have either an upper
 *  or lower bound, but need not have both bounds. The variables may have 
 *  no bounds; an upper bound; a lower bound or both an upper and lower
 *  bound.
*/

class NlpGenFactory : public ProblemFormulation{
protected:
  LinearAlgebraPackage * la;
  /** number of elements in x */
  long long nx;
	
  /** number of rows in C_e and b */
  long long my;
	
  /** number of rows in C_i */
  long long mz;

  int nnzQ;
  int nnzA;
  int nnzC;


public:

  NlpGenFactory();

  NlpGenFactory( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_);

  virtual Residuals     * makeResiduals( Data * prob_in );
  virtual Variables     * makeVariables( Data * prob_in );


//  virtual Data          * makeData     ( );
//  virtual Data          * makeData     (in);


  Data  * makeData();
  Data  * makeData( double    c[],
			 int    krowQ[],  int  jcolQ[],  double dQ[],
			 double  xlow[],  char ixlow[],
			 double  xupp[],  char ixupp[],
			 int    krowA[],  int  jcolA[],  double dA[],
			 double     b[],
			 int    krowC[],  int  jcolC[],  double dC[],
			 double  clow[],  char iclow[],
			 double  cupp[],  char icupp[] );



  virtual LinearSystem  * makeLinsys( Data * prob_in );

  virtual ~NlpGenFactory() {};
};

#endif
