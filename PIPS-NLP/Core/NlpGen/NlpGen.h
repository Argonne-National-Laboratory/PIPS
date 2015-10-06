/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef NLPGENFACTORY
#define NLPGENFACTORY

#include "ProblemFormulation.h"

class Data;
class Residuals;
class LinearSystem;
class Variables;
class LinearAlgebraPackage;
class OoqpVector;

/**
 * @defgroup NlpGen
 * 
 * OOQP's default general problem formulation:
 *
 *  <pre>
 *  minimize    f(x)         ; 
 *  subject to                      A(x)  = b    ;     <-- dual y
 *                         clow <=  C(x) <= cupp ; 
 *                         xlow <=    x <= xupp ;
 *  </pre> 
 *
 *  The general linear equality constraints must have either an upper
 *  or lower bound, but need not have both bounds. The variables may have 
 *  no bounds; an upper bound; a lower bound or both an upper and lower
 *  bound.
 *
 *  add slacks as:
 *  t: slack for lower bound:  C(x)-t=clow
 *  u: slack for lower bound:  C(x)+u=cupp 
 *  lambda: the dual var for the lower bound constraint  C(x)-t-clow=0
 *  pi: dual var for the upper bound C(x)+u-cupp=0
 *
 *  Or:  add slacks as:
 *  s:  C(x)=s;
 *  t: slack for lower bound: s-t=clow
 *  u: slack for lower bound:  s+u=cupp 
 *  lambda: the dual var for the lower bound constraint  s-t-clow=0
 *  pi: dual var for the upper bound s+u-cupp=0
 *
 *  xlow <=    x <= xupp  => x-v=xlow, x+w=xupp
 *  v: slack for lower bound:  x-v=xlow
 *  w: slack for lower bound:  x+w=xupp
 *  gamma: dual var for the lower bound  
 *  phi: dual var for the upper bound  
 *
 *
*/
class NlpGen : public ProblemFormulation {
protected:
  LinearAlgebraPackage * la;
  /** number of elements in x */
  long long nx;

  /** number of rows in A and b */
  long long my;

  /** number of rows in C */
  long long mz;

  NlpGen(){};

  NlpGen( long long nx_, long long my_, long long mz_ );
public:
  virtual Residuals     * makeResiduals( Data * prob_in );
  virtual Variables     * makeVariables( Data * prob_in );

  virtual void joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			OoqpVector& rhs2_in, OoqpVector& rhs3_in ) = 0;

  virtual void separateVars( OoqpVector& x_in, OoqpVector& y_in,
			     OoqpVector& z_in, OoqpVector& vars_in) = 0;

  virtual ~NlpGen() {};

  virtual void joinRHSXSYZ( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			OoqpVector& rhs2_in, OoqpVector& rhs3_in, OoqpVector& rhs4_in ) = 0;

  virtual void separateVarsXSYZ( OoqpVector& x_in, OoqpVector& s_in, 
  			OoqpVector& y_in, OoqpVector& z_in, OoqpVector& vars_in) = 0;

  virtual void copyXSYZ_fromArray( OoqpVector& vec_xsyz, double* array_in, const int nb_col);
  virtual void copyXSYZ_toArray( OoqpVector& vec_xsyz, double* array_in, const int nb_col);  
  
};

#endif



