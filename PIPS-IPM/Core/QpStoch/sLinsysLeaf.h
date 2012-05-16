/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHLEAFLINSYS
#define STOCHLEAFLINSYS

#include "sLinsys.h"

class StochTree;
class sFactory;
class sData;

/** This class solves the linear system corresponding to a leaf node.
 *  It just redirects the call to QpGenSparseLinsys.
 */
class sLinsysLeaf : public sLinsys
{
 public:
  //sLinsysLeaf(QpGenStoch * factory_, sData * prob_);
  sLinsysLeaf(sFactory* factory,
		       sData* prob_,				    
		       OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_,
		       OoqpVector* rhs_);

  virtual ~sLinsysLeaf();

  virtual void factor2( sData *prob, Variables *vars);
  virtual void Lsolve ( sData *prob, OoqpVector& x );
  virtual void Dsolve ( sData *prob, OoqpVector& x );
  virtual void Ltsolve( sData *prob, OoqpVector& x );

  //virtual void Lsolve2 ( OoqpVector& x );
  //virtual void Dsolve2 ( OoqpVector& x );
  virtual void Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp);

  virtual void putZDiagonal( OoqpVector& zdiag );
  //virtual void solveCompressed( OoqpVector& rhs );
  virtual void putXDiagonal( OoqpVector& xdiag_ );

  //void Ltsolve_internal(  sData *prob, StochVector& x, SimpleVector& xp);
  void sync();
  virtual void deleteChildren();
 protected:
  sLinsysLeaf() {};
}; 

#endif
