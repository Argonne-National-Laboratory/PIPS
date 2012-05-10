#ifndef QPGENSTOCHLEAFLINSYS
#define QPGENSTOCHLEAFLINSYS

#include "QpGenStochLinsys.h"

class StochTree;
class QpGenStoch;
class QpGenStochData;

/** This class solves the linear system corresponding to a leaf node.
 *  It just redirects the call to QpGenSparseLinsys.
 */
class QpGenStochLinsysLeaf : public QpGenStochLinsys
{
 public:
  //QpGenStochLinsysLeaf(QpGenStoch * factory_, QpGenStochData * prob_);
  QpGenStochLinsysLeaf(QpGenStoch* factory,
		       QpGenStochData* prob_,				    
		       OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_,
		       OoqpVector* rhs_);

  virtual ~QpGenStochLinsysLeaf();

  virtual void factor2( QpGenStochData *prob, Variables *vars);
  virtual void Lsolve ( QpGenStochData *prob, OoqpVector& x );
  virtual void Dsolve ( QpGenStochData *prob, OoqpVector& x );
  virtual void Ltsolve( QpGenStochData *prob, OoqpVector& x );

  //virtual void Lsolve2 ( OoqpVector& x );
  //virtual void Dsolve2 ( OoqpVector& x );
  virtual void Ltsolve2( QpGenStochData *prob, StochVector& x, SimpleVector& xp);

  virtual void putZDiagonal( OoqpVector& zdiag );
  //virtual void solveCompressed( OoqpVector& rhs );
  virtual void putXDiagonal( OoqpVector& xdiag_ );

  //void Ltsolve_internal(  QpGenStochData *prob, StochVector& x, SimpleVector& xp);
  void sync();
  virtual void deleteChildren();
 protected:
  QpGenStochLinsysLeaf() {};
}; 

#endif
