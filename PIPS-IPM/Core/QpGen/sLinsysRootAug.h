/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SAUGLINSYS
#define SAUGLINSYS

#include "sLinsysRoot.h"
class QpGenStochData;
/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class sLinsysRootAug : public sLinsysRoot {
 protected:
  sLinsysRootAug() {};

  virtual SymMatrix*   createKKT     (QpGenStochData* prob);
  virtual DoubleLinearSolver* 
                       createSolver  (QpGenStochData* prob, 
				      SymMatrix* kktmat);

  //virtual void         createChildren(QpGenStochData* prob) 
  //{sLinsysRootAug::createChildren(prob);};
 public:

  sLinsysRootAug(sFactory * factory_, QpGenStochData * prob_);
  sLinsysRootAug(sFactory* factory,
			     QpGenStochData* prob_,				    
			     OoqpVector* dd_, OoqpVector* dq_, 
			     OoqpVector* nomegaInv_,
			     OoqpVector* rhs_);
  virtual ~sLinsysRootAug();

 public:
  virtual void finalizeKKT(QpGenStochData* prob, Variables* vars);
 protected:
  virtual void solveReduced( QpGenStochData *prob, SimpleVector& b);


  SymMatrix* CtDC;
  SimpleVector* redRhs;
};


#endif

