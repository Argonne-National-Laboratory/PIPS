/* PIPS
   Authors: Miles Lubin and Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SAUGLINSYSEMTL
#define SAUGLINSYSEMTL

#include "sLinsysRoot.h"
#include "EmtlContext.h"

class QpGenStochData;
/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class sLinsysRootAugEmtl : public sLinsysRoot {
 protected:
  //sLinsysRootAugEmtl() {};

  virtual SymMatrix*   createKKT     (QpGenStochData* prob);
  virtual DoubleLinearSolver* 
                       createSolver  (QpGenStochData* prob, 
				      SymMatrix* kktmat);

  //virtual void         createChildren(QpGenStochData* prob) 
  //{sLinsysRoot::createChildren(prob);};
 public:

  sLinsysRootAugEmtl(sFactory * factory_, QpGenStochData * prob_, const EmtlContext &ctx_);
  sLinsysRootAugEmtl(sFactory* factory,
			     QpGenStochData* prob_,				    
			     OoqpVector* dd_, OoqpVector* dq_, 
			     OoqpVector* nomegaInv_,
			     OoqpVector* rhs_,
			     const EmtlContext &ctx_);
  virtual ~sLinsysRootAugEmtl();

 public:
  virtual void finalizeKKT(QpGenStochData* prob, Variables* vars);
  virtual void initializeKKT(QpGenStochData* prob, Variables* vars);
  virtual void factor2(QpGenStochData *prob, Variables *vars);
 protected:
  virtual void solveReduced( QpGenStochData *prob, SimpleVector& b);

  const EmtlContext &ctx;
  SymMatrix* CtDC;
  SimpleVector* redRhs;
  OoqpVector* emtlRhs;
};


#endif

