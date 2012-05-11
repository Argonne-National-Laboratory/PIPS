/* PIPS
   Authors: Miles Lubin and Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SAUGLINSYSEMTLSYM
#define SAUGLINSYSEMTLSYM

#include "sLinsysRoot.h"
#include "sLinsysRootAugEmtl.h"
#include "EmtlContext.h"


class QpGenStochData;
/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class sLinsysRootAugEmtlSym : public sLinsysRootAugEmtl {
 protected:
  //sLinsysRootAugEmtlSym() {};

  virtual DoubleLinearSolver* 
                       createSolver  (QpGenStochData* prob, 
				      SymMatrix* kktmat);

  //virtual void         createChildren(QpGenStochData* prob) 
  //{sLinsysRoot::createChildren(prob);};
 public:

  sLinsysRootAugEmtlSym(sFactory * factory_, QpGenStochData * prob_, const EmtlContext &ctx_);
  sLinsysRootAugEmtlSym(sFactory* factory,
			     QpGenStochData* prob_,				    
			     OoqpVector* dd_, OoqpVector* dq_, 
			     OoqpVector* nomegaInv_,
			     OoqpVector* rhs_,
			     const EmtlContext &ctx_);
  virtual ~sLinsysRootAugEmtlSym();

 public:
  virtual void factor2(QpGenStochData *prob, Variables *vars);
 protected:

};


#endif

