/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SAUGLINSYS
#define SAUGLINSYS

#include "sLinsysRoot.h"
class sData;
/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class sLinsysRootAug : public sLinsysRoot {
 protected:
  sLinsysRootAug() {};

  virtual SymMatrix*   createKKT     (sData* prob);
  virtual DoubleLinearSolver* 
                       createSolver  (sData* prob, 
				      SymMatrix* kktmat);

  //virtual void         createChildren(sData* prob) 
  //{sLinsysRootAug::createChildren(prob);};
 public:

  sLinsysRootAug(sFactory * factory_, sData * prob_);
  sLinsysRootAug(sFactory* factory,
		 sData* prob_,				    
		 OoqpVector* dd_, OoqpVector* dq_, 
		 OoqpVector* nomegaInv_,
		 OoqpVector* rhs_);
  virtual ~sLinsysRootAug();

 public:
  virtual void finalizeKKT( sData* prob, Variables* vars);
 protected:
  virtual void solveReduced( sData *prob, SimpleVector& b);


  SymMatrix* CtDC;
  SimpleVector* redRhs;
};


#endif

