/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef SAUGLINSYS
#define SAUGLINSYS

#include "sLinsysRoot.h"
class sData;
/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class sLinsysRootAug : public sLinsysRoot {
 protected:
  sLinsysRootAug();

  virtual SymMatrix*   createKKT     (sData* prob);
  virtual DoubleLinearSolver* 
                       createSolver  (sData* prob, 
				      SymMatrix* kktmat);
  virtual void initialize(sFactory* factory, sData* prob);
  //virtual void         createChildren(sData* prob) 
  //{sLinsysRootAug::createChildren(prob);};
 public:

  sLinsysRootAug(sFactory * factory_, sData * prob_);
  sLinsysRootAug(sFactory* factory,
		 sData* prob_,				    
		 OoqpVector* dd_, OoqpVector* dq_, 
		 OoqpVector* nomegaInv_,
		 OoqpVector* rhs_,
		 OoqpVector* additiveDiag_);
  virtual ~sLinsysRootAug();

 public:
  virtual void finalizeKKT( sData* prob, Variables* vars);
 protected:
  virtual void solveReduced( sData *prob, SimpleVector& b);
  void solveWithIterRef( sData *prob, SimpleVector& b);
  void solveWithBiCGStab( sData *prob, SimpleVector& b);

  /** y = beta*y - alpha* SC * x */
  void SCmult ( double beta, SimpleVector& y, double alpha, SimpleVector& x, sData* prob);

  SymMatrix* CtDC;
  SimpleVector* redRhs;
};

#endif

