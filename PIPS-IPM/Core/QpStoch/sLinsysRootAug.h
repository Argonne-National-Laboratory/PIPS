/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SAUGLINSYS
#define SAUGLINSYS

#include "sLinsysRoot.h"
//#include "pipschecks.h"

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
  virtual void solveReducedLinkCons( sData *prob, SimpleVector& b);
  virtual void finalizeKKTdense( sData* prob, Variables* vars);
  virtual void finalizeKKTsparse( sData* prob, Variables* vars);

  void solveWithIterRef( sData *prob, SimpleVector& b);
  void solveWithBiCGStab( sData *prob, SimpleVector& b);

  /** y = beta*y - alpha* SC * x */
  void SCmult ( double beta, SimpleVector& y, double alpha, SimpleVector& x, sData* prob);

  SymMatrix* CtDC;
  SimpleVector* redRhs;
};

#endif

