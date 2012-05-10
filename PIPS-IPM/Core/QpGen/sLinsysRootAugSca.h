/* PIPS
   Authors: Miles Lubin and Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SAUGLINSYSSCA
#define SAUGLINSYSSCA

#include "sLinsysRoot.h"
#include "scalapack.h"

class QpGenStochData;
/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class sLinsysRootAugSca : public sLinsysRoot {
 protected:
  sLinsysRootAugSca() {};

  virtual SymMatrix*   createKKT     (QpGenStochData* prob);
  virtual DoubleLinearSolver* 
                       createSolver  (QpGenStochData* prob, 
				      SymMatrix* kktmat);

  //virtual void         createChildren(QpGenStochData* prob) 
  //{sLinsysRoot::createChildren(prob);};
 public:

  sLinsysRootAugSca(sFactory * factory_, QpGenStochData * prob_, COMMINFO &cinfo_);
  sLinsysRootAugSca(sFactory* factory,
			     QpGenStochData* prob_,				    
			     OoqpVector* dd_, OoqpVector* dq_, 
			     OoqpVector* nomegaInv_,
			     OoqpVector* rhs_,
			     COMMINFO &cinfo_);
  virtual ~sLinsysRootAugSca();

 public:
  virtual void finalizeKKT(QpGenStochData* prob, Variables* vars);
  virtual void initializeKKT(QpGenStochData* prob, Variables* vars);
  virtual void factor2(QpGenStochData *prob, Variables *vars);
 protected:
  virtual void solveReduced( QpGenStochData *prob, SimpleVector& b);

  COMMINFO cinfo;
  SymMatrix* CtDC;
  SimpleVector* redRhs;
  OoqpVector* scaRhs;
};


#endif

