#ifndef QPGENSTOCHROOTAUGREDPRECLINSYS_CG
#define QPGENSTOCHROOTAUGREDPRECLINSYS_CG

#include "QpGenStochLinsysRootAugRedPrecond.h"
#include "RemoteMatTimesVec.h"
#include "CGSolver.h"

class QpGenStochData;

/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form that make use 
 * of a preconditioned BICGSTAB to solve the root system
 */
class QpStochLinsysRootAugRedPrCG : public QpGenStochLinsysRootAugRedPrecond {

 protected:
  QpStochLinsysRootAugRedPrCG() {};

  virtual SymMatrix*   createKKT     (QpGenStochData* prob);
  virtual DoubleLinearSolver* 
                       createSolver  (QpGenStochData* prob, 
				      SymMatrix* kktmat);
 public:

  QpStochLinsysRootAugRedPrCG(QpGenStoch * factory_, QpGenStochData * prob_);
  QpStochLinsysRootAugRedPrCG(QpGenStoch* factory,
			     QpGenStochData* prob_,				    
			     OoqpVector* dd_, OoqpVector* dq_, 
			     OoqpVector* nomegaInv_,
			     OoqpVector* rhs_);
  virtual ~QpStochLinsysRootAugRedPrCG();
};


#endif
