#ifndef QPGENSTOCHROOTAUGEXTLINSYS
#define QPGENSTOCHROOTAUGEXTLINSYS

#include "QpGenStochLinsysRootAug.h"

class QpGenStochData;
/** 
 * ROOT (= NON-leaf) linear system in big augmented form
 */
class QpGenStochLinsysRootAugExt : public QpGenStochLinsysRootAug {
 protected:
  QpGenStochLinsysRootAugExt() {};

  virtual SymMatrix*   createKKT     (QpGenStochData* prob);
  virtual DoubleLinearSolver* 
                       createSolver  (QpGenStochData* prob, 
				      SymMatrix* kktmat);
 public:

  QpGenStochLinsysRootAugExt(QpGenStoch * factory_, QpGenStochData * prob_);
  QpGenStochLinsysRootAugExt(QpGenStoch* factory,
			     QpGenStochData* prob_,				    
			     OoqpVector* dd_, OoqpVector* dq_, 
			     OoqpVector* nomegaInv_,
			     OoqpVector* rhs_);
  virtual ~QpGenStochLinsysRootAugExt();

 public:
  virtual void factor2 (QpGenStochData *prob, Variables *vars);

};

#endif

