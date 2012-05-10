#ifndef QPGENSTOCHROOTAUGREDPRECLINSYS
#define QPGENSTOCHROOTAUGREDPRECLINSYS

#include "QpGenStochLinsysRootAugRed.h"
#include "RemoteMatTimesVec.h"
#include "BiCGStabSolver.h"

class QpGenStochData;

/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form that make use 
 * of a preconditioned BICGSTAB to solve the root system
 */
class QpGenStochLinsysRootAugRedPrecond : public QpGenStochLinsysRootAugRed {
 protected:
  QpGenStochLinsysRootAugRedPrecond() {};

  virtual void         createChildren(QpGenStochData* prob);

  virtual SymMatrix*   createKKT     (QpGenStochData* prob);
  virtual DoubleLinearSolver* 
                       createSolver  (QpGenStochData* prob, 
				      SymMatrix* kktmat);
 public:

  QpGenStochLinsysRootAugRedPrecond(QpGenStoch * factory_, QpGenStochData * prob_);
  QpGenStochLinsysRootAugRedPrecond(QpGenStoch* factory,
			     QpGenStochData* prob_,				    
			     OoqpVector* dd_, OoqpVector* dq_, 
			     OoqpVector* nomegaInv_,
			     OoqpVector* rhs_);
  virtual ~QpGenStochLinsysRootAugRedPrecond();

 public:
  virtual void factor2 (QpGenStochData *prob, Variables *vars);
  // override Dsolve  of the parent
  virtual void Dsolve( QpGenStochData *prob, OoqpVector& x );
  //a different initialization of U'*V is needed for preconditioner
  void initializeUtV();
 protected:
  enum NodeType { eWorker,eSpecialWorker,ePrecond } ;

  RemoteMatTimesVec* Pmult;
  StoredMatTimesVec* Amult;

  int iErr;
 protected:
  NodeType me;
  NodeType whoAmI();
 public: 
  double* tmpVec1;
};


#endif

