#ifndef QPGENSTOCHROOTNRMLINSYSPRECOND
#define QPGENSTOCHROOTNRMLINSYSPRECOND

#include "QpGenStochLinsysRootNrmEqn.h"
#include "RemoteMatTimesVec.h"
#include "CGSolver.h"

class QpGenStochData;
/** 
 * ROOT (= NON-leaf) linear system in normal equations formthat make use 
 * of a preconditioned CG to solve the root system
 */
class QpGenStochLinsysRootNrmEqnPrecond : public QpGenStochLinsysRootNrmEqn {
 protected:
  QpGenStochLinsysRootNrmEqnPrecond() {};

  virtual void         createChildren(QpGenStochData* prob);
  virtual SymMatrix*   createKKT     (QpGenStochData* prob);
  virtual DoubleLinearSolver* 
                       createSolver  (QpGenStochData* prob, 
				      SymMatrix* kktmat);

 public:

  QpGenStochLinsysRootNrmEqnPrecond(QpGenStoch * factory_, QpGenStochData * prob_);
  QpGenStochLinsysRootNrmEqnPrecond(QpGenStoch* factory,
			     QpGenStochData* prob_,				    
			     OoqpVector* dd_, OoqpVector* dq_, 
			     OoqpVector* nomegaInv_,
			     OoqpVector* rhs_);
  virtual ~QpGenStochLinsysRootNrmEqnPrecond();

 public:
  virtual void factor2 (QpGenStochData *prob, Variables *vars);
  // override Dsolve  of the parent
  virtual void Dsolve( QpGenStochData *prob, OoqpVector& x );
  //a different initialization of U'*V is needed for preconditioner
  void initializeUtV();
 protected:
  enum NodeType { eWorker,eSpecialWorker,ePrecond } ;
 protected: 
  RemoteMatTimesVec* Pmult;
  StoredMatTimesVec* Amult;
  
  int iErr;
 private:
  NodeType me;
  NodeType whoAmI(); 
 public: 
  double* tmpVec1;
};


#endif

