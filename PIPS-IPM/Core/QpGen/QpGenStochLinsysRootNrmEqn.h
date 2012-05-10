#ifndef QPGENSTOCHROOTNRMLINSYS
#define QPGENSTOCHROOTNRMLINSYS

#include "QpGenStochLinsysRootAugRed.h"
class QpGenStochData;
/** 
 * ROOT (= NON-leaf) linear system in normal equations form
 */
class QpGenStochLinsysRootNrmEqn : public QpGenStochLinsysRootAugRed {
 protected:
  QpGenStochLinsysRootNrmEqn() {};

  virtual SymMatrix*   createKKT     (QpGenStochData* prob);
  virtual DoubleLinearSolver* 
                       createSolver  (QpGenStochData* prob, 
				      SymMatrix* kktmat);

  virtual void         createChildren(QpGenStochData* prob) 
  {QpGenStochLinsysRootAug::createChildren(prob);};
 public:

  QpGenStochLinsysRootNrmEqn(QpGenStoch * factory_, QpGenStochData * prob_);
  QpGenStochLinsysRootNrmEqn(QpGenStoch* factory,
			     QpGenStochData* prob_,				    
			     OoqpVector* dd_, OoqpVector* dq_, 
			     OoqpVector* nomegaInv_,
			     OoqpVector* rhs_);
  virtual ~QpGenStochLinsysRootNrmEqn();

 public:
  virtual void factor2 (QpGenStochData *prob, Variables *vars);

 protected:
  void updateKKT(QpGenStochData* prob, Variables* vars);
  void solveReduced( QpGenStochData *prob, SimpleVector& b);
  DenseSymMatrix* AQinvAt;
  DoubleLinearSolver* solver2;
  SimpleVector* redRhs;

 protected: 
  SimpleVector* rRhs2;
  //D = S';
  void myDenseFromSparseTrans(DenseGenMatrix& D, SparseGenMatrix& S);
};


#endif

