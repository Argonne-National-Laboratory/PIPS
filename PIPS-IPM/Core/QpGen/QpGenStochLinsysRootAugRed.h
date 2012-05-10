#ifndef QPGENSTOCHROOTAUGREDLINSYS
#define QPGENSTOCHROOTAUGREDLINSYS

#include "QpGenStochLinsysRootAug.h"
class QpGenStochData;
/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class QpGenStochLinsysRootAugRed : public QpGenStochLinsysRootAug {
 protected:
  QpGenStochLinsysRootAugRed() {};

  virtual SymMatrix*   createKKT     (QpGenStochData* prob);
  virtual DoubleLinearSolver* 
                       createSolver  (QpGenStochData* prob, 
				      SymMatrix* kktmat);

  virtual void         createChildren(QpGenStochData* prob) 
  {QpGenStochLinsysRootAug::createChildren(prob);};
 public:

  QpGenStochLinsysRootAugRed(QpGenStoch * factory_, QpGenStochData * prob_);
  QpGenStochLinsysRootAugRed(QpGenStoch* factory,
			     QpGenStochData* prob_,				    
			     OoqpVector* dd_, OoqpVector* dq_, 
			     OoqpVector* nomegaInv_,
			     OoqpVector* rhs_);
  virtual ~QpGenStochLinsysRootAugRed();

 public:
  virtual void factor2 (QpGenStochData *prob, Variables *vars);
  // override Dsolve  of th parent
  virtual void Dsolve( QpGenStochData *prob, OoqpVector& x );
 protected:
  virtual void updateKKT(QpGenStochData* prob, Variables* vars);
  virtual void solveReduced( QpGenStochData *prob, SimpleVector& b);

  void submatrixAllReduce(DenseSymMatrix* A, 
			  int row, int col, int drow, int dcol,
			  MPI_Comm comm);

  SymMatrix* CtDC;
  SimpleVector* redRhs;
};


#endif

