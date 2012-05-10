#ifndef QPGENSTOCHROOTAUGLINSYS
#define QPGENSTOCHROOTAUGLINSYS

#include "QpGenStochLinsysRoot.h"
class QpGenStochData;
/** 
 * ROOT (= NON-leaf) linear system in ABSTRACT augmented form 
 */
class QpGenStochLinsysRootAug : public QpGenStochLinsysRoot {
 protected:
  QpGenStochLinsysRootAug() {};

  virtual SymMatrix*   createKKT     (QpGenStochData* prob)=0;
  virtual DoubleLinearSolver* 
                       createSolver  (QpGenStochData* prob, 
				      SymMatrix* kktmat)=0;
  virtual void         createChildren(QpGenStochData* prob) 
  {QpGenStochLinsysRoot::createChildren(prob);};
 public:

  QpGenStochLinsysRootAug(QpGenStoch * factory_, QpGenStochData * prob_);
  QpGenStochLinsysRootAug(QpGenStoch* factory,
			     QpGenStochData* prob_,				    
			     OoqpVector* dd_, OoqpVector* dq_, 
			     OoqpVector* nomegaInv_,
			     OoqpVector* rhs_);

  virtual void factor2 (QpGenStochData *prob, Variables *vars) = 0;
  virtual void Lsolve  ( QpGenStochData *prob, OoqpVector& x );
  virtual void Dsolve  ( QpGenStochData *prob, OoqpVector& x );
  virtual void Ltsolve ( QpGenStochData *prob, OoqpVector& x );
  virtual void Ltsolve2( QpGenStochData *prob, StochVector& x, SimpleVector& xp);

 public:
  virtual ~QpGenStochLinsysRootAug();
 protected: //buffers
  double* dcolG; int ncolG;
  double* new_dcolG(int nDim);

  DenseSymMatrix* UtV;
  void allocUtV();
  virtual void initializeUtV();
  void deleteUtV();

  //faster than DenseSymMatrix::atPutZeros
  void myAtPutZeros(DenseSymMatrix* mat, 
		    int row, int col, 
		    int rowExtent, int colExtent);
};


#endif

