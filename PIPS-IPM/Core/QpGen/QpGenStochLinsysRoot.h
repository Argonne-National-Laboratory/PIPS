#ifndef QPGENSTOCHROOTLINSYS
#define QPGENSTOCHROOTLINSYS

#include "QpGenStochLinsys.h"
class StochTree;
class QpGenStoch;
class QpGenStochData;

/** 
 * ROOT (= NON-leaf) linear system
 */
class QpGenStochLinsysRoot : public QpGenStochLinsys {
 protected:
  QpGenStochLinsysRoot() {};

  virtual void         createChildren(QpGenStochData* prob);
  virtual void         deleteChildren();

  virtual SymMatrix*   createKKT     (QpGenStochData* prob) = 0;
  virtual DoubleLinearSolver* 
                       createSolver  (QpGenStochData* prob, 
				      SymMatrix* kktmat) = 0;
 public:
  std::vector<QpGenStochLinsys*> children;
  int iAmDistrib;
 public:

  QpGenStochLinsysRoot(QpGenStoch * factory_, QpGenStochData * prob_);
  QpGenStochLinsysRoot(QpGenStoch* factory,
		       QpGenStochData* prob_,				    
		       OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_,
		       OoqpVector* rhs_);

  virtual void factor2(QpGenStochData *prob, Variables *vars)=0;

  virtual void Lsolve ( QpGenStochData *prob, OoqpVector& x )=0;
  virtual void Dsolve ( QpGenStochData *prob, OoqpVector& x )=0;
  virtual void Ltsolve( QpGenStochData *prob, OoqpVector& x )=0;

  //virtual void Lsolve2 ( OoqpVector& x );
  //virtual void Dsolve2 ( OoqpVector& x );
  virtual void Ltsolve2( QpGenStochData *prob, StochVector& x, SimpleVector& xp)=0;


  virtual void putXDiagonal( OoqpVector& xdiag_ );
  virtual void putZDiagonal( OoqpVector& zdiag );
 
  virtual void AddChild(QpGenStochLinsys* child);

  void sync();
 public:
  virtual ~QpGenStochLinsysRoot();
 protected: //buffers

  OoqpVector* zDiag;
  OoqpVector* xDiag;

#ifdef STOCH_TESTING
 protected: 
  static void dumpRhs(int proc, const char* nameToken,  SimpleVector& rhs);
  static void dumpMatrix(int scen, int proc, const char* nameToken, DenseSymMatrix& M);
#endif
};

class DummyLinearSolver : public DoubleLinearSolver {

  void diagonalChanged( int idiag, int extent ) {};
  void matrixChanged() {};
  void solve ( OoqpVector& x ) {};
  void Lsolve  ( OoqpVector& x ) {};
  void Dsolve  ( OoqpVector& x ) {};
  void Ltsolve ( OoqpVector& x ) {};
};


#endif

