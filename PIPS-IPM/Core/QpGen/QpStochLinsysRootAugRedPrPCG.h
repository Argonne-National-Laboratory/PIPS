#ifndef QPSTOCHROOTAUGREDLINSYSPCG
#define QPSTOCHROOTAUGREDLINSYSPCG

#include "QpGenStochLinsysRootAugRed.h"
#include "PCGSolver.h"

class QpGenStochData;
class RemoteMatTimesVec;
class Ma57Solver;
/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class QpStochLinsysRootAugRedPrPCG : public QpGenStochLinsysRootAugRed {
 protected:
  QpStochLinsysRootAugRedPrPCG() {};

  virtual SymMatrix*   createKKT     (QpGenStochData* prob);
  virtual DoubleLinearSolver* 
                       createSolver  (QpGenStochData* prob, 
				      SymMatrix* kktmat);
  virtual void         updateKKT     (QpGenStochData* prob, 
				      Variables* vars);

 public:

  QpStochLinsysRootAugRedPrPCG(QpGenStoch * factory_, QpGenStochData * prob_);
  QpStochLinsysRootAugRedPrPCG(QpGenStoch* factory,
			     QpGenStochData* prob_,				    
			     OoqpVector* dd_, OoqpVector* dq_, 
			     OoqpVector* nomegaInv_,
			     OoqpVector* rhs_);
  virtual ~QpStochLinsysRootAugRedPrPCG();

 public:
  virtual void factor2 (QpGenStochData *prob, Variables *vars);
  // override Dsolve  of th parent
  virtual void Dsolve( QpGenStochData *prob, OoqpVector& x );

  void solveReduced( QpGenStochData *prob, SimpleVector& b);
 protected:
  //virtual void updateKKT(QpGenStochData* prob, Variables* vars);
  //virtual void solveReduced( QpGenStochData *prob, SimpleVector& b);

  enum NodeType { eWorker,eSpecialWorker,ePrecond } ;

  RemoteMatTimesVec* Pmult;
  StoredMatTimesVec* Qmult;
  StoredMatTransTimesVec* Atmult;
  int iErr;

 protected:
  NodeType me;
  NodeType whoAmI();

  void workerReduce_source(DenseSymMatrix* M, int target, MPI_Comm comm);
  void workerReduce_target(DenseSymMatrix* M, int target, MPI_Comm comm);
  void precndReduce_source(DenseSymMatrix* M, int target, MPI_Comm comm);
  void precndReduce_target(DenseSymMatrix* M, int target, MPI_Comm comm);
 public: 
  double* tmpVec1;
  SymMatrix* AAt;
  Ma57Solver* AAtSolver;
};


#endif

