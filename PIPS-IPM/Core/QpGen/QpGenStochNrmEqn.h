#ifndef STOCHQPGENFACTORYNRMEQN
#define STOCHQPGENFACTORYNRMEQN

#include "QpGenStoch.h"

class QpGenStochNrmEqn : public QpGenStoch {
 public:
  QpGenStochNrmEqn( StochInputTree* );
 private:
  QpGenStochNrmEqn( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  QpGenStochNrmEqn();
 public:
  virtual ~QpGenStochNrmEqn();

  virtual QpGenStochLinsysRoot* newLinsysRoot();
  virtual QpGenStochLinsysRoot* newLinsysRoot(QpGenStochData* prob,
					      OoqpVector* dd,OoqpVector* dq,
					      OoqpVector* nomegaInv, OoqpVector* rhs);
};
#endif
