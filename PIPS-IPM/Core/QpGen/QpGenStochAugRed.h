#ifndef STOCHQPGENFACTORYAUGRED
#define STOCHQPGENFACTORYAUGRED

#include "QpGenStoch.h"

class QpGenStochAugRed : public QpGenStoch {
 public:
  QpGenStochAugRed( StochInputTree* );
 private:
  QpGenStochAugRed( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  QpGenStochAugRed();
 public:
  virtual ~QpGenStochAugRed();

  virtual QpGenStochLinsysRoot* newLinsysRoot();
  virtual QpGenStochLinsysRoot* newLinsysRoot(QpGenStochData* prob,
					      OoqpVector* dd,OoqpVector* dq,
					      OoqpVector* nomegaInv, OoqpVector* rhs);
};
#endif
