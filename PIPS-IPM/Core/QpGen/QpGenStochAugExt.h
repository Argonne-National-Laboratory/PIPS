#ifndef STOCHQPGENFACTORYAUGEXT
#define STOCHQPGENFACTORYAUGEXT

#include "QpGenStoch.h"

class QpGenStochAugExt : public QpGenStoch {
 public:
  QpGenStochAugExt( StochInputTree* );
 private:
  QpGenStochAugExt( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  QpGenStochAugExt();
 public:
  virtual ~QpGenStochAugExt();

  virtual QpGenStochLinsysRoot* newLinsysRoot();
  virtual QpGenStochLinsysRoot* newLinsysRoot(QpGenStochData* prob,
					      OoqpVector* dd,OoqpVector* dq,
					      OoqpVector* nomegaInv, OoqpVector* rhs);
};
#endif
