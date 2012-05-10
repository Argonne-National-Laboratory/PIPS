#ifndef STOCHQPGENFACTORYAUGREDPR
#define STOCHQPGENFACTORYAUGREDPR

#include "QpGenStoch.h"

class QpGenStochAugRedPr : public QpGenStoch {
 public:
  QpGenStochAugRedPr( StochInputTree* );
 protected:
  QpGenStochAugRedPr( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  QpGenStochAugRedPr();
 public:
  virtual ~QpGenStochAugRedPr();

  virtual QpGenStochLinsysRoot* newLinsysRoot();
  virtual QpGenStochLinsysRoot* newLinsysRoot(QpGenStochData* prob,
					      OoqpVector* dd,OoqpVector* dq,
					      OoqpVector* nomegaInv, OoqpVector* rhs);
};
#endif
