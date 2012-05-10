#ifndef STOCHQPGENFACTORYAUGREDPR_CG
#define STOCHQPGENFACTORYAUGREDPR_CG

#include "QpGenStochAugRedPr.h"

class QpStochAugRedPrCG : public QpGenStochAugRedPr {
 public:
  QpStochAugRedPrCG( StochInputTree* );
 private:
  QpStochAugRedPrCG( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  QpStochAugRedPrCG();
 public:
  virtual ~QpStochAugRedPrCG();

  virtual QpGenStochLinsysRoot* newLinsysRoot();
  virtual QpGenStochLinsysRoot* newLinsysRoot(QpGenStochData* prob,
					      OoqpVector* dd,OoqpVector* dq,
					      OoqpVector* nomegaInv, OoqpVector* rhs);
};
#endif
