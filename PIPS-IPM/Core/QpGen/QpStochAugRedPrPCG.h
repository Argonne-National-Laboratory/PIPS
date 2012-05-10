#ifndef STOCHQPGENFACTORYAUGREDPR_CG_PROJ
#define STOCHQPGENFACTORYAUGREDPR_CG_PROJ

#include "QpGenStochAugRedPr.h"

class QpStochAugRedPrPCG : public QpGenStochAugRedPr {
 public:
  QpStochAugRedPrPCG( StochInputTree* );
 private:
  QpStochAugRedPrPCG( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  QpStochAugRedPrPCG();
 public:
  virtual ~QpStochAugRedPrPCG();

  virtual QpGenStochLinsysRoot* newLinsysRoot();
  virtual QpGenStochLinsysRoot* newLinsysRoot(QpGenStochData* prob,
					      OoqpVector* dd,OoqpVector* dq,
					      OoqpVector* nomegaInv, OoqpVector* rhs);
};
#endif
