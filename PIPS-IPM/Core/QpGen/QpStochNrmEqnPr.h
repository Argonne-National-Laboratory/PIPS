#ifndef STOCHQPGENFACTORYNRMEQNPRCND
#define STOCHQPGENFACTORYNRMEQNPRCND

#include "QpGenStoch.h"

class QpStochNrmEqnPr : public QpGenStoch {
 public:
  QpStochNrmEqnPr( StochInputTree* );
 private:
  QpStochNrmEqnPr( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  QpStochNrmEqnPr();
 public:
  virtual ~QpStochNrmEqnPr();

  virtual QpGenStochLinsysRoot* newLinsysRoot();
  virtual QpGenStochLinsysRoot* newLinsysRoot(QpGenStochData* prob,
					      OoqpVector* dd,OoqpVector* dq,
					      OoqpVector* nomegaInv, OoqpVector* rhs);
};
#endif
