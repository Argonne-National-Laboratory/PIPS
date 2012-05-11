#ifndef QPGENRESIDUALS2
#define QPGENRESIDUALS2

#include "QpGenResiduals.h"
#include <vector>

class StochTree;
class StochVector;
/** 
 * Class added to supply a more generic constructor for its parent, QpGenResiduals.
 * The default constructor of QpGenResiduals can not be always used since it assumes that vectors
 * are arrays. For stochastic problems, the vectors are trees of arrays.
 *
 * @ingroup QpGen
 */

class QpGenResiduals2 : public QpGenResiduals {
public:
  /**
   * Constructor
   */
  QpGenResiduals2(StochTree* tree,OoqpVector * rQ, 
		  OoqpVector * rA, OoqpVector * rC, 
		  OoqpVector * rz, 
		  OoqpVector * rt, OoqpVector * rlambda, 
		  OoqpVector * ru, OoqpVector * rpi, 
		  OoqpVector * rv, OoqpVector * rgamma, 
		  OoqpVector * rw, OoqpVector * rphi, 
		  OoqpVector * ixlow, double nxlowGlobal,
		  OoqpVector * ixupp, double nxuppGlobal,
		  OoqpVector * iclow, double mclowGlobal, 
		  OoqpVector * icupp, double mcuppGlobal );

  QpGenResiduals2(StochTree* tree,
		  OoqpVector * ixlow_, OoqpVector * ixupp_,
		  OoqpVector * iclow_, OoqpVector * icupp_ );

  virtual void sync();
 private:
  std::vector<QpGenResiduals2*> children;
  void createChildren();
  void destroyChildren();
  void AddChild(QpGenResiduals2* child);

  
 protected:
  //QpGenResiduals* residuals;
  StochTree* stochNode;
};

#endif
