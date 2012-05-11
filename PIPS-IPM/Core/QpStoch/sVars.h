#ifndef QPGENSTOCHVARS
#define QPGENSTOCHVARS

#include "QpGenVars.h"
#include "OoqpVectorHandle.h"
#include <vector>

class QpGen;
class QpGenData;
class LinearAlgebraPackage;
class StochTree;

class sVars : public QpGenVars {
 public:
  sVars(StochTree* tree,
	OoqpVector * ixlow_in, OoqpVector * ixupp_in,
	OoqpVector * iclow_in, OoqpVector * icupp_in);

  /** constructor in which the data and variable pointers are set to
      point to the given arguments */
  sVars( StochTree* tree,
	 OoqpVector * x_in, OoqpVector * s_in,
	 OoqpVector * y_in, OoqpVector * z_in,
	 OoqpVector * v_in, OoqpVector * gamma_in,
	 OoqpVector * w_in, OoqpVector * phi_in,
	 OoqpVector * t_in, OoqpVector * lambda_in,
	 OoqpVector * u_in, OoqpVector * pi_in,
	 OoqpVector * ixlow_in, int nxlowGlobal,
	 OoqpVector * ixupp_in, int nxuppGlobal,
	 OoqpVector * iclow_in, int mclowGlobal,
	 OoqpVector * icupp_in, int mcuppGlobal);

  virtual ~sVars();
  
  virtual void sync();
 protected:
  void createChildren();
  std::vector<sVars*> children;
  void AddChild(sVars* child);

  StochTree* stochNode;
};
#endif

