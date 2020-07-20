#ifndef STOCHRESIDUALS
#define STOCHRESIDUALS

#include "QpGenResiduals.h"
#include "StochVector_fwd.h"
#include <vector>

class sTree;

/** 
 * Class added to supply a more generic constructor for its parent, QpGenResiduals.
 * The default constructor of QpGenResiduals can not be always used since it assumes that vectors
 * are arrays. For stochastic problems, the vectors are trees of arrays.
 *
 * @ingroup QpGen
 */

class sResiduals : public QpGenResiduals {
public:
  /**
   * Constructor
   */
  sResiduals( sTree* tree,OoqpVector * rQ, 
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
  
  sResiduals( sTree* tree,
	      OoqpVector * ixlow_, OoqpVector * ixupp_,
	      OoqpVector * iclow_, OoqpVector * icupp_ );
  
  sResiduals( const sResiduals& res );

  sResiduals( const sResiduals& res,
        OoqpVectorHandle ixlow_, OoqpVectorHandle ixupp_,
        OoqpVectorHandle iclow_, OoqpVectorHandle icupp_ );

  virtual void sync();
 private:
  std::vector<sResiduals*> children;
  void createChildren();
  void destroyChildren();
  void AddChild(sResiduals* child);

  
 protected:
  sTree* stochNode;
};

#endif
