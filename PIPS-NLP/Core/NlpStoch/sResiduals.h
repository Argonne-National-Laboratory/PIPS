/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef STOCHRESIDUALS_NLP
#define STOCHRESIDUALS_NLP

#include "NlpGenResiduals.h"
#include <vector>

class sTree;
class StochVector;
/** 
 * Class added to supply a more generic constructor for its parent, NlpGenResiduals.
 * The default constructor of NlpGenResiduals can not be always used since it assumes that vectors
 * are arrays. For stochastic problems, the vectors are trees of arrays.
 *
 * @ingroup NlpGen
 */

class sResiduals : public NlpGenResiduals {
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
