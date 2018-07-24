/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHLEAFLINSYS_SCHURSLV
#define STOCHLEAFLINSYS_SCHURSLV

#include "sLinsysLeaf.h"

class StochTree;
class sFactory;
class sData;

/** This class solves the linear system corresponding to a leaf node.
 *  It just redirects the call to QpGenSparseLinsys.
 */
class sLinsysLeafSchurSlv : public sLinsysLeaf
{
 public:
  template<class LINSOLVER>
  sLinsysLeafSchurSlv(sFactory* factory,
		      sData* prob_,				    
		      OoqpVector* dd_, OoqpVector* dq_, 
		      OoqpVector* nomegaInv_,
		      OoqpVector* rhs_, LINSOLVER* solver);

  void factor2(sData *prob, Variables *vars);
  void addTermToDenseSchurCompl(sData *prob, 
				DenseSymMatrix& SC);
  void addTermToSparseSchurCompl(sData *prob,
            SparseSymMatrix& SC);

 private:
  bool switchedToSafeSlv;

}; 
template<class LINSOLVER>
sLinsysLeafSchurSlv::sLinsysLeafSchurSlv(sFactory* factory,
					 sData* prob,
					 OoqpVector* dd_, 
					 OoqpVector* dq_, 
					 OoqpVector* nomegaInv_,
					 OoqpVector* rhs_,
					 LINSOLVER* s)
: sLinsysLeaf(factory, prob, dd_, dq_, nomegaInv_, rhs_, s), 
  switchedToSafeSlv(false)
{
  
  
}


#endif
