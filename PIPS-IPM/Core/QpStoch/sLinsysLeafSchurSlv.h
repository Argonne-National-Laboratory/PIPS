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


  void addTermToDenseSchurCompl(sData *prob, 
				DenseSymMatrix& SC);

}; 
template<class LINSOLVER>
sLinsysLeafSchurSlv::sLinsysLeafSchurSlv(sFactory* factory,
					 sData* prob,
					 OoqpVector* dd_, 
					 OoqpVector* dq_, 
					 OoqpVector* nomegaInv_,
					 OoqpVector* rhs_,
					 LINSOLVER* s)
 : sLinsysLeaf(factory, prob, dd_, dq_, nomegaInv_, rhs_, s)
{
  
  
}


#endif
