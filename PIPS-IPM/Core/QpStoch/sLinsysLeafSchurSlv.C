/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysLeafSchurSlv.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"
#include "PardisoSolver.h"
#include "Ma57Solver.h"

static Ma57Solver* s=NULL;

//also templatize this constructor

sLinsysLeafSchurSlv::sLinsysLeafSchurSlv(sFactory* factory,
					 sData* prob,
					 OoqpVector* dd_, 
					 OoqpVector* dq_, 
					 OoqpVector* nomegaInv_,
					 OoqpVector* rhs_)
 : sLinsysLeaf(factory, prob, dd_, dq_, nomegaInv_, rhs_, s)
{
  
  
}


void sLinsysLeafSchurSlv::addTermToDenseSchurCompl( sData *prob, 
						    DenseSymMatrix& SC)
{

  cout  << " sLinsysLeafSchurSlv::addTermToDenseSchurComp is apparently called!" << endl;
  //sLinsysLeaf::addTermToDenseSchurCompl( prob, SC);
}
