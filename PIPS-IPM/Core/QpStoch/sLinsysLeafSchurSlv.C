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


void sLinsysLeafSchurSlv::addTermToDenseSchurCompl( sData *prob, 
						    DenseSymMatrix& SC)
{

  cout  << " sLinsysLeafSchurSlv::addTermToDenseSchurComp is apparently called!" << endl;
  sLinsysLeaf::addTermToDenseSchurCompl( prob, SC);
}
