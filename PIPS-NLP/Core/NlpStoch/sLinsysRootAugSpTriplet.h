/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SAUGLINSYSSPTRIP
#define SAUGLINSYSSPTRIP

#include "sLinsysRootAug.h"
#include "SparseSymMatrixRowMajList.h"
class sData;
/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class sLinsysRootAugSpTriplet : public sLinsysRootAug {
 protected:
  sLinsysRootAugSpTriplet() {};

  virtual SymMatrix*   createKKT(sData* prob);
  virtual DoubleLinearSolver* createSolver(sData* prob, 
					   SymMatrix* kktmat);
  
 public:

  sLinsysRootAugSpTriplet(sFactory * factory_, sData * prob_);
  sLinsysRootAugSpTriplet(sFactory* factory,
		 sData* prob_,				    
		 OoqpVector* dd_, OoqpVector* dq_, 
		 OoqpVector* nomegaInv_,
		 OoqpVector* rhs_,
		 OoqpVector* additiveDiag_);
  virtual ~sLinsysRootAugSpTriplet();

 public:
  virtual int factor2(sData *prob, Variables *vars);

  
  virtual void initializeKKT(sData* prob, Variables* vars);
  virtual void reduceKKT();
  //virtual int factorizeKKT();
  virtual void finalizeKKT(sData* prob, Variables* vars);

  virtual void UpdateMatrices( Data * prob_in,int const updateLevel=2);
 protected:
  
};

#endif

