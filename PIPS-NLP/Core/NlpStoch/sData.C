/* PIPS-IPM                                                             
 * Author: Cosmin G. Petra
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
 
/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "sData.h"
#include "sTree.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "NlpGenVars.h"
#include "SparseLinearAlgebraPackage.h"


#include "NlpInfo.h"

sData::sData(sTree* tree)
//  : NlpGenData(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
{
  stochNode = tree;
  H     = SymMatrixHandle ( tree->createQ() );
  grad     = OoqpVectorHandle    ( tree->createc() );

  blx   = OoqpVectorHandle    ( tree->createxlow()  );
  ixlow = OoqpVectorHandle    ( tree->createixlow() );
  bux   = OoqpVectorHandle    ( tree->createxupp()  );
  ixupp = OoqpVectorHandle    ( tree->createixupp() );


  Jeq = GenMatrixHandle        ( tree->createA() );
  bA  = OoqpVectorHandle       ( tree->createb() );


  Jineq = GenMatrixHandle     ( tree->createC() );  
  bl    = OoqpVectorHandle    ( tree->createclow() );
  iclow = OoqpVectorHandle    ( tree->createiclow() );
  bu    = OoqpVectorHandle    ( tree->createcupp()  );
  icupp = OoqpVectorHandle    ( tree->createicupp() );

  CeqBody  = OoqpVectorHandle ( tree->createCeqBody()	);
  CIneqBody  = OoqpVectorHandle ( tree->createCineqBody() );  


  trialBarrGrad_x =	OoqpVectorHandle ( tree->newPrimalVector()	); 
  trialBarrGrad_s =	OoqpVectorHandle ( tree->newDualZVector()	); 
  
  trialCeqBody  =	OoqpVectorHandle ( tree->newDualYVector()	);
  trialCIneqBody =	OoqpVectorHandle ( tree->newDualZVector()	);   

  dampind_xL_v  =	OoqpVectorHandle ( tree->newPrimalVector()	); 
  dampind_xU_w =	OoqpVectorHandle ( tree->newPrimalVector()	);  
  dampind_sL_t =	OoqpVectorHandle ( tree->newDualZVector()	);
  dampind_sU_u =	OoqpVectorHandle ( tree->newDualZVector()	);

  sc = OoqpVectorHandle ( tree->newPrimalVector() );

  nxlow = ixlow->numberOfNonzeros(); 
  nxupp = ixupp->numberOfNonzeros(); 
  mclow = iclow->numberOfNonzeros(); 
  mcupp = icupp->numberOfNonzeros();

  createChildren();
  global_nnz=-1;
  
}

sData::sData(int useMultiStage_, sTree* tree_, OoqpVector * c_in, SymMatrix * Q_in,
	     OoqpVector * xlow_in, OoqpVector * ixlow_in, long long nxlow_,
	     OoqpVector * xupp_in, OoqpVector * ixupp_in, long long nxupp_,
	     GenMatrix  * A_in, OoqpVector * bA_in,
	     GenMatrix  * C_in,
	     OoqpVector * clow_in, OoqpVector * iclow_in, long long mclow_,
	     OoqpVector * cupp_in, OoqpVector * icupp_in, long long mcupp_,
	     OoqpVector * CeqBody_in,OoqpVector * CIneqBody_in,
	     OoqpVector * trialBarrGrad_x_in,OoqpVector * trialBarrGrad_s_in,
	     OoqpVector * trialCeqBody, OoqpVector *trialCIneqBody,
	     OoqpVector * dampind_xL_v_in,OoqpVector * dampind_xU_w_in,
	     OoqpVector * dampind_sL_t_in, OoqpVector *dampind_sU_u_in)
  
  : NlpGenData(SparseLinearAlgebraPackage::soleInstance(),
	      c_in, Q_in, 
	      xlow_in, ixlow_in, nxlow_, 
	      xupp_in, ixupp_in, nxupp_,
	      A_in, bA_in,
	      C_in,
	      clow_in, iclow_in, mclow_,
	      cupp_in, icupp_in, mcupp_,
	      CeqBody_in, CIneqBody_in, 
	      trialBarrGrad_x_in, trialBarrGrad_s_in,
	      trialCeqBody, trialCIneqBody, 
	      dampind_xL_v_in, dampind_xU_w_in,
	      dampind_sL_t_in, dampind_sU_u_in)
{
  stochNode = tree_;

  useMultiStage=useMultiStage_;

  createChildren(useMultiStage);
  global_nnz=-1;

}

void sData::createChildren(int useMultiStage_)
{
  if(useMultiStage==0){
	createChildren();
	return;
  }
  assert("not yet"&&0); 
}


void sData::createChildren()
{
  //follow the structure of one of the tree objects and create the same
  //structure for this class, and link this object with the corresponding 
  //vectors and matrices
  StochVector& gSt     = dynamic_cast<StochVector&>(*grad);
  StochSymMatrix& QSt  = dynamic_cast<StochSymMatrix&>(*H); 
  StochVector& xlowSt  = dynamic_cast<StochVector&>(*blx); 
  StochVector& ixlowSt = dynamic_cast<StochVector&>(*ixlow); 
  StochVector& xuppSt  = dynamic_cast<StochVector&>(*bux); 
  StochVector& ixuppSt = dynamic_cast<StochVector&>(*ixupp);
  StochGenMatrix& ASt  = dynamic_cast<StochGenMatrix&>(*Jeq); 
  StochVector& bASt    = dynamic_cast<StochVector&>(*bA);
  StochGenMatrix& CSt  = dynamic_cast<StochGenMatrix&>(*Jineq);
  StochVector& clowSt  = dynamic_cast<StochVector&>(*bl); 
  StochVector& iclowSt = dynamic_cast<StochVector&>(*iclow);
  StochVector& cuppSt  = dynamic_cast<StochVector&>(*bu); 
  StochVector& icuppSt = dynamic_cast<StochVector&>(*icupp); 

  StochVector& CeqBodySt  = dynamic_cast<StochVector&>(*CeqBody); 
  StochVector& CIneqBodySt = dynamic_cast<StochVector&>(*CIneqBody);   


  StochVector& trialBarrGrad_xSt  = dynamic_cast<StochVector&>(*trialBarrGrad_x); 
  StochVector& trialBarrGrad_sSt = dynamic_cast<StochVector&>(*trialBarrGrad_s);   
  

  StochVector& trialCeqBodySt  = dynamic_cast<StochVector&>(*trialCeqBody); 
  StochVector& trialCIneqBodySt = dynamic_cast<StochVector&>(*trialCIneqBody);   

  StochVector& dampind_xL_vSt  = dynamic_cast<StochVector&>(*dampind_xL_v); 
  StochVector& dampind_xU_wSt = dynamic_cast<StochVector&>(*dampind_xU_w);   

  StochVector& dampind_sL_tSt  = dynamic_cast<StochVector&>(*dampind_sL_t); 
  StochVector& dampind_sU_uSt = dynamic_cast<StochVector&>(*dampind_sU_u); 


  for(size_t it=0; it<gSt.children.size(); it++) {
    AddChild(new sData(0,stochNode->children[it],
	       gSt.children[it], QSt.children[it],
	       xlowSt.children[it], ixlowSt.children[it], nxlow,
	       xuppSt.children[it], ixuppSt.children[it], nxupp,
	       ASt.children[it], bASt.children[it],
	       CSt.children[it],
	       clowSt.children[it], iclowSt.children[it], mclow,
	       cuppSt.children[it], icuppSt.children[it], mcupp,
	       CeqBodySt.children[it], CIneqBodySt.children[it],
	       trialBarrGrad_xSt.children[it], trialBarrGrad_sSt.children[it],
	       trialCeqBodySt.children[it], trialCIneqBodySt.children[it], 
	       dampind_xL_vSt.children[it], dampind_xU_wSt.children[it],
	       dampind_sL_tSt.children[it], dampind_sU_uSt.children[it]));
  }

}

void sData::destroyChildren()
{
  for(size_t it=0; it<children.size(); it++) {
    children[it]->destroyChildren();
    delete children[it];
  }
  children.clear();
}

void sData::AddChild(sData* child)
{
  children.push_back(child);
}

double sData::objectiveValue( NlpGenVars * vars )
{
  StochVector& x = dynamic_cast<StochVector&>(*vars->x);
  OoqpVectorHandle temp( x.clone() );

  this->getg( *temp );
  this->Qmult( 1.0, *temp, 0.5, *vars->x );

  return temp->dotProductWith( *vars->x );
}

void sData::createScaleFromQ()
{

  assert("Not implemented!" && 0);

  // Stuff the diagonal elements of H into the vector "sc"
  this->getDiagonalOfQ( *sc);
}



sData::~sData()
{
  for(size_t it=0; it<children.size(); it++)
    delete children[it];
}

int sData::getLocalnx()
{
  StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*H);
  return Qst.diag->size();
}

int sData::getLocalmy()
{
  long long my, nx;
  if(useMultiStage==0){
    StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*Jeq);
    Ast.Bmat->getSize(my, nx);
  }else{
	assert("not done" && 0);
  }
  return my;
}

int sData::getLocalmz()
{
  long long mz, nx;

  if(useMultiStage==0){
    StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*Jineq);
    Cst.Bmat->getSize(mz, nx);
  }else{
	assert("not done" && 0);
  }

  return mz;
}

int sData::getLocalSizes(int& nx, int& my, int& mz)
{
    long long nxll, myll, mzll;

  if(useMultiStage==0){
    StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*Jeq);
    Ast.Bmat->getSize(myll, nxll);  	
    StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*Jineq);
    Cst.Bmat->getSize(mzll, nxll);
  }else{
	assert("not done" && 0);
  }

  nx=nxll; my=myll; mz=mzll;
  return 0;
}


int sData::getLocalNnz(int& nnzQ, int& nnzB, int& nnzD)
{
  StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*H);
  //fixme: I think this is not correct!:
//  nnzQ = Qst.diag->getStorageRef().len + Qst.border->getStorageRef().len;
  nnzQ = Qst.diag->getStorageRef().len;

  if(useMultiStage==0){
    StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*Jeq);
    StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*Jineq);  	
    nnzB = Ast.Bmat->getStorageRef().len;
    nnzD = Cst.Bmat->getStorageRef().len;	
  }else{
	assert("not done" && 0);
  }
  
  return 0;
}


SparseSymMatrix& sData::getLocalQ()
{
  StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*H);
  return *Qst.diag;
}

SparseGenMatrix& sData::getLocalCrossHessian()
{
  StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*H);
  return *Qst.border;
}

// T_i x_0 + W_i x_i = b_i

// This is T_i
SparseGenMatrix& sData::getLocalA()
{
  if(useMultiStage==0){ 
  	StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*Jeq);
    return *Ast.Amat;
  }else{
	assert("not done" && 0);
  }
}

// This is W_i:
SparseGenMatrix& sData::getLocalB()
{
  if(useMultiStage==0){ 
  	StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*Jeq);
    return *Ast.Bmat;
  }else{
	assert("not done" && 0);
  }
}

// low_i <= C_i x_0 + D_i x_i <= upp_i

// This is C_i
SparseGenMatrix& sData::getLocalC()
{
  if(useMultiStage==0){ 
  	StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*Jineq);
    return *Cst.Amat;
  }else{
	assert("not done" && 0);
  }
}

// This is D_i
SparseGenMatrix& sData::getLocalD()
{
  if(useMultiStage==0){ 
  	StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*Jineq);
    return *Cst.Bmat;
  }else{
	assert("not done" && 0);
  }
}



void sData::sync()
{

  
  destroyChildren();
  
  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*grad)); 

//   int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
//   printf("vec -----------------------------------------------------\n");sleep(myRank+1);  
//   stochNode->displayVectorVsTreeStructure(dynamic_cast<StochVector&>(*grad), myRank);
//   printf("vec done ----------------------\n"); usleep(10000);

  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*blx));
  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*ixlow));
  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*bux));
  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*ixupp));
  stochNode->syncDualYVector(dynamic_cast<StochVector&>(*bA));
  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*bl));
  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*bu));
  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*iclow));
  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*icupp));

  stochNode->syncStochSymMatrix(dynamic_cast<StochSymMatrix&>(*H));

  if(useMultiStage==0){
    stochNode->syncStochGenMatrix(dynamic_cast<StochGenMatrix&>(*Jeq));
    stochNode->syncStochGenMatrix(dynamic_cast<StochGenMatrix&>(*Jineq));
  }else{
	assert("not done" && 0);
  }
  
  createChildren();
}




void sData::SetInputNlpPara(NlpInfo *updateNlp)
{
  long long dummy;
  
  updateNlp->nx = nx;
  updateNlp->my = my;
  updateNlp->mz = mz;

  updateNlp->A = Jeq;
  updateNlp->C = Jineq;
  updateNlp->Q = H;

  updateNlp->bA = bA;
  updateNlp->g = grad;

  nxLOri = updateNlp->nxL;
  nxUOri = updateNlp->nxU;
  nsLOri = updateNlp->nsL;
  nsUOri = updateNlp->nsU;

}







int sData::getGlobalNnz()
{ 
  if(global_nnz==-1) 
  	global_nnz = computeGlobalNnz();

  return global_nnz;
}  


// only support 2 stage problem. not for multi stage 
int sData::computeGlobalNnz()
{
  int totalNnz = 0;

  for (size_t it=0; it<children.size(); it++) {
  	totalNnz += children[it]->getGlobalNnz();
  } 
  
  if(stochNode->rankMe==0){
  	int tempNnz =0;
    MPI_Allreduce(&totalNnz, &tempNnz, 1, MPI_INT, MPI_SUM, stochNode->commWrkrs);
	totalNnz = tempNnz;
  }

  StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*H);
  totalNnz += Qst.diag->getStorageRef().len + Qst.border->getStorageRef().len;

  if(useMultiStage==0){
    StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*Jeq);
    StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*Jineq);  	
    totalNnz += Ast.Bmat->getStorageRef().len +  Ast.Amat->getStorageRef().len;
    totalNnz += Cst.Bmat->getStorageRef().len +  Cst.Amat->getStorageRef().len;	
  }else{
	assert("not done" && 0);
  }
  
  return totalNnz;
}

