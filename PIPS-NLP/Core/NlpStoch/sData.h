/* PIPS-IPM                                                             
 * Author: Cosmin G. Petra
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
 
/* 2015. Modified by Nai-Yuan Chiang for NLP*/


#ifndef DATANLPSTOCH
#define DATANLPSTOCH

#include "NlpGenData.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "DoubleMatrixHandle.h"
#include "sTree.h"


#include <vector>

//class sTree;
class LinearAlgebraPackage;
class NlpGenVars;
class NlpInfo;


class sData : public NlpGenData {
 public:
  std::string datarootname;
  std::string datalocalname;
  
  /** constructor that makes data objects of the specified dimensions */
  sData( sTree* tree);

  /** constructor that sets up pointers to the data objects that are
      passed as arguments */

  sData(int useMultiStage, sTree* tree_, OoqpVector * c_in, SymMatrix * Q_in,
	     OoqpVector * xlow_in, OoqpVector * ixlow_in, long long nxlow_,
	     OoqpVector * xupp_in, OoqpVector * ixupp_in, long long nxupp_,
	     GenMatrix  * A_in, OoqpVector * bA_in,
	     GenMatrix  * C_in,
	     OoqpVector * clow_in, OoqpVector * iclow_in, long long mclow_,
	     OoqpVector * cupp_in, OoqpVector * icupp_in, long long mcupp_,
	     OoqpVector * CeqBody_in,OoqpVector * CIneqBody_in,
	     OoqpVector * trialBarrGrad_x,OoqpVector * trialBarrGrad_s,
	     OoqpVector * trialCeqBody, OoqpVector *trialCIneqBody,
	     OoqpVector * dampind_xL_v_in,OoqpVector * dampind_xU_w_in,
	     OoqpVector * dampind_sL_t_in, OoqpVector *dampind_sU_u_in);

  
  std::vector<sData*> children;
  void AddChild(sData* child);
  sTree* stochNode;
 public:
//  long long nxlow, nxupp, mclow, mcupp;

  int getLocalnx();
  int getLocalmy();
  int getLocalmz();
  int getLocalSizes(int& nx, int& my, int& mz);

  int getLocalNnz(int& nnzQ, int& nnzB, int& nnzD);
  int getGlobalNnz();  

  SparseSymMatrix& getLocalQ();
  SparseGenMatrix& getLocalCrossHessian();
  SparseGenMatrix& getLocalA();
  SparseGenMatrix& getLocalB();
  SparseGenMatrix& getLocalC();
  SparseGenMatrix& getLocalD();

  void sync();
 public:


  virtual double objectiveValue( NlpGenVars * vars );
  virtual void createScaleFromQ();
  virtual void datainput() {};

  virtual ~sData();

  virtual void SetInputNlpPara(NlpInfo *updateNlp);


  virtual long long getGlobalNx(){return stochNode->N;};
  virtual long long getGlobalMy(){return stochNode->MY;};
  virtual long long getGlobalMz(){return stochNode->MZ;};


 protected:
  void createChildren();
  void createChildren(int useMultiStage);
  void destroyChildren();

  int useMultiStage;

  int global_nnz;

  // only support 2 stage problem. not for multi stage 
  int computeGlobalNnz();  
};


#endif
