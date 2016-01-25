/* PIPS-NLP                                                             
 * Author: Nai-Yuan Chiang
 * (C) 2015 Argonne National Laboratory
 */

#ifndef NLPINFOFROMNL
#define NLPINFOFROMNL


#include "sInfo.h"
#include "OoqpVectorHandle.h"
#include "OoqpVector.h"

#include "AmplData_NL.hpp"

class NlpGenVars;
class sData;
class stochasticInput;
class amplGenStochInput;
class amplGenStochInput_AddSlack;


struct ASL_pfgh;

class sNlpInfoFromNL : public sInfo
{

protected: 

  std::string datarootname;

  ASL_pfgh* asl_local;  

  double ObjScal;
  
  int iAmDistrib;


  int *amplRowMap;

  std::map<int,int> *LocGloVarMap;
  std::map<int,int> *LocLocVarMap;  
  




  std::map<int,int> *LocWmatJacGoffMap;	
  std::map<int,int> *LocTmatJacGoffMap;	

  std::map<int,int> *LocQAmatHesGoffMap;
  std::map<int,int> *LocQWmatHesGoffMap;
  std::map<int,int> *LocQTmatHesGoffMap; 

  int nnzQDiag, nnzQCross, nnzQParent;

  std::map<int,int> *LocAeqLinkJacGoffMap;	
  std::map<int,int> *LocBeqLocJacGoffMap;	
  std::map<int,int> *LocCineqLinkJacGoffMap;	
  std::map<int,int> *LocDineqLocJacGoffMap;	

  int nnzAeqLink, nnzCineqLink, nnzBeqLoc, nnzDineqLoc;
  
public:


  
	
  sNlpInfoFromNL();    
  virtual ~sNlpInfoFromNL();

  
  sNlpInfoFromNL(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in);
  sNlpInfoFromNL( int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in,
				  int nxL_in,int nxU_in,int nsL_in,int nsU_in);  

  sNlpInfoFromNL(sData *data_in){}
  sNlpInfoFromNL(sData *data_in, stochasticInput& in);
  sNlpInfoFromNL(sData *data_in, amplGenStochInput& in, const int ChildIdx);
  sNlpInfoFromNL(sData *data_in, amplGenStochInput_AddSlack& in, const int ChildIdx);

  void createChildren(sData *data_in, amplGenStochInput& in);
  void createChildren(sData *data_in, amplGenStochInput_AddSlack& in);
  
  virtual double ObjValue( NlpGenVars * vars) ;
  
  virtual void ConstraintBody( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq);

  virtual int ObjGrad( NlpGenVars * vars, OoqpVector *grad );

  

  virtual void Hessian( NlpGenVars * vars, SymMatrix *Hess );

  virtual void JacFull( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC); 

  virtual void get_InitX0(OoqpVector* vX);

  virtual void Hessian_FromSon( NlpGenVars * vars, double *tempFromParH );

  virtual void ObjGrad_FromSon( NlpGenVars * vars, OoqpVector *grad,double *tempFromParGrad );

  virtual void writeSolution( NlpGenVars * vars_);
private:
  double ObjValue_General( NlpGenVars * vars_);
  double ObjValue_DummyCon( NlpGenVars * vars_);
  
  void ConstraintBody_General( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq);
  void ConstraintBody_DummyCon( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq);

  void ObjGrad_General( NlpGenVars * vars, OoqpVector *grad );
  void ObjGrad_DummyCon( NlpGenVars * vars, OoqpVector *grad );

  void Hessian_General( NlpGenVars * vars, SymMatrix *Hess );
  void Hessian_DummyCon( NlpGenVars * vars, SymMatrix *Hess );  

  void JacFull_General( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC); 
  void JacFull_DummyCon( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC); 

  void get_InitX0_General(OoqpVector* vX);
  void get_InitX0_DummyCon(OoqpVector* vX);
 
};




#endif

