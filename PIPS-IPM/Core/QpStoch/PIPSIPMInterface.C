#include "PIPSIPMInterface.h"
#include "MehrotraStochSolver.h"
#include "sFactoryAug.h"
#include "mpi.h"

class CallbackData
{
public:
  CallbackData(stochasticInput &in_): in(in_){};
  stochasticInput &in;
};

CallbackData* callbackData=NULL;



extern "C" {
int fnnzQ(void* user_data, int id, int* nnz);
int fnnzA(void* user_data, int id, int* nnz);
int fnnzB(void* user_data, int id, int* nnz);
int fnnzC(void* user_data, int id, int* nnz);
int fnnzD(void* user_data, int id, int* nnz);
int fQ(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fA(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fB(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fC(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fD(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fc(void* user_data, int id, double* vec, int len);
int fb(void* user_data, int id, double* vec, int len);
int fclow  (void* user_data, int id, double* vec, int len);
int ficlow (void* user_data, int id, double* vec, int len);
int fcupp  (void* user_data, int id, double* vec, int len);
int ficupp (void* user_data, int id, double* vec, int len);
int fxlow  (void* user_data, int id, double* vec, int len);
int fixlow (void* user_data, int id, double* vec, int len);
int fxupp  (void* user_data, int id, double* vec, int len);
int fixupp (void* user_data, int id, double* vec, int len);
}


PIPSIPMInterface::PIPSIPMInterface(stochasticInput &in)
  : inputTree(NULL)
{
  callbackData=new CallbackData(in);
}

PIPSIPMInterface::~PIPSIPMInterface()
{
  if(callbackData)
    delete callbackData;
  if(inputTree)
    delete inputTree;
}

void PIPSIPMInterface::loadData()
{  
  int nx0,my0,mz0;
  nx0=callbackData->in.nFirstStageVars();
  getNum1stStgEqIneq(my0,mz0);

  StochInputTree::
    StochInputNode data(&callbackData, 0, 
			nx0, my0, mz0,
			fQ, fnnzQ, fc,
			fA, fnnzA,
			fB,  fnnzB,
			fb,
			fC, fnnzC,
			fD, fnnzD,
			fclow, ficlow, fcupp, ficupp,
			fxlow, fixlow, fxupp, fixupp );
  inputTree = new StochInputTree(data);
  
  
  for(int i=0; i<callbackData->in.nScenarios(); i++) {

    int my,mz;
    getNum2ndStgEqIneq(i,my,mz);

    StochInputTree::StochInputNode data(&callbackData, i+1,
			     callbackData->in.nSecondStageVars(i), my, mz, 
			     fQ, fnnzQ, fc,
			     fA, fnnzA,
			     fB, fnnzB,
			     fb,
			     fC, fnnzC,
			     fD, fnnzD,
			     fclow, ficlow, fcupp, ficupp,
			     fxlow, fixlow, fxupp, fixupp );
    
    inputTree->AddChild(new StochInputTree(data));
  }
}

void PIPSIPMInterface::go()
{
  double t = MPI_Wtime();
  StochRunParams* params = defaultStochRunParams();

  MehrotraStochSolver* method=NULL;
  sFactoryAug* formulation=NULL;

  qpstoch_solve(inputTree,params,method,formulation);

  t = MPI_Wtime() - t;
}

void PIPSIPMInterface::getNum1stStgEqIneq(int& my, int &mz)
{
  vector<double> lb=callbackData->in.getFirstStageRowLB();
  vector<double> ub=callbackData->in.getFirstStageRowUB();
  my=0; mz=0;

  for(size_t i=0;i<lb.size(); i++) {
    if(lb[i]==ub[i]) my++;
    else if(lb[i]>-1e+20 || ub[i]<1e+20) //should be always true
      mz++;
    else assert(false);
  }
  //printf("1stStg:%d %d\n", my,mz);
}

void PIPSIPMInterface::getNum2ndStgEqIneq(int scen, int& my, int &mz)
{
  vector<double> lb=callbackData->in.getSecondStageRowLB(scen);
  vector<double> ub=callbackData->in.getSecondStageRowUB(scen);
  my=0; mz=0;

  for(size_t i=0;i<lb.size(); i++) {
    if(lb[i]==ub[i]) my++;
    else if(lb[i]>-1e+20 || ub[i]<1e+20) //should be always true
      mz++;
    else assert(false);
  }
  //printf("2ndStg:%d %d\n", my,mz);
}


extern "C" int fnnzQ(void* user_data, int id, int* nnz)
{
  //CallbackData* data=(CallbackData*)user_data;
  (*nnz)=0;
  return 0;
}

extern "C" int fnnzA(void* user_data, int id, int* nnz){
  if(id==0) {
    CoinPackedMatrix A=callbackData->in.getFirstStageConstraints();
  } else {
    CoinPackedMatrix W=callbackData->in.getLinkingConstraints(id-1);
    // W.getVectorSize(0); nnz on vector (row or col) 0
    // ooqp  row-major
  }
  (*nnz)=0;
  return 0;
}

extern "C" int fnnzB(void* user_data, int id, int* nnz){
  CallbackData* data=(CallbackData*)user_data;
  printf("%d  ", id);
  return 0;
}

extern "C" int fnnzC(void* user_data, int id, int* nnz){
  CallbackData* data=(CallbackData*)user_data;

  return 0;
}

extern "C" int fnnzD(void* user_data, int id, int* nnz){
  CallbackData* data=(CallbackData*)user_data;

  return 0;
}

extern "C" int fQ(void* user_data, int id, int* krowM, int* jcolM, double* M){
  int n;
  if(id==0)
    n=callbackData->in.nFirstStageVars();
  else
    n=callbackData->in.nSecondStageVars(id-1);

  for(int i=0;i<n;i++)
    krowM[i]=0;

  return 0;
}

extern "C" int fA(void* user_data, int id, int* krowM, int* jcolM, double* M){

  if(id==0) {
    CoinPackedMatrix A=callbackData->in.getFirstStageConstraints();
  } else {
    CoinPackedMatrix W=callbackData->in.getLinkingConstraints(id-1);
  }
  
  return 0;
}

extern "C" int fB(void* user_data, int id, int* krowM, int* jcolM, double* M){

  return 0;
}

extern "C" int fC(void* user_data, int id, int* krowM, int* jcolM, double* M){
  CallbackData* data=(CallbackData*)user_data;

  return 0;
}

extern "C" int fD(void* user_data, int id, int* krowM, int* jcolM, double* M){
  CallbackData* data=(CallbackData*)user_data;

  return 0;
}

extern "C" int fc(void* user_data, int id, double* vec, int len){

  CallbackData* data=(CallbackData*) user_data;


  if(id==0) {
    vector<double> c = callbackData->in.getFirstStageObj();
    copy(c.begin(), c.end(), vec);
  }  else {
    vector<double> c = callbackData->in.getSecondStageObj(id-1);
    copy(c.begin(), c.end(), vec);
  }
  return 0;
}

extern "C" int fb(void* user_data, int id, double* vec, int len){
  CallbackData* data=(CallbackData*)user_data;

  return 0;
}

extern "C" int fclow  (void* user_data, int id, double* vec, int len){
  CallbackData* data=(CallbackData*)user_data;

  return 0;
}

extern "C" int ficlow (void* user_data, int id, double* vec, int len){
  CallbackData* data=(CallbackData*)user_data;

  return 0;
}

extern "C" int fcupp  (void* user_data, int id, double* vec, int len){
  CallbackData* data=(CallbackData*)user_data;

  return 0;
}

extern "C" int ficupp (void* user_data, int id, double* vec, int len){
  CallbackData* data=(CallbackData*)user_data;

  return 0;
}

extern "C" int fxlow  (void* user_data, int id, double* vec, int len){
  vector<double> x;
  if(id==0) {
    x=callbackData->in.getFirstStageColLB();
  } else {
    x=callbackData->in.getSecondStageColLB(id-1);
  }
  assert(len-x.size()==0);
  for(size_t i=0; i<x.size(); i++)
    if(x[i]>-1e+20) 
      vec[i]=x[i];
    else 
      vec[i]=0.0;
  
  return 0;
}

extern "C" int fixlow (void* user_data, int id, double* vec, int len){
  vector<double> x;
  if(id==0) {
    x=callbackData->in.getFirstStageColLB();
  } else {
    x=callbackData->in.getSecondStageColLB(id-1);
  }
  assert(len-x.size()==0);
  for(size_t i=0; i<x.size(); i++)
    if(x[i]>-1e+20) 
      vec[i]=1.0;
    else 
      vec[i]=0.0;

  return 0;
}

extern "C" int fxupp  (void* user_data, int id, double* vec, int len){
  vector<double> x;
  if(id==0) {
    x=callbackData->in.getFirstStageColUB();
  } else {
    x=callbackData->in.getSecondStageColUB(id-1);
  }

  assert(len-x.size()==0);
  for(size_t i=0; i<x.size(); i++)
    if(x[i]<1e+20) 
      vec[i]=x[i];
    else 
      vec[i]=0.0;


  return 0;
}

extern "C" int fixupp (void* user_data, int id, double* vec, int len){
  vector<double> x;
  if(id==0) {
    x=callbackData->in.getFirstStageColUB();
  } else {
    x=callbackData->in.getSecondStageColUB(id-1);
  }

  assert(len-x.size()==0);
  for(size_t i=0; i<x.size(); i++)
    if(x[i]<1e+20) 
      vec[i]=1.0;
    else 
      vec[i]=0.0;

  return 0;
}

