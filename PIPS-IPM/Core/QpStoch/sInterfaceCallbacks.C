/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#include "sInterfaceCallbacks.h"
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


sInterfaceCallbacks::sInterfaceCallbacks(stochasticInput &in)
  : inputTree(NULL)
{
  callbackData=new CallbackData(in);
  loadData();
}

sInterfaceCallbacks::~sInterfaceCallbacks()
{
  if(callbackData)
    delete callbackData;
  if(inputTree)
    delete inputTree;
}

void sInterfaceCallbacks::loadData()
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

void sInterfaceCallbacks::go()
{
  double t = MPI_Wtime();
  StochRunParams* params = defaultStochRunParams();

  MehrotraStochSolver* method=NULL;
  sFactoryAug* formulation=NULL;

  qpstoch_solve(inputTree,params,method,formulation);

  t = MPI_Wtime() - t;
}

void sInterfaceCallbacks::getNum1stStgEqIneq(int& my, int &mz)
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

void sInterfaceCallbacks::getNum2ndStgEqIneq(int scen, int& my, int &mz)
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

namespace {

class eq_comp
{
public:
  inline bool   operator()(const double& lb, const double& ub) const { return (lb==ub); }
};
class ineq_comp
{
public:
  inline bool   operator()(const double& lb, const double& ub) const { return (lb!=ub); }
};

/** Counts the nnz in Mcol's rows corresponding to entries in lb and ub 
 * satisfying 'compFun' condition. */
template<typename Compare>
int countNNZ(const CoinPackedMatrix& Mcol, 
	     const vector<double>& lb, const vector<double>& ub, 
	     const Compare& compFun)
{
  int nnz=0;

  //convert to row-major, probably getting rows is faster 
  CoinPackedMatrix M; M.reverseOrderedCopyOf(Mcol); 

  assert(false==M.hasGaps());  

  size_t R=lb.size();
  for(size_t i=0; i<R; i++) {
    if (compFun(lb[i],ub[i])) {
      nnz += M.getVectorSize(i);
    }
  }
  return nnz;
}		 

/** Extracts the in Mcol's rows corresponding to entries in lb and ub 
 * satisfying 'compFun' condition. Mcol is in column-major format, the
 * output krowM,jcolM,dM represent a row-major submatrix of Mcol. */
template<typename Compare>
void extractRows(const CoinPackedMatrix& Mcol, 
		 const vector<double>& lb, const vector<double>& ub, 
		 const Compare& compFun,
		 int* krowM, int* jcolM, double* dM)
{
  //convert to row-major, extracting rows is probably faster 
  CoinPackedMatrix M; M.reverseOrderedCopyOf(Mcol); 

  size_t R=lb.size();

  int nRow=0;
  int* indRow=new int[R];//overallocated, but that's fine
  
  for(size_t i=0; i<R; i++) {
    if (compFun(lb[i],ub[i])) {
      indRow[nRow++]=i;
    }
  }

  if (nRow==0) {
    krowM[0]=0;
    delete[] indRow;
    return;
  }

  CoinPackedMatrix Msub;

  Msub.submatrixOf(M, nRow, indRow); //this seems to crash if nRow==0

  assert(false==Msub.hasGaps());

  memcpy(krowM,Msub.getVectorStarts(),(Msub.getNumRows()+1)*sizeof(int));
  memcpy(jcolM,Msub.getIndices(),Msub.getNumElements()*sizeof(int));
  memcpy(dM,Msub.getElements(),Msub.getNumElements()*sizeof(double));

  /*printf("Rows=%d Cols=%d nnz=%d\n", Msub.getNumRows(), Msub.getNumCols(), Msub.getNumElements());
  for(int r=0;r<=Msub.getNumRows(); r++) printf("%d ", krowM[r]); printf("\n");
  for(int e=0; e<Msub.getNumElements(); e++) printf("%d ", jcolM[e]); printf("\n");
  for(int e=0; e<Msub.getNumElements(); e++) printf("%12.5e ", dM[e]); printf("\n");
  printf("---------------\n");
  */

  //r.getElements() returns a vector containing:
  //    3 1 -2 -1 -1 2 1.1 1 1 2.8 -1.2 5.6 1 1.9
  //  r.getIndices() returns a vector containing:
  //    0 1  3  4  7 1 2   2 5 3    6   0   4 7
  //  r.getVectorStarts() returns a vector containing:
  //    0 5 7 9 11 14
}		 

}

extern "C" int fnnzA(void* user_data, int id, int* nnz){
  eq_comp compFun;
  if(id==0) {
    CoinPackedMatrix Mcol=callbackData->in.getFirstStageConstraints();
    cout << "StochInput constraints: nnz=" << Mcol.getNumElements() << endl;
    
    (*nnz)=countNNZ(Mcol,
		    callbackData->in.getFirstStageRowLB(), 
		    callbackData->in.getFirstStageRowUB(), 
		    compFun);
  } else {
    int scen=id-1;
    CoinPackedMatrix Mcol=callbackData->in.getLinkingConstraints(scen);

    (*nnz)=countNNZ(Mcol,
		    callbackData->in.getSecondStageRowLB(scen), 
		    callbackData->in.getSecondStageRowUB(scen), 
		    compFun);
  }
  printf("nnzA = %d\n", *nnz);
  return 0;
}

extern "C" int fnnzB(void* user_data, int id, int* nnz){

  assert(id>=1);    
  int scen=id-1;
  CoinPackedMatrix Mcol=callbackData->in.getSecondStageConstraints(scen);
  eq_comp compFun;
  (*nnz)=countNNZ(Mcol,
		  callbackData->in.getSecondStageRowLB(scen),
		  callbackData->in.getSecondStageRowUB(scen),
		  compFun);
  printf("nnzB = %d\n", *nnz);

  return 0;
}

extern "C" int fnnzC(void* user_data, int id, int* nnz){

  if(id==0) {
    CoinPackedMatrix Mcol=callbackData->in.getFirstStageConstraints();
    (*nnz)=countNNZ(Mcol,
		    callbackData->in.getFirstStageRowLB(), 
		    callbackData->in.getFirstStageRowUB(), 
		    ineq_comp());
  } else {
    int scen=id-1;
    CoinPackedMatrix Mcol=callbackData->in.getLinkingConstraints(scen);

    (*nnz)=countNNZ(Mcol,
		    callbackData->in.getSecondStageRowLB(scen), 
		    callbackData->in.getSecondStageRowUB(scen), 
		    ineq_comp());
  }
  //printf("nnzC = %d\n", *nnz);
  return 0;
}

extern "C" int fnnzD(void* user_data, int id, int* nnz){
  assert(id>=1);    
  int scen=id-1;
  CoinPackedMatrix Mcol=callbackData->in.getSecondStageConstraints(scen);

  (*nnz)=countNNZ(Mcol,
		  callbackData->in.getSecondStageRowLB(scen),
		  callbackData->in.getSecondStageRowUB(scen),
		  ineq_comp());
  //printf("nnzD = %d\n", *nnz);

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
  //printf("fA %d\n",id);
  if(id==0) {
    CoinPackedMatrix Mcol=callbackData->in.getFirstStageConstraints();
    extractRows(Mcol,
		callbackData->in.getFirstStageRowLB(), 
		callbackData->in.getFirstStageRowUB(), 
		eq_comp(),
		krowM, jcolM, M);

  } else {
    int scen=id-1;
    CoinPackedMatrix Mcol=callbackData->in.getLinkingConstraints(scen);

    extractRows(Mcol,
		callbackData->in.getSecondStageRowLB(scen), 
		callbackData->in.getSecondStageRowUB(scen), 
		eq_comp(),
		krowM, jcolM, M);
  }  
  return 0;
}

extern "C" int fB(void* user_data, int id, int* krowM, int* jcolM, double* M){
  //printf("fB %d\n",id);
  assert(id>=1);
  CoinPackedMatrix Mcol=callbackData->in.getSecondStageConstraints(id-1);
  extractRows(Mcol,
	      callbackData->in.getSecondStageRowLB(id-1),
	      callbackData->in.getSecondStageRowUB(id-1),
	      eq_comp(),
	      krowM, jcolM, M);
  return 0;
}

extern "C" int fC(void* user_data, int id, int* krowM, int* jcolM, double* M){
  //printf("fC %d\n",id);
  if(id==0) {
    CoinPackedMatrix Mcol=callbackData->in.getFirstStageConstraints();
    extractRows(Mcol,
		callbackData->in.getFirstStageRowLB(), 
		callbackData->in.getFirstStageRowUB(), 
		ineq_comp(),
		krowM, jcolM, M);

  } else {
    int scen=id-1;
    CoinPackedMatrix Mcol=callbackData->in.getLinkingConstraints(scen);

    extractRows(Mcol,
		callbackData->in.getSecondStageRowLB(scen), 
		callbackData->in.getSecondStageRowUB(scen), 
		ineq_comp(),
		krowM, jcolM, M);
  }
  return 0;
}

extern "C" int fD(void* user_data, int id, int* krowM, int* jcolM, double* M){
  assert(id>=1);

  CoinPackedMatrix Mcol=callbackData->in.getSecondStageConstraints(id-1);
  extractRows(Mcol,
	      callbackData->in.getSecondStageRowLB(id-1),
	      callbackData->in.getSecondStageRowUB(id-1),
	      ineq_comp(),
	      krowM, jcolM, M);
  return 0;
}

extern "C" int fc(void* user_data, int id, double* vec, int len){

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
  vector<double> lb, ub;

  if(id==0) {
    lb = callbackData->in.getFirstStageRowLB();
    ub = callbackData->in.getFirstStageRowUB();
  }  else {
    lb = callbackData->in.getSecondStageRowLB(id-1);
    ub = callbackData->in.getSecondStageRowUB(id-1);
  }
  int eq_cnt=0;
  for(size_t i=0; i<lb.size(); i++)
    if(lb[i]==ub[i])
      vec[eq_cnt++]=lb[i];  
  assert(eq_cnt==len);

  return 0;
}

extern "C" int fclow  (void* user_data, int id, double* vec, int len){
  vector<double> lb, ub;

  if(id==0) {
    lb = callbackData->in.getFirstStageRowLB();
    ub = callbackData->in.getFirstStageRowUB();
  }  else {
    lb = callbackData->in.getSecondStageRowLB(id-1);
    ub = callbackData->in.getSecondStageRowUB(id-1);
  }
  int ineq_cnt=0;
  for(size_t i=0; i<lb.size(); i++)
    if(lb[i]!=ub[i]) {
      if(lb[i]>-1e20)
	vec[ineq_cnt]=lb[i];  
      else
	vec[ineq_cnt]=0.0;
      ineq_cnt++;
    }
  assert(ineq_cnt==len);

  return 0;
}

extern "C" int ficlow (void* user_data, int id, double* vec, int len){
  vector<double> lb, ub;

  if(id==0) {
    lb = callbackData->in.getFirstStageRowLB();
    ub = callbackData->in.getFirstStageRowUB();
  }  else {
    lb = callbackData->in.getSecondStageRowLB(id-1);
    ub = callbackData->in.getSecondStageRowUB(id-1);
  }
  int ineq_cnt=0;
  for(size_t i=0; i<lb.size(); i++)
    if(lb[i]!=ub[i]) {
      if(lb[i]>-1e20)
	vec[ineq_cnt]=1.0;
      else
	vec[ineq_cnt]=0.0;
      ineq_cnt++;
    }
  assert(ineq_cnt==len);
  return 0;
}

extern "C" int fcupp  (void* user_data, int id, double* vec, int len){
  vector<double> lb, ub;

  if(id==0) {
    lb = callbackData->in.getFirstStageRowLB();
    ub = callbackData->in.getFirstStageRowUB();
  }  else {
    lb = callbackData->in.getSecondStageRowLB(id-1);
    ub = callbackData->in.getSecondStageRowUB(id-1);
  }
  int ineq_cnt=0;
  for(size_t i=0; i<lb.size(); i++)
    if(lb[i]!=ub[i]) {
      if(ub[i]<1e20)
	vec[ineq_cnt]=ub[i];  
      else
	vec[ineq_cnt]=0.0;
      ineq_cnt++;
    }
  return 0;
}

extern "C" int ficupp (void* user_data, int id, double* vec, int len){
  vector<double> lb, ub;

  if(id==0) {
    lb = callbackData->in.getFirstStageRowLB();
    ub = callbackData->in.getFirstStageRowUB();
  }  else {
    lb = callbackData->in.getSecondStageRowLB(id-1);
    ub = callbackData->in.getSecondStageRowUB(id-1);
  }
  int ineq_cnt=0;
  for(size_t i=0; i<lb.size(); i++)
    if(lb[i]!=ub[i]) {
      if(ub[i]<1e20)
	vec[ineq_cnt]=1.0;
      else
	vec[ineq_cnt]=0.0;
      ineq_cnt++;
    }
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

