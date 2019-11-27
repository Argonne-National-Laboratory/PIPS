#include "StochInputTree.h"
#include "QpGenStochDriver.h"
#include "MehrotraSolver.h"
#include "mpi.h"
#include "pipsport.h"

#include <math.h>
#include <vector>
#include <assert.h>
#include <cstdlib>

#include <string>
#include <stdio.h>
#include <iostream>
using namespace std;

extern "C" {
/**** The dynamically generated example. */
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
int fclow (void* user_data, int id, double* vec, int len);
int ficlow(void* user_data, int id, double* vec, int len);
int fcupp (void* user_data, int id, double* vec, int len);
int ficupp(void* user_data, int id, double* vec, int len);
int fxlow (void* user_data, int id, double* vec, int len);
int fixlow(void* user_data, int id, double* vec, int len);
int fxupp (void* user_data, int id, double* vec, int len);
int fixupp(void* user_data, int id, double* vec, int len);

/***** The callbacks for 2 level fixed example */
int fnnzQ0(void* user_data, int id, int* nnz);
int fnnzQ1(void* user_data, int id, int* nnz);
int fnnzQ2(void* user_data, int id, int* nnz);

int fnnzA0(void* user_data, int id, int* nnz);
int fnnzA1(void* user_data, int id, int* nnz);
int fnnzA2(void* user_data, int id, int* nnz);

int fnnzB0(void* user_data, int id, int* nnz);
int fnnzB1(void* user_data, int id, int* nnz);
int fnnzB2(void* user_data, int id, int* nnz);

int fnnzC0(void* user_data, int id, int* nnz);
int fnnzC1(void* user_data, int id, int* nnz);
int fnnzC2(void* user_data, int id, int* nnz);

int fnnzD0(void* user_data, int id, int* nnz);
int fnnzD1(void* user_data, int id, int* nnz);
int fnnzD2(void* user_data, int id, int* nnz);

int fQ0(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fQ1(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fQ2(void* user_data, int id, int* krowM, int* jcolM, double* M);

int fA0(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fA1(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fA2(void* user_data, int id, int* krowM, int* jcolM, double* M);

int fB0(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fB1(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fB2(void* user_data, int id, int* krowM, int* jcolM, double* M);

int fC0(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fC1(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fC2(void* user_data, int id, int* krowM, int* jcolM, double* M);

int fD0(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fD1(void* user_data, int id, int* krowM, int* jcolM, double* M);
int fD2(void* user_data, int id, int* krowM, int* jcolM, double* M);

int fc0(void* user_data, int id, double* vec, int len);
int fb0(void* user_data, int id, double* vec, int len);

int fclow0 (void* user_data, int id, double* vec, int len);
int ficlow0(void* user_data, int id, double* vec, int len);
int fcupp0 (void* user_data, int id, double* vec, int len);
int ficupp0(void* user_data, int id, double* vec, int len);

int fxlow0 (void* user_data, int id, double* vec, int len);
int fixlow0(void* user_data, int id, double* vec, int len);
int fxupp0 (void* user_data, int id, double* vec, int len);
int fixupp0(void* user_data, int id, double* vec, int len);

int fc1(void* user_data, int id, double* vec, int len);
int fb1(void* user_data, int id, double* vec, int len);

int fclow1 (void* user_data, int id, double* vec, int len);
int ficlow1(void* user_data, int id, double* vec, int len);
int fcupp1 (void* user_data, int id, double* vec, int len);
int ficupp1(void* user_data, int id, double* vec, int len);

int fxlow1 (void* user_data, int id, double* vec, int len);
int fixlow1(void* user_data, int id, double* vec, int len);
int fxupp1 (void* user_data, int id, double* vec, int len);
int fixupp1(void* user_data, int id, double* vec, int len);

int fc2(void* user_data, int id, double* vec, int len);
int fb2(void* user_data, int id, double* vec, int len);

int fclow2 (void* user_data, int id, double* vec, int len);
int ficlow2(void* user_data, int id, double* vec, int len);
int fcupp2 (void* user_data, int id, double* vec, int len);
int ficupp2(void* user_data, int id, double* vec, int len);

int fxlow2 (void* user_data, int id, double* vec, int len);
int fixlow2(void* user_data, int id, double* vec, int len);
int fxupp2 (void* user_data, int id, double* vec, int len);
int fixupp2(void* user_data, int id, double* vec, int len);
}

struct MyNodeData
{ 
  int level;
  int id;
  int n;
};

StochInputTree* buildStochNode(const int& globalID,
			  const int& level, const int& noSubnodes, const int& nOfEachLevel,
			  vector<MyNodeData*>& vecMyNodeData)
{
  MyNodeData* mydata = new MyNodeData;
  mydata->level=level;
  mydata->id = globalID;
  mydata->n = nOfEachLevel;

  vecMyNodeData.push_back(mydata);

  StochInputTree::StochInputNode data( mydata, globalID, 
			    nOfEachLevel, nOfEachLevel/5, nOfEachLevel/2,
			    fQ, fnnzQ, fc,
			    fA, fnnzA,
			    fB,  fnnzB,
			    fb,
			    fC, fnnzC,
			    fD, fnnzD,
			    fclow, ficlow, fcupp, ficupp,
			    fxlow, fixlow, fxupp, fixupp );

  StochInputTree* node = new StochInputTree(data);

  //!log
  /*string strIdent = "";
  for(int it=0; it<level; it++)
    strIdent = strIdent + "  ";
  strIdent = strIdent + "- ";
  cout << strIdent << "node with id " << globalID << " allocated." << endl;
  */
  return node;
}

void generateStochTree(int level, 
		       const int& levels, const int& noSubnodes, const int& nOfEachLevel, 
		       StochInputTree* root, int& globalID,
		       vector<MyNodeData*>& vecMyNodeData)
{

  for(int i=0; i<noSubnodes; i++) {
    StochInputTree* node = buildStochNode(globalID++,level,noSubnodes,nOfEachLevel,vecMyNodeData);

    if(level<levels)
      generateStochTree(level+1, levels, noSubnodes, nOfEachLevel, node, globalID, vecMyNodeData);

    root->AddChild(node);    
  }
}

int generateAndSolveArtificalProblem( int argc, char* argv[], 
				      int levels, int noSubnodes, int nOfEachLevel)
{

  MPI_Init(&argc, &argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int size; MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(0==rank) {
    //cout << "---------------------------------------------------------------------------\n";
    cout << "Generating the problem .... " << endl;
  }

  vector<MyNodeData*> vecMyNodeData;

  //build the root
  int level = 0;
  StochInputTree* root = buildStochNode(0,level,noSubnodes,nOfEachLevel,vecMyNodeData);

  level=1;
  int globalId = 1;
  generateStochTree(level, levels, noSubnodes, nOfEachLevel, root, globalId, vecMyNodeData);

  if(0==rank) {
    cout << "=======================" << endl;
    cout << "Stochastic test problem" << endl;
    cout << "=======================" << endl;

    cout << "  - A problem with " << levels << " levels "
       << " and " << noSubnodes << " subnodes per node was generated." << endl;

    cout << "  - Each node has " << nOfEachLevel << " variables." << endl;

    cout << "  - The stoch tree has a total of " 
	 << (pow(1.0*noSubnodes, levels+1.0)-1)/(noSubnodes-1) << " nodes." << endl;

    cout << "  - The problem has a total of " 
	 << nOfEachLevel*(pow(1.0*noSubnodes, levels+1.0)-1)/(noSubnodes-1) 
	 << " variables." << endl;
    //cout << "-----------------------------------------------------------------------------\n";
    cout << endl << "Solving with " << size << " processes... " << endl;
  }
#ifdef TESTINGG
  root->runTestNGP();
  return 0;
#endif
  StochRunParams* params = defaultStochRunParams();
  //  params->printx = 1;
 
  MehrotraSolver* method=nullptr;
  qpgenstoch_solve( root, method, params);
  
  delete params;  
  
  for(int i=0; i<vecMyNodeData.size(); i++) delete vecMyNodeData[i];
  vecMyNodeData.clear();

  MPI_Finalize();

  return 0;
}

int main( int argc, char *argv[] )
{ 
  if(argc>1) {
    int levels = atoi(argv[1]);
    
    int noSubnodes = 4;
    if(argc>2) noSubnodes = atoi(argv[2]);

    int n=50;
    if(argc>3) n = atoi(argv[3]);
    
    generateAndSolveArtificalProblem( argc, argv, levels, noSubnodes, n);

  } else {
    void* usrdata=nullptr;
    StochInputTree::StochInputNode node( usrdata, 0, 
			      5, 2, 2,
			      fQ0, fnnzQ0, fc0,
			      fA0, fnnzA0,
			      fB0,  fnnzB0,
			      fb0,
			      fC0, fnnzC0,
			      fD0, fnnzD0,
			      fclow0, ficlow0, fcupp0, ficupp0,
			      fxlow0, fixlow0, fxupp0, fixupp0 );
    
    StochInputTree::StochInputNode node1( usrdata, 1, 
			     2, 1, 1,
			     fQ1, fnnzQ1, fc1,
			     fA1, fnnzA1,
			     fB1, fnnzB1,
			     fb1,
			     fC1, fnnzC1,
			     fD1, fnnzD1,
			     fclow1, ficlow1, fcupp1, ficupp1,
			     fxlow1, fixlow1, fxupp1, fixupp1 );
  
    StochInputTree::StochInputNode node2( usrdata, 2, 
			     2, 1, 2,
			     fQ2, fnnzQ2, fc2,
			     fA2, fnnzA2,
			     fB2, fnnzB2,
			     fb2,
			     fC2, fnnzC2,
			     fD2, fnnzD2,
			     fclow2, ficlow2, fcupp2, ficupp2,
			     fxlow2, fixlow2, fxupp2, fixupp2 );
  
    StochInputTree root(node);
    root.AddChild(node1);
    root.AddChild(node2);
    
    StochRunParams* params = defaultStochRunParams();
    
    MehrotraSolver* method=nullptr;
    qpgenstoch_solve( argc, argv, &root, method, params);
    
    delete params;
  }
  return 0;
}



/**************  CALBACKS   *********************/

extern "C"
int fQ0(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==0);
  int n=5; int nnz=7;
  
//   krowM[0]=0; krowM[1]=2;krowM[2]=4;krowM[3]=5;krowM[4]=6;krowM[5]=7;
  
//   jcolM[0]=0; jcolM[1]=1;jcolM[2]=0;jcolM[3]=1;jcolM[4]=2;jcolM[5]=3;
//   jcolM[6]=4;

//   M[0]=2.0;M[1]=0.1;
//   M[2]=0.1;M[3]=2.0;
//   M[4]=1.0;
//   M[5]=1.0;
//   M[6]=1.0;
  for(int i=0; i<n; i++) {
    krowM[i] = i;
    jcolM[i] = i;
    M[i] = 1.0;
  }
  krowM[n] = n;

  return 0;
}
extern "C"
int fQ1(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==1);
  int n=2;
  for(int i=0; i<n; i++) {
    krowM[i] = i;
    jcolM[i] = i;
    M[i] = 1.0;
  }
  krowM[n] = n;
  return 0;
}
extern "C"
int fQ2(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==2);
  int n=2;
  for(int i=0; i<n; i++) {
    krowM[i] = i;
    jcolM[i] = i;
    M[i] = 1.0;
  }
  krowM[n] = n;
  return 0;
}
extern "C" 
int fnnzQ0(void* user_data, int id, int* nnz)
{
  *nnz=5; return 0;
}
extern "C" 
int fnnzQ1(void* user_data, int id, int* nnz)
{
  *nnz=2; return 0;
}
extern "C" 
int fnnzQ2(void* user_data, int id, int* nnz)
{
  *nnz=2; return 0;
}
extern "C"
int fA0(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==0);
  int m=2, n=5, nnz=4;
  
  krowM[0]=0; krowM[1]=2; krowM[2]=4;
  jcolM[0]=0; jcolM[1]=4; jcolM[2]=1; jcolM[3]=3;
  M[0]=M[1]=M[2]=M[3]=1.0;
  return 0;
}
extern "C"
int fA1(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==1);
  int m=1, n=5, nnz=2;
  
  krowM[0]=0; krowM[1]=2;
  jcolM[0]=2; jcolM[1]=3;
  M[0]=M[1]=1.0;
  return 0;
}
extern "C"
int fA2(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==2);
  int m=1, n=5, nnz=3;
  
  krowM[0]=0; krowM[1]=3;
  jcolM[0]=1; jcolM[1]=2; jcolM[2]=3;
  M[0]=M[1]=M[2]=1.0;
  return 0;
}
extern "C" 
int fnnzA0(void* user_data, int id, int* nnz)
{
  *nnz=4; return 0;
}
extern "C" 
int fnnzA1(void* user_data, int id, int* nnz)
{
  *nnz=2; return 0;
}
extern "C" 
int fnnzA2(void* user_data, int id, int* nnz)
{
  *nnz=3; return 0;
}
extern "C"
int fB0(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==0);
  return 0;
}
extern "C"
int fB1(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==1);
  int m=1, n=2, nnz=1;
  
  krowM[0]=0; krowM[1]=1;
  jcolM[0]=1;
  M[0]=1.0;
  return 0;
}
extern "C"
int fB2(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==2);
  int m=1, n=2, nnz=1;
  
  krowM[0]=0; krowM[1]=1;
  jcolM[0]=0;
  M[0]=-1.0;
  return 0;
}
extern "C" 
int fnnzB0(void* user_data, int id, int* nnz)
{
  *nnz=0; return 0;
}

extern "C" 
int fnnzB1(void* user_data, int id, int* nnz)
{
  *nnz=1; return 0;
}
extern "C" 
int fnnzB2(void* user_data, int id, int* nnz)
{
  *nnz=1; return 0;
}


extern "C"
int fC0(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==0);
  int m=2, n=5, nnz=7;

  krowM[0]=0; krowM[1]=3; krowM[2]=7;
  jcolM[0]=0; jcolM[1]=3; jcolM[2]=4; jcolM[3]=0; jcolM[4]=1; jcolM[5]=3; jcolM[6]=4; 
  M[0] =-2.0; M[1] = M[2] = 1.0;      M[3] =-1.0; M[4] = 1.0; M[5] = 1.0; M[6] = 1.0;
  return 0;
}
extern "C"
int fC1(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==1);
  int m=1, n=5, nnz=2;

  krowM[0]=0; krowM[1]=2;
  jcolM[0]=2; jcolM[1]=3; 
  M[0] = M[1] = 1.0;
  return 0;
}
extern "C"
int fC2(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==2);
  int m=2, n=5, nnz=2;
  
  krowM[0]=0; krowM[1]=1; krowM[2]=2;
  jcolM[0]=0; jcolM[1]=4;
  M[0] = M[1] = 1.0;
  return 0;
}
extern "C" 
int fnnzC0(void* user_data, int id, int* nnz)
{
  *nnz=7; return 0;
}
extern "C" 
int fnnzC1(void* user_data, int id, int* nnz)
{
  *nnz=2; return 0;
}
extern "C" 
int fnnzC2(void* user_data, int id, int* nnz)
{
  *nnz=2; return 0;
}




extern "C" 
int fD0(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==0);
  return 0;
}
extern "C" 
int fD1(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==1); 
  int m=1, n=2, nnz=2;

  krowM[0]=0; krowM[1]=2;
  jcolM[0]=0; jcolM[1]=1;
  M[0] = M[1] = 1.0;
  return 0;
}
extern "C" 
int fD2(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  assert(id==2);
  int m=2, n=2, nnz=3;

  krowM[0]=0; krowM[1]=1; krowM[2]=3;
  jcolM[0]=0; jcolM[1]=0; jcolM[2]=1;
  M[0] = 7.0; M[1] = 1.0; M[2] = 17.0;
  return 0;
}
extern "C" 
int fnnzD0(void* user_data, int id, int* nnz) 
{
  *nnz=0; return 0;
}
extern "C" 
int fnnzD1(void* user_data, int id, int* nnz)
{
  *nnz=2; return 0;
}
extern "C"
int fnnzD2(void* user_data, int id, int* nnz)
{
  *nnz=3; return 0;
}

/***** vectors functions **************/

extern "C"
int fc0(void* user_data, int id, double* vec, int len)
{
  int n=5;
  assert(n==len);
  for(int i=0; i<n; i++) vec[i] = 1.0;
  return 0;
}
extern "C"
int fb0(void* user_data, int id, double* vec, int len)
{
  int my=2;
  assert(my=len);
  vec[0] = 1.0; vec[1] = 2.0;
  return 0;
}

extern "C"
int fclow0 (void* user_data, int id, double* vec, int len)
{
  int mz=2;
  assert(len==mz);
  vec[0] = 0.0; vec[1] = 1.0; 
  return 0;
}
extern "C"
int ficlow0(void* user_data, int id, double* vec, int len)
{
  int mz=2;
  assert(len==mz);
  vec[0] = 1; vec[1] = 1; // both active
  return 0;
}
extern "C"
int fcupp0 (void* user_data, int id, double* vec, int len)
{
  int mz=2;
  assert(len==mz);
  vec[0] = 2.0; vec[1] = 0.0;
  return 0;
}
extern "C"
int ficupp0(void* user_data, int id, double* vec, int len)
{
  int mz=2;
  assert(len==mz);
  vec[0] = 1; vec[1] = 0; //only the first constraint is active  
  return 0;
}

extern "C"
int fxlow0 (void* user_data, int id, double* vec, int len)
{
  int nx=5;
  assert(len==nx);
  
  for(int i=0; i<nx; i++) vec[i] = 0.0;
  return 0;
}
extern "C"
int fixlow0(void* user_data, int id, double* vec, int len)
{
  int nx=5;
  assert(len==nx);
  for(int i=0; i<nx; i++) vec[i] = 1;
  return 0;
}
extern "C"
int fxupp0 (void* user_data, int id, double* vec, int len)
{
  int nx=5;
  assert(len==nx);

  //nothing needed here, there is no upper limit
 
  for(int i=0; i<nx; i++) vec[i] = 0.0;
  return 0;
}
extern "C"
int fixupp0(void* user_data, int id, double* vec, int len)
{
  int nx=5;
  assert(len==nx);
  for(int i=0; i<nx; i++) vec[i] = 0;
  return 0;
}




extern "C"
int fc1(void* user_data, int id, double* vec, int len)
{
  int nx=2;
  assert(len==nx);
  vec[0]=vec[1]=1.0;
  return 0;
}
extern "C"
int fb1(void* user_data, int id, double* vec, int len)
{
  int my=1;
  assert(len==my);
  vec[0] = 4;
  return 0;
}

extern "C"
int fclow1 (void* user_data, int id, double* vec, int len)
{
  int mz=1;
  assert(len==mz);
  vec[0]=5.0;
  return 0;
}
extern "C"
int ficlow1(void* user_data, int id, double* vec, int len)
{
  int mz=1;
  assert(len==mz);
  vec[0]=1;  
  return 0;
}
extern "C"
int fcupp1 (void* user_data, int id, double* vec, int len)
{
  int mz=1;
  assert(len==mz);
  vec[0]=10.0;  
  return 0;
}
extern "C"
int ficupp1(void* user_data, int id, double* vec, int len)
{
  int mz=1;
  assert(len==mz);
  vec[0]=1;  
  return 0;
}

extern "C"
int fxlow1 (void* user_data, int id, double* vec, int len)
{
  int nx=2;
  assert(nx==len);
  vec[0]=vec[1]=0.0;
  return 0;
}
extern "C"
int fixlow1(void* user_data, int id, double* vec, int len)
{
  int nx=2;
  assert(nx==len);
  vec[0]=vec[1]=1;
  return 0;
}
extern "C"
int fxupp1 (void* user_data, int id, double* vec, int len)
{
  //not needed, there are no upper limits
  int nx=2;
  assert(nx==len);
  vec[0]=vec[1]=0.0;   
  return 0;
}
extern "C"
int fixupp1(void* user_data, int id, double* vec, int len)
{
  int nx=2;
  assert(nx==len);
  vec[0]=vec[1]=0; //not active
  return 0;
}





extern "C"
int fc2(void* user_data, int id, double* vec, int len)
{
  int n=2;
  assert(len==n);
  vec[0]=vec[1]=1.0;
  return 0;
}
extern "C"
int fb2(void* user_data, int id, double* vec, int len)
{
  int my=1;
  assert(len==my);
  vec[0] = -2;
  return 0;
}

extern "C"
int fclow2 (void* user_data, int id, double* vec, int len)
{
  int mz=2;
  assert(mz==len);
  vec[0]=0.0; vec[1]=1.0; //both active
  return 0;
}
extern "C"
int ficlow2(void* user_data, int id, double* vec, int len)
{
  int mz=2;
  assert(mz==len);
  vec[0]=vec[1]=1.0; //both active
  return 0;
}
extern "C"
int fcupp2 (void* user_data, int id, double* vec, int len)
{
  int mz=2;
  assert(mz==len);
  vec[0]=35.0; vec[1]=7.0; //both active
  return 0;
}
extern "C"
int ficupp2(void* user_data, int id, double* vec, int len)
{
  int mz=2;
  assert(mz==len);
  vec[0]=vec[1]=1; //both active
  return 0;
}

extern "C"
int fxlow2 (void* user_data, int id, double* vec, int len)
{
  //does not have to do anything since the ixlow is all 0;
  //anyhow, let's set it to 0
  vec[0]=vec[1]=0.0;
  return 0;
}
extern "C"
int fixlow2(void* user_data, int id, double* vec, int len)
{
  int n=2;
  assert(n==len);
  vec[0]=vec[1]=0; //not active
  return 0;
  return 0;
}
extern "C"
int fxupp2 (void* user_data, int id, double* vec, int len)
{
  int n=2;
  assert(n==len);
  vec[0]=vec[1]=6.0;
  return 0;
  return 0;
}
extern "C"
int fixupp2(void* user_data, int id, double* vec, int len)
{
  int n=2;
  assert(n==len);
  vec[0]=vec[1]=1; //both active
  return 0;
}

/*********************************************
 *          callback for dynamical problem   *
 *********************************************/

extern "C"
int fnnzQ(void* user_data, int id, int* nnz)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  *nnz = mydata->n;
}

extern "C"
int fnnzA(void* user_data, int id, int* nnz)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  *nnz = 5*(mydata->n/5);
}

extern "C"
int fnnzB(void* user_data, int id, int* nnz)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  *nnz = 5*(mydata->n/5);
}

extern "C"
int fnnzC(void* user_data, int id, int* nnz)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  *nnz = 2*(mydata->n/2);
}

extern "C"
int fnnzD(void* user_data, int id, int* nnz)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  *nnz = 2*(mydata->n/2);
}

extern "C"
int fQ(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  MyNodeData* mydata = (MyNodeData*)user_data;

  int n=mydata->n;
  for(int i=0; i<n; i++) {
    krowM[i] = i;
    jcolM[i] = i;
    M[i] = 1.0;
  }
  krowM[n] = n;

  return 0;
}

#define ERE 5 //= nonzero elements per row in eq. constraints
#define IRE 2 //= nonzero elements per row in ineq. constraints

extern "C"
int fA(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  int n=mydata->n;

  assert(n == ERE*(n/ERE));

  for(int i=0; i<n/ERE; i++) {
    krowM[i]=ERE*i;

    for(int j=0; j<ERE; j++) {      
      jcolM[ERE*i+j] = ERE*i+j;
      M    [ERE*i+j] = 1.0;
    }
  }
  krowM[n/ERE]= n;

  return 0;
}

extern "C"
int fB(void* user_data, int id_, int* krowM, int* jcolM, double* M)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  int id = mydata->id;

  if(id==0) return 0;

  int n  = mydata->n;

  for(int i=0; i<n/ERE; i++) {
    krowM[i]=ERE*i;

    for(int j=0; j<ERE; j++) {      
      jcolM[ERE*i+j] = ERE*i+j;
      M    [ERE*i+j] = 1.0;
    }
  }
  krowM[n/ERE]= n;

  return 0;
}

extern "C"
int fC(void* user_data, int id_, int* krowM, int* jcolM, double* M)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  int n=mydata->n;

  assert(n == IRE*(n/IRE));

  for(int i=0; i<n/IRE; i++) {
    krowM[i]=IRE*i;

    for(int j=0; j<IRE; j++) {      
      jcolM[IRE*i+j] = IRE*i+j;
      M    [IRE*i+j] = 1.0;
    }
  }
  krowM[n/IRE]= n;

  return 0;
}

extern "C"
int fD(void* user_data, int id_, int* krowM, int* jcolM, double* M)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  int id = mydata->id;

  if(id==0) return 0;

  int n  = mydata->n;

  for(int i=0; i<n/IRE; i++) {
    krowM[i]=IRE*i;

    for(int j=0; j<IRE; j++) {      
      jcolM[IRE*i+j] = IRE*i+j;
      M    [IRE*i+j] = 1.0;
    }
  }
  krowM[n/IRE]= n;

  return 0;
}

extern "C"
int fc(void* user_data, int id, double* vec, int len)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  int n  = mydata->n;

  assert(n==len);
  for(int i=0; i<len; i++) vec[i]=0.0;//n/100.0;

  //vec[len-1]=-1.5;

  return 0;
}

extern "C"
int fb(void* user_data, int id, double* vec, int len)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  int n  = mydata->n;

  double rhs=ERE;
  if(mydata->id>0) rhs = 2*rhs;

  //if(mydata->id==2) rhs = rhs/1.5;

  assert(len==n/ERE);
  for(int i=0; i<len; i++) vec[i] = rhs; 

  return 0;
}

extern "C"
int fclow (void* user_data, int id, double* vec, int len)
{
  MyNodeData* mydata = (MyNodeData*)user_data;

  int n  = mydata->n;
  //none active 
  assert(len==n/IRE);
  for(int i=0; i<len; i++) vec[i] = 0.0;

  return 0;
}

extern "C"
int ficlow(void* user_data, int id, double* vec, int len)
{
  MyNodeData* mydata = (MyNodeData*)user_data;

  int n  = mydata->n;
  assert(len==n/IRE);
  for(int i=0; i<len; i++) vec[i] = 0.0;
  return 0;
}

extern "C"
int fcupp (void* user_data, int id, double* vec, int len)
{
  MyNodeData* mydata = (MyNodeData*)user_data;

  int n  = mydata->n;

  double rhs = IRE;
  if(mydata->id>0) rhs = rhs*2.0; 

  //all active
  assert(len==n/IRE);
  for(int i=0; i<len; i++) vec[i] = rhs;

  return 0;
}

extern "C"
int ficupp(void* user_data, int id, double* vec, int len)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
 
  int n  = mydata->n; 
  //all active
  assert(len==n/IRE);
  for(int i=0; i<len; i++) vec[i] = 1;
  return 0;
}

extern "C"
int fxlow (void* user_data, int id, double* vec, int len)
{
  MyNodeData* mydata = (MyNodeData*)user_data;

  int n  = mydata->n;
  //all active
  assert(len==n);
  for(int i=0; i<len; i++) vec[i] = 0.0;
  return 0;
}

extern "C"
int fixlow(void* user_data, int id, double* vec, int len)
{
  MyNodeData* mydata = (MyNodeData*)user_data;

  int n  = mydata->n;
  assert(len==n);
  for(int i=0; i<len; i++) vec[i] = 1;
  return 0;
}

extern "C"
int fxupp (void* user_data, int id, double* vec, int len)
{
  MyNodeData* mydata = (MyNodeData*)user_data;

  int n  = mydata->n;
  //none active
  assert(len==n);
  for(int i=0; i<len; i++) vec[i] = 0.0;
  return 0;
}

extern "C"
int fixupp(void* user_data, int id, double* vec, int len)
{
  MyNodeData* mydata = (MyNodeData*)user_data;

  int n  = mydata->n;
  //none active
  assert(len==n);
  for(int i=0; i<len; i++) vec[i] = 0;
  return 0;
}
