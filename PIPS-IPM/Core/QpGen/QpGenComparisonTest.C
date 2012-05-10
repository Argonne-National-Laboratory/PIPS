#include "StochTree.h"
#include "QpGenStochDriver.h"
#include "QpGenDriver.h"
#include "QpGenSparseMa57.h"
#include "MehrotraSolver.h"

int solveQPSparse(int levels, int printx);
int solveQPStoch(int levels, int printx);


extern "C" {
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

int main( int argc, char *argv[] )
{
  int goodUsage = 1; // Assume good. Why not be optimistic
  int iargv    = 1;
  int solveWith = -1, levels = -1, printx=0;
  while( iargv < argc ) {
    if( argv[iargv][0] != '-' ) {break;} // done with the options
    int iopt = 1;
    if( argv[iargv][1] == '-' ) iopt = 2; // it is a --arg
    
    if( 0 == strcmp( "solve", &argv[iargv][iopt] ) ) {

      iargv++;
      if(iargv>=argc) break;

      char* endptr;
      solveWith = strtol( argv[iargv], &endptr, 10 );
      if( '\0' != *endptr ) { 
	goodUsage = 0;assert(0);
	break;
      }
    } else if( 0 == strcmp("levels", &argv[iargv][iopt]) ) {

      iargv++;
      if(iargv>=argc) break;

      char* endptr;
      levels = strtol( argv[iargv], &endptr, 10 );
      if( '\0' != *endptr ) { 
	goodUsage = 0;
	break;
      }
    } else if( 0 == strcmp( "printx", &argv[iargv][iopt] ) ) {
      printx=1;
    } else {
      goodUsage=0;
      break;
    }
    iargv++;
  } // end while

  if( argc>1 && iargv == argc && goodUsage ) {

    if(solveWith!=0 && solveWith!=1 && solveWith!=2) {
      cerr << "Solver(s) is not recognized/specified.\n";
      return 1;
    }

    if(levels!=2 && levels!=3) {
      cerr << "Number of levels is not valid/specified.\n";
      return 1;
    }
    
  } else {
    
    cerr << "Usage: " << argv[0] << "  --solve num --levels num [--printx]\n\n";
    cerr << "\t--solve num  - solver used: 0 STOCH, 1 SPARSE, 2 BOTH\n";
    cerr << "\t--levels num - number of levels the problem should have, 2 or 3 supported.\n";
    cerr << "\t--printx     - print the 'x' portion of the solution.\n";
    cerr <<"\n";
    return 1;
  }

  if( solveWith==0 ) {
    solveQPStoch(levels,printx);
  } else if( solveWith==1 ) {
    solveQPSparse(levels,printx);
  } else {
    solveQPStoch(levels,printx);
    solveQPSparse(levels,printx);
  }

  return 0;
}

int solveQPStoch(int levels, int printx)
{
  void* usrdata=NULL;
  StochTree::TreeNode node( usrdata, 0, 
			    5, 2, 2,
			    fQ0, fnnzQ0, fc0,
			    fA0, fnnzA0,
			    fB0,  fnnzB0,
			    fb0,
			    fC0, fnnzC0,
			    fD0, fnnzD0,
			    fclow0, ficlow0, fcupp0, ficupp0,
			    fxlow0, fixlow0, fxupp0, fixupp0 );
  
  StochTree::TreeNode node1( usrdata, 1, 
			     2, 1, 1,
			     fQ1, fnnzQ1, fc1,
			     fA1, fnnzA1,
			     fB1, fnnzB1,
			     fb1,
			     fC1, fnnzC1,
			     fD1, fnnzD1,
			     fclow1, ficlow1, fcupp1, ficupp1,
			     fxlow1, fixlow1, fxupp1, fixupp1 );
  
  StochTree::TreeNode node2( usrdata, 2, 
			     2, 1, 2,
			     fQ2, fnnzQ2, fc2,
			     fA2, fnnzA2,
			     fB2, fnnzB2,
			     fb2,
			     fC2, fnnzC2,
			     fD2, fnnzD2,
			     fclow2, ficlow2, fcupp2, ficupp2,
			     fxlow2, fixlow2, fxupp2, fixupp2 );
  
  StochTree root(node);
  root.AddChild(node1);
  root.AddChild(node2);

  StochRunParams* params = defaultStochRunParams();
  params->printx = printx;

  MehrotraSolver* method=NULL;
  qpgenstoch_solve(&root, method, params);

  delete params;
  return 0;
}


int solveQPSparse(int levels, int printx)
{
  int n1, m1, m2, nnzQ, nnzA, nnzC;
  int nx=9; int my=4; int mz=5;
    
  nnzQ = 9; nnzA = 4+ 2+1 + 3+1; nnzC=7+  2+2 + 2+3; 

  //Q
  int* krowQ=new int[nx+1]; int*jcolQ=new int[nnzQ]; double* dQ = new double[nnzQ];
  for(int i=0; i<=nx; i++)  krowQ[i]=i;
  for(int i=0; i<nnzQ; i++) {jcolQ[i]=i;dQ[i]=1.0;}

  //A
  int* krowA=new int[my+1]; int*jcolA=new int[nnzA]; double* dA = new double[nnzA];
  krowA[0]=0; jcolA[0]=0; jcolA[1]=4;                       dA[0]=1.0;dA[1]=1.0;
  krowA[1]=2; jcolA[2]=1; jcolA[3]=3;                       dA[2]=1.0;dA[3]=1.0;
  krowA[2]=4; jcolA[4]=2; jcolA[5]=3;jcolA[6]=6;            dA[4]=1.0;dA[5]=1.0;dA[6]=1.0;
  krowA[3]=7; jcolA[7]=1; jcolA[8]=2;jcolA[9]=3;jcolA[10]=7;dA[7]=1.0;dA[8]=1.0;dA[9]=1.0;dA[10]=-1.0;
  krowA[4]=11;
  
  //C
  int* krowC=new int[mz+1]; int*jcolC=new int[nnzC]; double* dC = new double[nnzC];
  krowC[0]=0; jcolC[0]=0;jcolC[1]=3;jcolC[2]=4;
  dC[0]=-2.0;dC[1]=1.0; dC[2]=1.0;
  
  krowC[1]=3;jcolC[3]=0;jcolC[4]=1;jcolC[5]=3;jcolC[6]=4;
  dC[3]=-1.0;dC[4]=1.0;dC[5]=1.0;dC[6]=1.0;
  
  krowC[2]=7;jcolC[7]=2;jcolC[8]=3;jcolC[9]=5;jcolC[10]=6;
  dC[7]=1.0;dC[8]=1.0;dC[9]=1.0;dC[10]=1.0;
  
  krowC[3]=11; jcolC[11]=0;jcolC[12]=7; dC[11]=1.0;dC[12]=7.0;
  
  krowC[4]=13; jcolC[13]=4;jcolC[14]=7;jcolC[15]=8; dC[13]=1.0;dC[14]=1.0;dC[15]=17.0;
  krowC[5]=16;

  //xlow
  char* ixlow=new char[nx];double* xlow=new double[nx];
  for(int i=0; i<nx; i++) {ixlow[i]=1;xlow[i]=0.0;}
  ixlow[7]=0;ixlow[8]=0;
  
  //xupp
  char* ixupp=new char[nx];double* xupp=new double[nx];
  for(int i=0; i<nx; i++) {ixupp[i]=0;xupp[i]=0.0;}
  ixupp[7]=1; ixupp[8]=1;
  xupp[7]=6.0;xupp[8]=6.0;
  
  //clow
  char* iclow=new char[mz];double* clow=new double[mz];
  iclow[0]=1; iclow[1]=1; iclow[2]=1; iclow[3]=1; iclow[4]=1;
  clow[0]=0.0;clow[1]=1.0;clow[2]=5.0;clow[3]=0.0;clow[4]=1.0;
  
  //cupp
  char* icupp=new char[mz];double* cupp=new double[mz];
  icupp[0]=1; icupp[1]=0; icupp[2]=1; icupp[3]=1; icupp[4]=1; 
  cupp[0]=2.0;cupp[1]=0.0;cupp[2]=10;cupp[3]=35.0;cupp[4]=7.0;
  
  
  //b
  double* b = new double[my];
  b[0]=1.0;b[1]=2.0;b[2]=4.0;b[3]=-2.0;
  
  //c 
  double* c = new double[nx];
  for(int i=0;i<nx; i++)c[i]=1.0;
  
  QpGenSparseMa57 * qp = new QpGenSparseMa57( nx, my, mz, nnzQ, nnzA, nnzC );
  QpGenData * prob = (QpGenData * ) qp->makeData(c, 
						 krowQ, jcolQ, dQ,
						 xlow, ixlow,
						 xupp, ixupp,
						 krowA, jcolA, dA,
						 b,
						 krowC, jcolC, dC,
						 clow, iclow,
						 cupp, icupp);
  
  
  
  QpGenVars     * vars  = (QpGenVars * ) qp->makeVariables( prob );
  Residuals     * resid = qp->makeResiduals( prob );
  MehrotraSolver * s            = new MehrotraSolver( qp, prob );
  s->monitorSelf();
  
  int result = s->solve(prob,vars, resid);
  
  double objective = prob->objectiveValue(vars);
  
  cout << " " << nx << " variables, " 
       << my  << " equality constraints, " 
       << mz  << " inequality constraints.\n";
  
  cout << " Iterates: " << s->iter
       <<",    Optimal Solution:  " << objective;
  cout << ".\n";

  if(printx) {
    cout << "The x-solution is:\n";
    vars->x->writeToStream(cout);
  }
  delete s;
  delete vars;  
  delete resid;
  delete prob;
  delete qp;
  
  return result;
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

