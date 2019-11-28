#include "StochInputTree.h"
#include "QpGenDriver.h"
#include "QpGenStochDriver.h"
#include "MehrotraSolver.h"
#include "MehrotraStochSolver.h"
#include "GondzioSolver.h"
#include "QpGenSparseMa57.h"
#include "QpGenStochAugExt.h"
#include "pipsport.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
using namespace std;

#include "mpi.h"

/** User defined node data */
struct MyNodeData 
{
public:
  MyNodeData(int nt_, int nx_) : nx(nx_), nt(nt_) {}
  int nx, nt;
};


extern "C" {
/**** callbacks for BUILDING model */
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
int fixlowc(void* user_data, int id, char*   vec, int len);
int fxupp  (void* user_data, int id, double* vec, int len);
int fixupp (void* user_data, int id, double* vec, int len);
int fixuppc(void* user_data, int id, char*   vec, int len);
};

/** Declarations of utilities
 */
string g_strDataDir = "../../apps/building/";
void loadWeatherData(const char* szFile, int scenario, int tSize, double* data);
void loadElectricityPrices(const char* szFile, int tSize, double* data);
void stochSolve(int argc, char *argv[], int ns, int nt, int nx, int printx);
void qpgenSolve(int ns, int nt, int nx, int printx);
/** 
 * Usage: building.exe -solve 0 -ns 100 -nt 100 -nx 100 
 */

int main( int argc, char *argv[] )
{ 
  int ns, nt, nx;
  ns = 100; nt = 12; nx = 10;  //! 120 100

  int goodUsage = 1; // Assume good. Why not be optimistic
  int iargv    = 1;
  int solveWith = -1, printx=0;
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
    } else if( 0 == strcmp("ns", &argv[iargv][iopt]) ) {

      iargv++;
      if(iargv>=argc) break;

      char* endptr;
      ns = strtol( argv[iargv], &endptr, 10 );
      if( '\0' != *endptr ) { 
	goodUsage = 0;
	break;
      }
    } else if( 0 == strcmp("nt", &argv[iargv][iopt]) ) {

      iargv++;
      if(iargv>=argc) break;

      char* endptr;
      nt = strtol( argv[iargv], &endptr, 10 );
      if( '\0' != *endptr ) { 
	goodUsage = 0;
	break;
      }
    } else if( 0 == strcmp("nx", &argv[iargv][iopt]) ) {

      iargv++;
      if(iargv>=argc) break;

      char* endptr;
      nx = strtol( argv[iargv], &endptr, 10 );
      if( '\0' != *endptr ) { 
	goodUsage = 0;
	break;
      }
    }else if( 0 == strcmp( "printx", &argv[iargv][iopt] ) ) {
      printx=1;
    } else if( 0 == strcmp("datadir", &argv[iargv][iopt]) ) {

      iargv++;
      if(iargv>=argc) break;

      char* endptr;
      g_strDataDir = argv[iargv];//strtol( argv[iargv], &endptr, 10 );
      if('/' != g_strDataDir[g_strDataDir.length()-1]) g_strDataDir += "/";
	
    } else {
      goodUsage=0;
      break;
    }
    iargv++;
  } // end while

  if( argc>1 && iargv == argc && goodUsage ) {

    if(solveWith!=0 && solveWith!=1) {
      cerr << "Solver(s) is not recognized/specified.\n";
      return 1;
    }
    //!if(ns!=-100) {
    //  cerr << "Currently the only supported value for 'ns' is 100.\n";
    //   return 1;
    // }
  } else {
    
    cerr << "Usage: " << argv[0] << "  --solve num --nx num --nt num --ns num [--printx]\n\n";
    cerr << "\t--solve num    - solver used: 0 STOCH, 1 SPARSE\n";
    cerr << "\t--ns num       - number of scenarios to use\n";
    cerr << "\t--nt num       - time discretization\n";
    cerr << "\t--nx num       - spatial 2D discretization of the wall\n";
    cerr << "\t--datadir path - path the directory containing data of building problem\n";
    cerr << "\t--printx       - print the 'x' portion of the solution.\n";
    cerr <<"\n";
    return 1;
  }




  if(solveWith==0) {
    //stochastic solve
    stochSolve(argc, argv, ns, nt, nx, printx);
  } else {
    //QpGen solve
    qpgenSolve(ns, nt, nx, printx);
  }


  return 0;
}



void stochSolve(int argc, char *argv[], int ns, int nt, int nx, int printx)
{
  //initialize MPI env
  MPI_Init(&argc, &argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int size; MPI_Comm_size(MPI_COMM_WORLD, &size);

  //generate input stochastic(scenario) tree
  int globalID=0; 


  if(rank==0) {
    printf("1st Stage: %d variables, 0 eq, 0 ineq.\n", 3*(nt-1));
    printf("2nd Stage: %d variables, %d eq, 0 ineq.\n", nt+nt*nx, nt*(nx+1));
  }

  StochInputTree::StochInputNode data(
		       new MyNodeData(nt,nx), globalID, 
		       3*(nt-1), 0, 0, 
		       fQ, fnnzQ, fc,
		       fA, fnnzA,
		       fB,  fnnzB,
		       fb,
		       fC, fnnzC,
		       fD, fnnzD,
		       fclow, ficlow, fcupp, ficupp,
		       fxlow, fixlow, fxupp, fixupp );
  globalID++;
  StochInputTree* root = new StochInputTree(data);
  
  for(int i=0; i<ns; i++) {
    
    //int my = (nt-1) + (nt-1)*(nx-2) + (nt-1) + (nt-1) + 1 + nx;
    int my = nt*(nx+1);
    int n  = nt+nt*nx;
    int mz = 0;
    StochInputTree::StochInputNode data(
                         new MyNodeData(nt,nx), globalID, 
			 n, my, mz, 
			 fQ, fnnzQ, fc,
			 fA, fnnzA,
			 fB,  fnnzB,
			 fb,
			 fC, fnnzC,
			 fD, fnnzD,
			 fclow, ficlow, fcupp, ficupp,
			 fxlow, fixlow, fxupp, fixupp );
    globalID++;
    
    root->AddChild(new StochInputTree(data));
  }


  //call the stochastic IPM solver
  StochRunParams* params = defaultStochRunParams();
  params->printx = printx;

  MehrotraStochSolver* method=nullptr;
  QpGenStochAugExt* formulation=nullptr;
  qpgenstoch_solve( root, params, method, formulation);
  
  delete params;
  delete root;

  MPI_Finalize();

}

void qpgenSolve(int ns, int nt, int nx, int printx)
{
  int NX, MY, MZ=0;
  int nnzA=0, nnzQ=0, nnzC=0;

  NX = 3*(nt-1) + ns*nt*(nx+1);
  MY = ns*nt*(nx+1);

  MyNodeData* nd = new MyNodeData(nt,nx);
  int tmp;
  for(int i=1; i<=ns; i++) {
    fnnzA(nd, i, &tmp); nnzA+= tmp;
    fnnzB(nd, i, &tmp); nnzA+= tmp;
  }

  //Q
  int* krowQ=new int[NX+1]; int*jcolQ=new int[nnzQ]; double* dQ = new double[nnzQ];
  for(int i=0; i<=NX; i++) krowQ[i] = 0;

  //A and b
  int* krowA=new int[MY+1]; int*jcolA=new int[nnzA]; double* dA = new double[nnzA];
  double* b = new double[MY];
  // populate A
  // part of A corresponding to root is empty-> do nothing
  int shift=0; int jcol=0;
  for(int s=1; s<=ns; s++) {
    int NNZA; fnnzA(nd, s, &NNZA);
    int NNZB; fnnzB(nd, s, &NNZB);
    
    int* krowAs=new int[nt*(nx+1)+1]; int*jcolAs=new int[NNZA]; double* dAs = new double[NNZA];
    int* krowBs=new int[nt*(nx+1)+1]; int*jcolBs=new int[NNZB]; double* dBs = new double[NNZB];

    fA(nd, s, krowAs, jcolAs, dAs);
    fB(nd, s, krowBs, jcolBs, dBs);
    
    for(int i=0; i<nt*(nx+1); i++) {
      krowA[i+shift] = jcol;

      for(int k=krowAs[i]; k<krowAs[i+1]; k++) {
	jcolA[jcol] = jcolAs[k];
	dA   [jcol] = dAs[k];
	jcol++;
      }
      
      for(int k=krowBs[i]; k<krowBs[i+1]; k++) {
	jcolA[jcol] = jcolBs[k] + 3*(nt-1) + (s-1)*nt*(nx+1);
	dA   [jcol] = dBs   [k];
	jcol++;
      }

      krowA[i+shift+1] = jcol;
    }
    delete[] krowAs; delete[] jcolAs; delete[] dAs;
    delete[] krowBs; delete[] jcolBs; delete[] dBs;

    //b
    fb(nd, s, &b[shift], nt*(nx+1));

    shift += nt*(nx+1);
    assert(shift==s*nt*(nx+1));
  }
  assert(shift==ns*nt*(nx+1));
  assert(jcol==nnzA);
  //~ done with A and b

   //C
  int* krowC=new int[MZ+1]; int*jcolC=new int[nnzC]; double* dC = new double[nnzC];
  for(int i=0; i<=MZ; i++) krowC[i]=0;

  //xlow, ixlow, xupp, ixupp
  char* ixlow=new char[NX];double* xlow=new double[NX];
  char* ixupp=new char[NX];double* xupp=new double[NX];
  shift=0;
  fxlow  (nd, 0, &xlow[shift] , 3*(nt-1));
  fixlowc(nd, 0, &ixlow[shift], 3*(nt-1));
  fxupp  (nd, 0, &xupp[shift] , 3*(nt-1));
  fixuppc(nd, 0, &ixupp[shift], 3*(nt-1));
  shift += 3*(nt-1);
  for(int s=1; s<=ns; s++) {
    fxlow  (nd, s, &xlow[shift] , nt*(nx+1));
    fixlowc(nd, s, &ixlow[shift], nt*(nx+1));
    fxupp  (nd, s, &xupp[shift] , nt*(nx+1));
    fixuppc(nd, s, &ixupp[shift], nt*(nx+1));
    shift += nt*(nx+1);
  } 

  //clow, iclow, cupp, icupp
  char* iclow=new char[MZ];double* clow=new double[MZ];
  char* icupp=new char[MZ];double* cupp=new double[MZ];

  //c 
  double* c = new double[NX];
  shift=0;
  fc(nd, 0, &c[shift], 3*(nt-1));
  shift += 3*(nt-1);
  for(int s=1; s<=ns; s++) {
    fc(nd, s, &c[shift], nt*(nx+1));
    shift += nt*(nx+1);
  }
  delete nd;


  ////////////////////////////////////////////////////
  // Solve the problem
  ////////////////////////////////////////////////////
  QpGenSparseMa57 * qp = new QpGenSparseMa57( NX, MY, MZ, nnzQ, nnzA, nnzC );
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
  //GondzioSolver * s            = new GondzioSolver( qp, prob );
  MehrotraSolver * s            = new MehrotraSolver( qp, prob );
  s->monitorSelf();
  
  rusage before_solve;
  getrusage( RUSAGE_SELF, &before_solve );

  int result = s->solve(prob,vars, resid);

  rusage  after_solve;
  getrusage( RUSAGE_SELF, &after_solve );  

  double solve_time =
    (after_solve.ru_utime.tv_sec - before_solve.ru_utime.tv_sec)
    + (after_solve.ru_utime.tv_usec - before_solve.ru_utime.tv_usec)
    / 1000000.0;
  
  double objective = prob->objectiveValue(vars);

  if(printx) {
    cout << "The x-solution is:\n";
    vars->x->writeToStream(cout);
  }  

  cout << " " << NX << " variables, " 
       << MY  << " equality constraints, " 
       << MZ  << " inequality constraints.\n";
  
  cout << " Iterates: " << s->iter
       <<",    Optimal Solution:  " << objective;
  cout << ".\n";

  cout  << "QP solved in " << solve_time << " seconds.\n";

  delete s;
  delete vars;  
  delete resid;
  delete prob;
  delete qp;
}

//model constants
static const double alphap  =       4.0;        //internal convective heat transfer coefficient (kcal/m2-hr-C)
static const double alphapp =      10.0;        //external convective heat transfer coefficient (kcal/m2-hr-C)
static const double       L =     0.3;          //wall thickness (m)
static const double      K  =     0.1;          //wall conductivity (kcal/m-hr-C)
static const double    beta =   0.001;          //wall thermal diffusivity -concrete (m2/hr) (a in AMPL model)
static const double       S =    3500.0;        //wall surface area (m2)
static const double       V =   10000.0;        //building volume (m3)
static const double      CI =   34800.0;        //internal air heat capacity (kcal/C)
static const double     CEL = 0.12*4.18/3600;   // on-peak price (0.12 $/kW) -> $/(kcal/hr)
static const double    CGAS = 0.10*4.18/3600;   // gas price (0.10 $/kW) -> $/(kcal/hr) 


extern "C"
int fnnzQ(void* user_data, int id, int* nnz)
{ *nnz = 0; }

extern "C"
int fnnzA(void* user_data, int id, int* nnz)
{ 
  if(id==0) {
    //no equalities
    *nnz = 0;
  } else {
    
    MyNodeData* mydata = (MyNodeData*)user_data;
    int nx = mydata->nx; int nt = mydata->nt;
    *nnz = 3*(nt-1); 
  }
}

extern "C"
int fnnzB(void* user_data, int id, int* nnz)
{ 
  if(id==0) {
    //no equalities
    *nnz = 0;
  } else {
    MyNodeData* mydata = (MyNodeData*)user_data;
    int nx = mydata->nx; int nt = mydata->nt;
    *nnz = (nt-1)*4*nx + nx+1;
  }
}

extern "C"
int fnnzC(void* user_data, int id, int* nnz)
{ *nnz = 0; }

extern "C"
int fnnzD(void* user_data, int id, int* nnz)
{ *nnz = 0; }

extern "C"
int fQ(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  int nx=mydata->nx; int nt = mydata->nt;
 
  if(id==0) {
    for(int i=0; i<=3*(nt-1);  i++) krowM[i]=0;
  } else {
    for(int i=0; i<=nt*(nx+1); i++) krowM[i]=0;
  }

  return 0;
}

extern "C"
int fc(void* user_data, int id, double* vec, int len)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  int nx=mydata->nx; int nt = mydata->nt;

  if(id==0) {
    assert(len==3*(nt-1));
    
    for(int i=0; i<nt-1; i++) vec[i] = CGAS;
    
    loadElectricityPrices("temp_tx_2008.dat", nt-1, &vec[nt-1]);
    
    for(int i=0; i<nt-1; i++) {
      vec[nt-1+i] *= CEL;
      vec[2*(nt-1)+i] = vec[nt-1+i];
    }

  } else {
    assert(len==nt*(nx+1));
    for(int i=0; i<len; i++) vec[i]=0.0;
  }  
}

extern "C"
int fA(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  int nx=mydata->nx; int nt = mydata->nt;
  if(id==0) {
    //no equalities
  } else {

    // nt(nx+1)-by-3(nt-1) matrix
    krowM[0]=0;
    for(int row=0; row<nt-1; row++) {
      jcolM[3*row]   = row;        M[3*row]   = -1.0;
      jcolM[3*row+1] = row+nt-1;   M[3*row+1] =  1.0; 
      jcolM[3*row+2] = row+2*nt-2; M[3*row+2] = -1.0;

      krowM[row+1]=3*row+3; 
    }
    for(int row=nt-1; row<nt*(nx+1); row++) 
      krowM[row+1] = krowM[row];
  }
  return 0;
}

extern "C"
int fB(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  int nx=mydata->nx; int nt = mydata->nt;
  
  if(id==0) {
    //no equalities
    for(int i=0; i<(nt-1)*4*nx + nx+1; i++) krowM[i]=0;

  } else {
    
    // nt(nx+1)-by-nt(nx+1) matrix

    ///////////////////////////////////////////
    //first nt-1 rows      -- 1 --
    ///////////////////////////////////////////
    double dt=1.0;
    for(int i=0; i<nt-1; i++) {
      int colind=3*i;
      krowM[i]=colind;

      jcolM[colind]=i; jcolM[colind+1]=i+1; jcolM[colind+2]=nt+i+1;
      M[colind]=-CI/dt; 
      M[colind+1]=CI/dt+S*alphap;
      M[colind+2]=-S*alphap;
    }
    int rowShift=nt-1; int colShift=3*(nt-1);
    ///////////////////////////////////////
    // rows  nt to (nt-1)(nx-2)  -- 2 --
    ///////////////////////////////////////
    double dx=L/nx; double dx2=dx*dx;
    for(int j=0; j<nx-2; j++) {
      for(int k=0; k<nt-1; k++) {
	int colind = colShift + 4*(j*(nt-1)+k);
	krowM[rowShift+j*(nt-1)+k] = colind;

	jcolM[colind]   = j*nt + k + nt+1;       M[colind  ]=  beta/dx2;
	jcolM[colind+1] = j*nt + k + nt+nt;      M[colind+1]=  1/dt;
	jcolM[colind+2] = j*nt + k + nt+nt+1;    M[colind+2]= -2*beta/dx2-1/dt;
	jcolM[colind+3] = j*nt + k + nt+nt+nt+1; M[colind+3]=  beta/dx2;
      }
    }
    rowShift+=(nt-1)*(nx-2); colShift+= 4*(nt-1)*(nx-2);
    ///////////////////////////////////////////////////////////////////
    // rows nt+(nt-1)(nx-2)  to   nt-1 + (nt-1)(nx-2) + nt-1   -- 3 --
    ///////////////////////////////////////////////////////////////////
    for(int k=0; k<nt-1; k++) {
      int colind = colShift + 3*k;
      krowM[rowShift+k] = colind;

      jcolM[colind]   = 1+k;       M[colind  ] =  alphap;
      jcolM[colind+1] = nt+1+k;    M[colind+1] = -alphap-K/dx;
      jcolM[colind+2] = nt+1+nt+k; M[colind+2] =  K/dx;
    }
    rowShift+=nt-1; colShift+=3*(nt-1);
    ///////////////////////////////////////////////////////////////////
    // rows nt-1 + (nt-1)(nx-2) + nt    to 
    //    nt-1 + (nt-1)(nx-2) + nt + nt -1   --- 4 --- 
    ///////////////////////////////////////////////////////////////////
    for(int k=0; k<nt-1; k++) {
      int colind = colShift + 2*k;
      krowM[rowShift+k] = colind;
      
      int shift=nt+(nx-2)*nt; 
      jcolM[colind]   = shift+1+k;         M[colind  ] = -K/dx;
      jcolM[colind+1] = shift+nt+1+k;      M[colind+1] = alphapp+K/dx;
    }
    rowShift+=nt-1; colShift+=2*(nt-1);
    /////////////////////////////////////////////////////////////////
    // one row                -- 5 -- 
    /////////////////////////////////////////////////////////////////
    krowM[rowShift+0] = colShift;
    jcolM[colShift+0] = 0; M[colShift+0] = 1.0;
    rowShift++; colShift++;

    /////////////////////////////////////////////////////////////////
    // last rows              -- 6 --
    /////////////////////////////////////////////////////////////////
    for(int k=0; k<nx; k++) {
      int colind = colShift + k;
      krowM[rowShift+k] = colind;
      jcolM[colind] = nt + k*nt; M[colind] = 1.0;
    }
    assert(rowShift+nx == nt*(nx+1));
    assert(colShift+nx == (nt-1)*4*nx + nx+1);
    krowM[rowShift+nx] = colShift+nx;
    
  }
  return 0;
}

extern "C"
int fC(void* user_data, int id, int* krowM, int* jcolM, double* M)
{ return 0; }

extern "C"
int fD(void* user_data, int id, int* krowM, int* jcolM, double* M)
{ return 0; }


extern "C"
int fb(void* user_data, int id, double* vec, int len)
{
  if(id==0) {
    //no equalities
  } else {
    MyNodeData* mydata = (MyNodeData*)user_data;
    int nx=mydata->nx; int nt = mydata->nt;
    
    // nt(nx+1)-by-nt(nx+1) matrix

    ///////////////////////////////////////////
    //first nt-1       -- 1 --
    ///////////////////////////////////////////
    double dt=1.0;
    for(int i=0; i<nt-1; i++) {
      vec[i] = 20000.0;
    }
    int rowShift=nt-1;
    ///////////////////////////////////////
    // rows  nt to (nt-1)(nx-2)  -- 2 --
    ///////////////////////////////////////
    for(int j=0; j<nx-2; j++) {
      for(int k=0; k<nt-1; k++) {
	vec[rowShift+j*(nt-1)+k] = 0.0;
      }
    }
    rowShift+=(nt-1)*(nx-2);
    ///////////////////////////////////////////////////////////////////
    // rows nt+(nt-1)(nx-2)  to   nt-1 + (nt-1)(nx-2) + nt-1   -- 3 --
    ///////////////////////////////////////////////////////////////////
    for(int k=0; k<nt-1; k++) {
      vec[rowShift+k] = 0.0;
    }
    rowShift+=nt-1; 
    ///////////////////////////////////////////////////////////////////
    // rows nt-1 + (nt-1)(nx-2) + nt    to 
    //    nt-1 + (nt-1)(nx-2) + nt + nt -1   --- 4 --- 
    ///////////////////////////////////////////////////////////////////
    for(int k=0; k<nt-1; k++) vec[rowShift+k] = 0.0;
    loadWeatherData("samples_weather_model.dat", id-1, nt-1, &vec[rowShift]);
    //loadWeatherData("samples_weather_model.dat", 50, nt-1, &vec[rowShift]);    
    for(int k=0; k<nt-1; k++) vec[rowShift+k] *= alphapp;
    
    rowShift+=nt-1; 
    /////////////////////////////////////////////////////////////////
    // one row                -- 5 -- 
    /////////////////////////////////////////////////////////////////
    vec[rowShift+0] = 25.0; //TI_ini;
    rowShift++; 

    /////////////////////////////////////////////////////////////////
    // last rows              -- 6 --
    /////////////////////////////////////////////////////////////////
    for(int k=0; k<nx; k++) {
      vec[rowShift+k] = 25.0;//TW_ini[k];
    }
    assert(rowShift+nx == nt*(nx+1));
  }
  return 0;
}

extern "C"
int fclow (void* user_data, int id, double* vec, int len)
{  return 0; }

extern "C"
int ficlow(void* user_data, int id, double* vec, int len)
{ return 0; }

extern "C"
int fcupp (void* user_data, int id, double* vec, int len)
{ return 0; }


extern "C"
int ficupp(void* user_data, int id, double* vec, int len)
{ return 0; }


extern "C"
int fxlow (void* user_data, int id, double* vec, int len)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  int nx=mydata->nx; int nt = mydata->nt;

  if(id==0) {
    assert(len==3*(nt-1));
    for(int i=0; i<len; i++) vec[i]= 0.0;
  } else {
    assert(len==nt*(nx+1));
    for(int i=0;  i<nt;        i++) vec[i] =  20.0;
    for(int i=nt; i<nt*(nx+1); i++) vec[i] = -20.0;
  }
  return 0;
}


extern "C"
int fixlow(void* user_data, int id, double* vec, int len)
{
  for(int i=0; i<len; i++) vec[i]=1.0;
  return 0;
}
extern "C"
int fixlowc(void* user_data, int id, char* vec, int len)
{
  for(int i=0; i<len; i++) vec[i]=1;
  return 0;
}


extern "C"
int fxupp (void* user_data, int id, double* vec, int len)
{
  MyNodeData* mydata = (MyNodeData*)user_data;
  int nx=mydata->nx; int nt = mydata->nt;

  if(id==0) {
    assert(len==3*(nt-1));
    for(int i=0;      i<nt-1;   i++) vec[i]= 100000.0;
    for(int i=nt-1;   i<2*nt-2; i++) vec[i]= 100000.0;
    for(int i=2*nt-2; i<3*nt-3; i++) vec[i]=  30000.0;
  } else {
    assert(len==nt*(nx+1));
    for(int i=0;  i<nt;        i++) vec[i] = 25.0;
    for(int i=nt; i<nt*(nx+1); i++) vec[i] = 50.0;
  }
  return 0;
}

extern "C"
int fixupp(void* user_data, int id, double* vec, int len)
{
  for(int i=0; i<len; i++) vec[i]=1.0;
  return 0;
}
extern "C"
int fixuppc(void* user_data, int id, char* vec, int len)
{
  for(int i=0; i<len; i++) vec[i]=1;
  return 0;
}

void loadWeatherData(const char* szFile, int scenario, int tSize, double* data)
{
  string strLine, strFile; int nLine=0;
  strFile = g_strDataDir+szFile;
  fstream file_op(strFile.c_str(), ios::in);
  while(!file_op.eof()) {
    getline(file_op, strLine);
    if(nLine==scenario) break;
    nLine++;
  }
  file_op.close();
  assert(strLine.length()>0);

  istringstream sstrLine(strLine); strLine.clear();
  double temper;   int i=0;
  //dump the first temperature <- not needed
  sstrLine >> temper;

  while( sstrLine >> temper && i<tSize) data[i++]=temper;
  assert(i==tSize);
}

void loadElectricityPrices(const char* szFile, int tSize, double* data)
{
  int nLine=0; double val;

  string strFile = g_strDataDir+szFile;
  fstream file_op(strFile.c_str(), ios::in);
  assert(file_op.is_open());
  while(!file_op.eof()) {
    file_op >> val;                       //time
    assert(val==nLine+1.0);
    file_op >> val;                       //TAV
    file_op >> data[nLine];               //price
    
    nLine++;
    if(nLine>=tSize) break;    
  }
  file_op.close();
}
