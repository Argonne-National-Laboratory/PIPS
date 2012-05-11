#include "StochInputTree.h"
#include "QpGenDriver.h"
#include "QpGenStochDriver.h"
#include "MehrotraSolver.h"
#include "GondzioSolver.h"
#include "QpGenSparseMa57.h"
#include "QpGenSparseMa27.h"
#include "MehrotraStochSolver.h"
#include "sFactoryAug.h"
//#include "sFactoryAugSca.h"
//#include "sFactoryAugEmtl.h"
//#include "sFactoryAugEmtlSym.h"
#include "QpGenStochAugExt.h"
#include "QpGenStochAugRed.h"
#include "QpGenStochAugRedPr.h"
#include "QpGenStochNrmEqn.h"
#include "QpStochAugRedPrCG.h"
#include "QpStochAugRedPrPCG.h"


#include "StochResourcePlanner.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <errno.h>
using namespace std;

#include "mpi.h"

#define SCALE_OBJ 100.0 //1000.0
#define SCALE_C   1000.0 //500.0 //5000.0

/** Problem parameters and data */
struct PbData 
{
public:
  int N,T,NW,NWT,ND;
  int S;
  int idx,num_days;
  double *Pmax, *Pmin;
  double *SU, *SD;
  double *RU, *RD;
  int *UT, *DT;
  int *v0;
  double *a, *b, *c;
  double *hc, *cc, *tcold, *C;
  int *initstate;

  double *demand, *R, *Rp;

  int *S0, *U0;
  int *G, *L;
  double **K;
  
public: //methods
  PbData(const string& strFile);
  virtual ~PbData();

  void compute1stStageSizes(int& n, int& my, int&mz);
  void compute2ndStageSizes(int& n, int& my, int&mz);

  int get_nx(int stage);
  int get_my(int stage);
  int get_mz(int stage);

  void loadWindScenario(int scenario, int tmStep, 
			double* data, int noWindUnits);
#ifdef TIMING
  static int OVERunits, OVERscen;
#endif
protected:
  void loadFromFile(const string& strFile);
  void paramFromLine(const string& strLine, double* param);
  void paramFromLine(const string& strLine, int* param);
  int getNextNonEmptyLine(fstream& file, string& line);

  double*** wind_scen;
  int scen_loaded;

  PbData() {};
};

static int myRank = 0;

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
int fxupp  (void* user_data, int id, double* vec, int len);
int fixupp (void* user_data, int id, double* vec, int len);
};

/** Declarations of utilities
 */
//string g_strDataDir = "../../apps/unitcommitment/";
#ifndef __bg__
string g_strDataDir = "/homes/mlubin/stoch/apps/unitcommitment/";
#else
string g_strDataDir = "/home/mlubin/STOCHPROG/apps/unitcommitment/";
//string g_strDataDir = "/pvfs-surveyor/mlubin/STOCHPROG-apps/unitcommitment/";
#endif

void stochSolve(int argc, char *argv[], PbData& pbData, int printx, int solveWith);
void qpgenSolve(PbData& data, int printx);
/** 
 * Usage: building.exe -solve 0 -ns 100 -nt 100 -nx 100 
 */
// config override variables
#ifdef TIMING
int PbData::OVERunits = -1, PbData::OVERscen = -1;
#endif

//#include <TAU.h>

int main( int argc, char *argv[] )
{ 
  //initialize MPI env
  MPI_Init(&argc, &argv);
  int mype; MPI_Comm_rank(MPI_COMM_WORLD, &mype);
 
  //TAU_TRACK_MEMORY();

  int goodUsage = 1; // Assume good. Why not be optimistic
  int iargv    = 1;
  int solveWith = -1, printx=0, scaprocs=-1;
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
    } else  if( 0 == strcmp( "scaprocs", &argv[iargv][iopt] ) ) {

      iargv++;
      if(iargv>=argc) break;

      char* endptr;
      scaprocs = strtol( argv[iargv], &endptr, 10 );
      if( '\0' != *endptr ) { 
        goodUsage = 0;assert(0);
        break;
      }
    }
    else  if( 0 == strcmp( "printx", &argv[iargv][iopt] ) ) {
      printx=1;
    } else if( 0 == strcmp("datadir", &argv[iargv][iopt]) ) {

      iargv++;
      if(iargv>=argc) break;

      char* endptr;
      g_strDataDir = argv[iargv];//strtol( argv[iargv], &endptr, 10 );
      if('/' != g_strDataDir[g_strDataDir.length()-1]) g_strDataDir += "/";
#ifdef TIMING
    } else if ( 0 == strcmp("units", &argv[iargv][iopt]) ) {
      iargv++;
      if(iargv>=argc) break;

      char* endptr;
      PbData::OVERunits = strtol( argv[iargv], &endptr, 10 );
      if( '\0' != *endptr ) { 
      	goodUsage = 0;assert(0);
      	break;
      }
    } else if ( 0 == strcmp("scen", &argv[iargv][iopt]) ) {
      iargv++;
      if(iargv>=argc) break;

      char* endptr;
      PbData::OVERscen = strtol( argv[iargv], &endptr, 10 );
      if( '\0' != *endptr ) { 
      	goodUsage = 0;assert(0);
      	break;
      }
#endif
    } else {
      goodUsage=0;
      break;
    }
    iargv++;
  } // end while
  
  // find a better way to handle this
  StochResourcePlanner::noScaProcesses = scaprocs;
  //printf("ARGS PROC %d, %d %d %d %d\n", mype, solveWith, scaprocs, PbData::OVERunits, PbData::OVERscen);

  if( argc>1 && iargv == argc && goodUsage ) {

    if(solveWith!=0 && 
       solveWith!=1 && solveWith!=2 && solveWith!=3 &&
       solveWith!=4 && solveWith!=5 && solveWith!=6 &&
       solveWith!=11 && solveWith!=12 && solveWith!=13 && 
       solveWith!=20 ) {
      cerr << "Solver(s) is not recognized/specified.\n";
      return 1;
    }
  } else {
    cerr << "Usage: " << argv[0] << "  --solve num [--datadir path] [--printx] [--scaprocs num]\n";
    cerr << "\t--solve num    - solver used: \n";
    cerr << "\t\t\t -  0 OOQP's default\n";
    cerr << "\t\t\t -  1 STOCH augmented (inequalities reduced)\n";
    cerr << "\t\t\t -  2 STOCH augmented (inequalities reduced) memory-friendly\n";
    cerr << "\t\t\t -  3 STOCH augmented\n";
    cerr << "\t\t\t - 11 STOCH reduced-augmented with BICGSTAB\n";
    cerr << "\t\t\t - 12 STOCH reduced-augmented with Projected CG\n";
    cerr << "\t\t\t - 13 STOCH reduced-augmented with CG (root node should have no eq. constraints)\n";
    cerr << "\t\t\t - 20 STOCH normal equations.\n";
    cerr << "\t--datadir path - path to the directory containing data of Unit Commitment problem\n";
    cerr << "\t--printx       - print the 'x' portion of the solution.\n";
    cerr << "\t--scaprocs num - number of processors to use with scalapack, can't be more than total.\n";
    cerr <<"\n";
    
    return 1;
  }


  string strFilename = "UnitCommitment.dat";
  PbData params(strFilename);
  if(solveWith>0) {

    //stochastic solve
    stochSolve(argc, argv, params, printx, solveWith);
  } else {
    //QpGen solve
    qpgenSolve(params, printx);
  }

  MPI_Finalize();

  return 0;
}



void stochSolve(int argc, char *argv[], PbData& pbData, int printx, int solveWith)
{

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  myRank = rank;
  int size; MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(size<2 && (solveWith==11 || solveWith==12 || solveWith==13) ) {
    cerr << "Could not employ preconditioner since it needs at least 2 processes.\n";
    return;
  }


  //generate input stochastic(scenario) tree
  int globalID=0; 
  int n,my,mz,nx0,my0,mz0;
  pbData.compute1stStageSizes(nx0,my0,mz0);

  StochInputTree::StochInputNode data(&pbData, globalID, 
			   nx0,my0,mz0,
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
  
  for(int i=0; i<pbData.S; i++) {
    pbData.compute2ndStageSizes(n,my,mz);

    StochInputTree::StochInputNode data(&pbData, globalID, 
			     n, my, mz, 
			     fQ, fnnzQ, fc,
			     fA, fnnzA,
			     fB, fnnzB,
			     fb,
			     fC, fnnzC,
			     fD, fnnzD,
			     fclow, ficlow, fcupp, ficupp,
			     fxlow, fixlow, fxupp, fixupp );
    globalID++;
    
    root->AddChild(new StochInputTree(data));
  }
  if(myRank==0) {
#ifdef TIMING
    char *match = strstr(argv[0],"mkl");
    if (match != 0) {
      printf("VERSION MKL\n");
    } else {
      printf("VERSION NETLIB\n");
    }
    printf("ALL NPROCS %d\n", size);
    printf("UNITS %d T %d SCEN %d  ND %d  NW %d NWT %d\n", pbData.N,pbData.T,pbData.S,
	   pbData.ND,pbData.NW,pbData.NWT);
    printf("L1 VAR %4d EQ %4d INEQ %4d\n", nx0, my0, mz0);
    printf("L2 VAR %4d EQ %4d INEQ %4d\n", n, my, mz);
#else
    printf("%d processes are used.\n", size);
   

    printf("N=%d  T=%d S=%d  ND=%d  NW=%d NWT=%d\n", pbData.N,pbData.T,pbData.S,
	   pbData.ND,pbData.NW,pbData.NWT);

    printf("Level 1: %4d variables, %4d eq constraints, %5d ineq constr\n",
	   nx0, my0, mz0);
    printf("Level 2: %4d variables, %4d eq constraints, %5d ineq constr\n",
	   n, my, mz);
#endif
  }
  
  //call the stochastic IPM solver
  StochRunParams* params = defaultStochRunParams();
  params->printx = printx;
  

  MehrotraStochSolver* method=NULL;
  if(solveWith==3) {
    QpGenStochAugExt* formulation=NULL;
    qpgenstoch_solve( root, params, method, formulation);
  } else if(solveWith==11) {
    QpGenStochAugRedPr* formulation=NULL;
    qpgenstoch_solve( root, params, method, formulation);
  } else if(solveWith==12) {
    QpStochAugRedPrPCG* formulation=NULL;
    qpgenstoch_solve( root, params, method, formulation);
  } else if(solveWith==13) {
    QpStochAugRedPrCG* formulation=NULL;
    qpgenstoch_solve( root, params, method, formulation);
  } else if(solveWith==20) {
    QpGenStochNrmEqn* formulation=NULL;
    qpgenstoch_solve( root, params, method, formulation);
  } else if(solveWith==1) {
    QpGenStochAugRed* formulation=NULL;
    qpgenstoch_solve( root, params, method, formulation);
  } else if(solveWith==2) {
    sFactoryAug* formulation=NULL;
    qpgenstoch_solve( root, params, method, formulation);
  } /*else if(solveWith==4) {
    sFactoryAugSca* formulation=NULL;
    qpgenstoch_solve( root, params, method, formulation);
  } else if(solveWith==5) {
    sFactoryAugEmtl* formulation=NULL;
    qpgenstoch_solve( root, params, method, formulation);
  } else if(solveWith==6) {
    sFactoryAugEmtlSym* formulation=NULL;
    qpgenstoch_solve( root, params, method, formulation);
    }*/
  
  
  delete params;
  delete root;


}


extern "C"
int fnnzQ(void* user_data, int id, int* nnz)
{ *nnz = 0; }

extern "C"
int fnnzA(void* user_data, int id, int* nnz)
{ 
  PbData* data = (PbData*)user_data;
  int N=data->N; int ND=data->ND; int T=data->T; int S=data->S; int NW=data->NW;

  if(id==0) {
    *nnz = 0;
    //!2 *nnz = 3*N; //1 
    //for 16 & 19
    for(int j=0; j<N; j++) 
      *nnz += (data->G[j]+data->L[j]);

  } else { //2nd stage -> eq matrix for 1st stage variables
    *nnz = 0;
    //!2 *nnz = N*(T-1);        // 1
  }
}

extern "C"
int fnnzB(void* user_data, int id, int* nnz)
{ 
  if(id==0) {
    //nothing in B at 1st stage
    *nnz = 0;
  } else {

    *nnz = 0;
    //!2
    //PbData* data = (PbData*)user_data;
    //int N=data->N; int T=data->T;
    //*nnz = 2*N*(T-1);        // 1
  }
}

extern "C"
int fnnzC(void* user_data, int id, int* nnz)
{ 
  PbData* data = (PbData*)user_data;
  int N=data->N; int ND=data->ND; int T=data->T; int NW=data->NW; int S=data->S;
  
  if(id==0) {

    // 2
    *nnz = N*ND*T+N*ND*T + N*((2*T-1)*ND*(ND+1)/2 - ND*(ND+1)*(2*ND+1)/6) / 2;

    // 3 
    *nnz += (N*(2*T-1)+N*T);
    //    printf("3nnz=%d\n", *nnz);
    // 5 & 6
    //*nnz += (N+NW + N+NW);
    *nnz += (N+NW + N+NW);
    
    // 8, 9 & 10
    //*nnz += (2*N + 2*N + 2*N);
    *nnz += (2*N + 2*N + 2*N);

    // 17
    for(int j=0; j<N; j++) {
      *nnz += ( (1+data->UT[j]) * (T - data->UT[j] - data->G[j] + 1) );
      if(data->G[j]==0) {(*nnz)--;}

      assert(T - data->UT[j] - data->G[j] + 1 >=0);
    }

    // 18  
    for(int j=0; j<N; j++) {
      assert(data->T-data->UT[j]+2>=2);
      
      for(int k=data->T-data->UT[j]+2; k<=data->T; k++) {
	*nnz += ( data->T - k + 2 );
      }
    }

    // 20
    for(int j=0; j<N; j++) {
      *nnz += ( (1+data->DT[j]) * (T - data->DT[j] - data->L[j] + 1) );
      if(data->L[j]==0) {(*nnz)--;}

      assert(T - data->DT[j] - data->L[j] + 1 >=0);
    }

    // 21
    for(int j=0; j<N; j++) {
      assert(data->T-data->DT[j]+2>=2);
      for(int k=data->T-data->DT[j]+2; k<=data->T; k++)
	*nnz += ( data->T - k + 2 );
    }
  } else { //2nd stage -> ineq matrix for 1st stage variables
    // 5 -> nothing 6 -> nothing
    *nnz = N*(T-1);         // 8 
    // 9->nothing
    *nnz += N*(T-1);         // 10
    *nnz += 2*N*(T-1)+N;     // 11
    *nnz += 2*N*(T-1)+N;     // 12
  }
}

extern "C"
int fnnzD(void* user_data, int id, int* nnz)
{ 
  PbData* data = (PbData*)user_data;
  int N=data->N; int ND=data->ND; int T=data->T; int NW=data->NW;

  if(id==0) {
    *nnz = 0;
  } else {
    *nnz =  (N+NW)*(T-1);    // 5
    *nnz += (N+NW)*(T-1);    // 6
    *nnz += N*(T-1);         // 8
    *nnz += 2*N*(T-1);       // 9
    *nnz += N*(T-1);         //10
    *nnz += N*(T-1)+N*(T-2); //11
    *nnz += N*(T-1)+N*(T-2); //12
  }
  
}

extern "C"
int fQ(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  PbData* data = (PbData*)user_data;
  if(id==0) {
    int n = data->get_nx(1);
    for(int i=0; i<=n;  i++) krowM[i]=0;
  } else {
    int n = data->get_nx(2);
    for(int i=0; i<=n; i++) krowM[i]=0;
  }
  return 0;
}

extern "C"
int fc(void* user_data, int id, double* vec, int len)
{
  PbData* data = (PbData*)user_data;
  int N=data->N; int T=data->T; int NW=data->NW; int S=data->S;
  double* a = data->a; double* b = data->b;
  int NT=N*T; 

  if(id==0) {
    for(int k=0; k<T; k++) {
      for(int j=0; j<N; j++) {
	vec[k*N+j]          = a[j];   // v
	vec[k*N+j + N*T]    = 1.0;    //cu
	vec[k*N+j + 2*NT]   = 1.0;    //cd
      }
    }
    for(int j=0; j<N; j++) {
      vec[3*NT + j]    = b[j];   // pc[j,1,s] 
      vec[3*NT+N +j]   = 0.0;   // pcmax
    }
 
    assert(3*NT+2*N+NW==len); 
    for(int i=3*NT+2*N; i<3*NT+2*N+NW; i++)
      vec[i]      = 0.0;   // pwind[j,1]
  } else {

    for(int k=0; k<T-1; k++) {
      for(int j=0; j<N; j++) {
	vec[k*N+j]      = b[j]/S;    //pc
	//!2 vec[i+NT-N] = 1.0/S;  //cp
      }
    }
    for(int i=N*(T-1); i<N*(T-1)+NW*(T-1); i++)
      vec[i]      = 0.0;   //pwind

    assert(N*(T-1)+NW*(T-1)+N*(T-1) == len);
    for(int i=(N+NW)*(T-1); i<(2*N+NW)*(T-1); i++)
      vec[i]      = 0.0;   //pcmax
  }  
  for(int i=0; i<len; i++) vec[i] /= SCALE_OBJ;
}

extern "C"
int fA(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  PbData* data = (PbData*)user_data;
  int N=data->N; int T=data->T; int NW=data->NW; int S=data->S;
  int NT = N*T;
  double* a = data->a; double* b = data->b;
  int* G = data->G; int* L = data->L;

  if(id==0) {
    //1st stage

    // 1 -> N  x  3*NT + 3*N + NW
    /*
    //!2
    krowM[0]=0;
    for(int j=0; j<N; j++) {
      jcolM[3*j]   = j;        M[3*j]   = a[j];    //A
      jcolM[3*j+1] = j+3*NT;   M[3*j+1] = b[j];    //B
      jcolM[3*j+2] = j+3*NT+N; M[3*j+2] = -1.0;    //C

      krowM[j+1] = 3*j+3; 
    }
    int ind_krowM=N; int ind_jcolM=3*N; assert(ind_jcolM==krowM[ind_krowM]);
    */
    int ind_krowM=0; int ind_jcolM=0;
    //16
    for(int j=0; j<N; j++) {
      if(G[j]==0) continue;
      for(int c=0; c<G[j]; c++) {
	jcolM[ind_jcolM+c] = c*N+j; M[ind_jcolM+c] = 1.0; //A
      }
      ind_jcolM += G[j];

      krowM[ind_krowM+1] = ind_jcolM;
      ind_krowM++;
    }


    //19
    for(int j=0; j<N; j++) {
      if(L[j]==0) continue;
      for(int c=0; c<L[j]; c++) {
	jcolM[ind_jcolM+c] = c*N*j; M[ind_jcolM+c] = 1.0; //A
      }
      ind_jcolM += L[j];
      krowM[ind_krowM+1] = ind_jcolM;
      ind_krowM++;
    }

    if(ind_krowM==0) {
      int my = data->get_my(1); 
      for(int i=0; i<=my; i++) {
	krowM[i]=0;
      }
    }

#ifdef DEBUG
    int my = data->get_my(1); 
    assert(my==ind_krowM);

    int nnz; fnnzA(data, 0, &nnz);
    assert(nnz==krowM[ind_krowM]);
    
#endif
  } else {

    //////////////////////////////////////////
    //2nd stage - for 1st stage vars
    //////////////////////////////////////////
    // N*(T-1) equalities
    int colind=0, rowind=0; krowM[rowind] = colind;

    // ------- 1 --------
    /* //!2
    for(int k=0; k<T-1; k++) {
      for(int j=0; j<N; j++) {
	
	jcolM[colind] = (k+1)*N+j; //!1
	M[colind]    = a[j];
	colind++;
	krowM[rowind+1] = colind; rowind++;
      }
    }
    */

    if(rowind==0) {
      int my = data->get_my(2); 
      for(int i=0; i<=my; i++) {
	krowM[i]=0;
      }
    }
#ifdef DEBUG  
    int my = data->get_my(2); 
    assert(rowind==my);
    int nnz; fnnzA(data, id, &nnz);
    assert(colind==nnz);
#endif

  }
  return 0;
}

extern "C"
int fB(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
  PbData* data = (PbData*)user_data;
  int N=data->N; int T=data->T; int NW=data->NW;
 
  if(id==0) {
    //empty matrix 
    assert(false);

  } else {
    //////////////////////////////////////////////////
    // 2nd stage -> 2nd stage variables
    //////////////////////////////////////////////////
    int colind=0, rowind=0; krowM[rowind] = colind;

    // ------ 1 ------
    /* //!2
    double* b = data->b;
    for(int k=0; k<T-1; k++) {
      for(int j=0; j<N; j++) {
	jcolM[colind] = k*N+j;           M[colind] = b[j]; colind++;
	jcolM[colind] = k*N+j + N*(T-1); M[colind] =-1.0;  colind++;
	krowM[rowind+1] = colind; rowind++;
      }
    }
    */
    if(rowind==0) {
      int my = data->get_my(2); 
      for(int i=0; i<=my; i++) {
	krowM[i]=0;
      }
    }
#ifdef DEBUG  
    int my = data->get_my(2); 
    assert(rowind==my);
    int nnz; fnnzB(data, id, &nnz);
    assert(colind==nnz);
#endif

  }

  return 0;
}

extern "C"
int fC(void* user_data, int id, int* krowM, int* jcolM, double* M)
{ 
  PbData* data = (PbData*)user_data;
  int N=data->N; int T=data->T; int NW=data->NW; int ND=data->ND; int S = data->S;
  double** K = data->K;
  int NT=N*T;
  int nnz=0;

  if(id==0) {
    ////////////////////////////////////
    //   ------ 1st stage ---------
    ////////////////////////////////////

    int checkcolind=0;

    int rowind=0; int colind=0;
    krowM[rowind]=0;
    /////////////////////////
    // -------2 ------- 
    /////////////////////////
    for(int t=1; t<=ND; t++) {
      // for each t we have T*N equations = T blocks of N eqn
      
      //first (t+1) blocks
      
      for(int k=0; k<t+1; k++) {

	for(int line=0; line<N; line++) {
	  krowM[rowind] = colind;

	  //first k-1 blocks, k=0,..,t
	  for(int block=0; block<k; block++) {
	    jcolM[colind] = block*N+line;
	    M[colind] = K[line][t-1];
	    colind++;
	  }
	  //block k
	  jcolM[colind] = k*N+line;
	  M[colind] = -K[line][t-1];
	  colind++;

	  //entry for cu <- identity NT x NT
	  jcolM[colind] = N*T + k*N+line; assert(N*T + k*N+line<2*N*T);
	  M[colind] = 1.0;
	  colind++;

	  krowM[rowind+1] = colind;
	  rowind++;
	} // end line=0,..,N-1
      } // end k=0,t-1

      assert(rowind==(t-1)*N*T + (t+1)*N);
      assert(colind==checkcolind+N*(t+1)*(t+2)/2+(t+1)*N);
 
      //next (last) (T-t-1) blocks
      for(int k=t+1; k<T; k++) {
	for(int line=0; line<N; line++) {
	  krowM[rowind] = colind;

	  //first t elements of the line 'line' in the k-th block
	  for(int block=0; block<t; block++) {
	    jcolM[colind] = block*N+line + (k-t)*N;
	    M[colind]     = K[line][t-1];
	    colind++;
	  }
	  // last block ( (t+1)-th )
	  jcolM[colind] = t*N+line + (k-t)*N; 
	  M[colind] = -K[line][t-1];
	  colind++;

	  //entry for cu <- identity NT x NT shifted NT columns to the right
	  jcolM[colind] = N*T + k*N+line; //assert(N*T + k*N+line<2*N*T);
	  M[colind] = 1.0;
	  colind++;

	  krowM[rowind+1] = colind;
	  rowind++;
	} // end line=0,..,N-1
      } // end for k=t, T
      assert(colind==checkcolind+N*(t+1)*(t+2)/2+(t+1)*N  +  (t+1)*N*(T-t-1)+(T-t-1)*N); 
      checkcolind = colind;
      assert(rowind==t*N*T);

      //printf("fC t=%d  --- done\n", t);
    } // end outer 'for' t=0,...,ND-1

#ifdef DEBUG
    nnz = 2*N*ND*T + N*( (2*T-1)*ND*(ND+1)/2 - ND*(ND+1)*(2*ND+1)/6 ) / 2;
    for(int j=0; j<nnz; j++) if(jcolM[j]>=2*N*T) assert(false);
    assert(colind==nnz); int nnz0=nnz;
#endif

#ifdef DEBUG
    int row_save=rowind; int col_save=colind;
#endif

    /////////////////////////
    // ------- 3 ---------
    /////////////////////////
    double* C = data->C;
    for(int k=0; k<T; k++) {
      for(int j=0; j<N; j++) {
	krowM[rowind] = colind;
	
	if(k>0) {
	  jcolM[colind] = (k-1)*N+j;
	  M[colind] = -C[j];
	  colind++;
	}
	jcolM[colind] = k*N+j;
	M[colind] = C[j]; 
	colind++;

	//identity block corresponding to 'cd'
	jcolM[colind] = k*N+j + 2*N*T;
	M[colind] = 1.0;
	colind++;

	krowM[rowind+1] = colind;
	rowind++;
      }
    }

#ifdef DEBUG
    assert(rowind-row_save==T*N);
    row_save = rowind;
    nnz = N*(3*T-1);
    assert(colind-col_save==nnz);
    col_save = colind;
    for(int i=nnz0; i<nnz; i++) if(jcolM[i]>=N*T && jcolM[i]<2*N*T) assert(false);
    nnz0=nnz;
#endif   

    ////////////////////////
    // ------- 5 --------
    ////////////////////////
    krowM[rowind] = colind;
    for(int i=0; i<N; i++) {
      jcolM[colind] = 3*NT+i;
      M[colind] = 1.0;
      colind++;
    }
    for(int i=0; i<NW; i++) {
      jcolM[colind] = 3*NT+2*N+i; //!2
      M[colind] = -data->NWT;
      colind++;
    }
    krowM[rowind+1] = colind;
    rowind++;


    // ------- 6 ------- 
    for(int i=0; i<N; i++) {
      jcolM[colind] = 3*NT+N+i;//!2
      M[colind] = 1.0;
      colind++;
    }
    for(int i=0; i<NW; i++) {
      jcolM[colind] = 3*NT+2*N+i;//!2
      M[colind] = -data->NWT;
      colind++;
    }
    krowM[rowind+1] = colind;
    rowind++;

#ifdef DEBUG
    nnz = (N+NW + N+NW);
    assert(rowind-row_save==2);
    row_save = rowind;
    assert(colind-col_save==nnz);
    col_save = colind;
#endif

    // 8, 9 & 10
    //*nnz += S*(2*N + 2*N + 2*N);

    // ------ 8 ------
    for(int i=0; i<N; i++) {
      jcolM[colind] = i;
      M[colind] = -data->Pmin[i];
      colind++;

      jcolM[colind] = i+3*NT;
      M[colind] = 1.0;
      colind++;

      krowM[rowind+1] = colind; rowind++;
    }
    // ------- 9 ------ 
    for(int i=0; i<N; i++) {
      jcolM[colind] = 3*NT+i;
      M[colind] = -1.0; //pc
      colind++;

      jcolM[colind] = i+3*NT+N; //!2
      M[colind] = 1.0; //pcmax
      colind++; 

      krowM[rowind+1] = colind; rowind++;
    }

    // ------ 10 -------
    for(int i=0; i<N; i++) {
      jcolM[colind] = i;
      M[colind] = data->Pmax[i];
      colind++;

      jcolM[colind] = i+3*NT+N;//!2
      M[colind] = -1.0;
      colind++;

      krowM[rowind+1] = colind; rowind++;
    }

#ifdef DEBUG
    nnz = (6*N);
    assert(rowind-row_save==3*N);
    row_save = rowind;
    assert(colind-col_save==nnz);
    col_save = colind;
    fnnzC(data, id, &nnz);
#endif   
    
    // ------- 17 -------  (T - UT[j] - G[j] + 1)
    int* UT = data->UT; int* G = data->G;
    for(int j=0; j<N; j++) {
      //the below implem should work for G[j]>0 but is not tested.
      assert(G[j]==0);

      for(int k=G[j]; k<T-UT[j]+1; k++) {
	assert(k<T);
	//[j,k-1]
	if(k>0) {
	  jcolM[colind] = (k-1)*N+j;
	  M[colind] = UT[j];
	  colind++;
	}
	//[j,k]
	jcolM[colind] = k*N+j;
	M[colind] = 1.0-UT[j];
	colind++;
	//j,k+1 .. k+UT[j]-1]
	for(int i=k+1; i<=k+UT[j]-1; i++) {
	  jcolM[colind] = i*N+j;
	  M[colind] = 1.0;
	  colind++;
	}

	krowM[rowind+1] = colind; rowind++;
#ifdef DEBUG
	assert(colind<nnz);	
#endif
      }
    }
    
    // ------ 18 ------
    for(int j=0; j<N; j++) {

      for(int k=T-UT[j]+1; k<T; k++) {

	assert(k>0);
	//[j,k-1]
	jcolM[colind] = (k-1)*N+j;
	M[colind] = 1.0+T-(k+1); //!1
	colind++;
	//[j,k]
	jcolM[colind] = k*N+j;
	M[colind] = 0.0-T+(k+1); //!1
	colind++;
	//[j,..]
	for(int i=k+1; i<T; i++) {
	  jcolM[colind] = i*N+j;
	  M[colind] = 1.0;
	  colind++;
	}
	krowM[rowind+1] = colind; rowind++;
#ifdef DEBUG
	assert(colind<nnz);
#endif
      }
    }

    // ------- 20 -------
    int* DT = data->DT; int* L = data->L;
    for(int j=0; j<N; j++) {
     //the below implem should work for G[j]>0 but is not tested.
      assert(L[j]==0);
      for(int k=L[j]; k<T-DT[j]+1; k++) {
	assert(k<T);
	//[j,k-1]
	if(k>0) {
	  jcolM[colind] = (k-1)*N+j;
	  M[colind] = 0.0-DT[j];
	  colind++;
	}
	//[j,k]
	jcolM[colind] = k*N+j;
	M[colind] = -1.0+DT[j];
	colind++;
	//j,k+1 .. k+UT[j]-1]
	for(int i=k+1; i<=k+DT[j]-1; i++) {
	  jcolM[colind] = i*N+j;
	  M[colind] = -1.0;
	  colind++;
	}

	krowM[rowind+1] = colind; rowind++;
#ifdef DEBUG
	assert(colind<nnz);
#endif
      }
    }
    // ------- 21 -------
    for(int j=0; j<N; j++) {
      for(int k=T-DT[j]+1; k<T; k++) {
	assert(k>0);
	//[j,k-1]
	jcolM[colind] = (k-1)*N+j;
	M[colind] = 0.0-T+ (k+1) +1.0; //!1
	colind++;
	//[j,k]
	jcolM[colind] = k*N+j;
	M[colind] = 0.0+T-(k+1); //!1
	colind++;
	//[j,..]
	for(int i=k+1; i<T; i++) {
	  jcolM[colind] = i*N+j;
	  M[colind] = -1.0;
	  colind++;
	}
	krowM[rowind+1] = colind; rowind++;
#ifdef DEBUG
	assert(colind<=nnz);
#endif
      }
    }
    fnnzC(data, id, &nnz);

#ifdef DEBUG
    assert(colind==nnz);
    int mz = data->get_mz(1);
    assert(rowind==mz);

    //a whole new level of checking
    int nx = data->get_nx(1);
    for(int i=0; i<nnz; i++) {
      if(jcolM[i]>=nx) {
	assert(false);
      }
    }
#endif    
    for(int i=0; i<nnz; i++) M[i] /= SCALE_C;

  } else {
    ////////////////////////////////////
    // -------- 2nd stage ---------
    ////////////////////////////////////
    double* Pmin = data->Pmin; double* Pmax = data->Pmax;

    krowM[0]=0; int rowind=0; int colind=0;

    // ------- 5 & 6 -------
    // 1st stage variables are  not part of these equations 
    while(rowind<2*(T-1)) krowM[rowind++]=0;
    assert(rowind==2*(T-1));

    // ------- 8 -------
    for(int k=0; k<T-1; k++) {
      for(int j=0; j<N; j++) {
	jcolM[colind] = (k+1)*N+j; //!1
	M[colind]=-Pmin[j];
	colind++;

	krowM[rowind+1] = colind; rowind++;
      }
    }
    
#ifdef DEBUG
    assert(rowind==2*(T-1) + N*(T-1));
    assert(colind==N*(T-1)); 
#endif

    // ------ 9 ------
    for(int i=1; i<=N*(T-1); i++) {
      krowM[rowind+1] = colind; rowind++;
    }

    // ------- 10 -------
    for(int k=0; k<T-1; k++) {
      for(int j=0; j<N; j++) {
	jcolM[colind] = (k+1)*N+j;//!1
	M[colind] = Pmax[j];
	colind++;

	krowM[rowind+1] = colind; rowind++;
      }
    }
#ifdef DEBUG
    assert(rowind==2*(T-1) + 3*N*(T-1));
    assert(colind==2*N*(T-1)); 
#endif

    // ------ 11 -------
    double* RU = data->RU; double* SU = data->SU; 
    int ss = id-1;//assert(ss>=0 && ss <S);
    // --- k=2 
    for(int j=0; j<N; j++) {
      jcolM[colind] = j;       M[colind] = RU[j]-SU[j];   colind++; //v[j,1]
      jcolM[colind] = N+j;     M[colind] = SU[j]-Pmax[j]; colind++; //v[j,2]
      //jcolM[colind] = 3*N*T+j; M[colind] = 1.0;           colind++; 
      jcolM[colind] = 3*N*T +j; M[colind] = 1.0;           colind++; //pc[j,1,1] 
      krowM[rowind+1] = colind; rowind++;
    }
    // --- k>2 -> do it T-2 times. We start at k=1 to be compatible with the indices.
    for(int k=1; k<T-1; k++) {
      for(int j=0; j<N; j++) {
	jcolM[colind] = k*N+j;      M[colind] = RU[j]-SU[j];   colind++;//v[j,k-1]
	jcolM[colind] = (k+1)*N+j;  M[colind] = SU[j]-Pmax[j]; colind++;//v[j,k]

	krowM[rowind+1] = colind; rowind++;
      }
    }
#ifdef DEBUG
    assert(rowind==2*(T-1) + 3*N*(T-1) + N*(T-1));
    assert(colind==2*N*(T-1) + 2*N*(T-1)+N); 
#endif

    // ------- 12 -------
    double* RD = data->RD; double* SD = data->SD;
    // --- k=2 
    for(int j=0; j<N; j++) {
      jcolM[colind] = j;            M[colind] = SD[j]-Pmax[j]; colind++;
      jcolM[colind] = N+j;          M[colind] = RD[j]-SD[j];   colind++; 
      //jcolM[colind] = 3*N*T+N+NW+j; M[colind] = -1.0;          colind++; //pcmax[j,1,s]
      jcolM[colind] = 3*N*T+N + j; M[colind] = -1.0;          colind++; //pcmax[j,1,s]
      krowM[rowind+1] = colind; rowind++;
    }
    // --- k>2 -> do it T-2 times. We start at k=1 to be compatible with the indices.
    for(int k=1; k<T-1; k++) {
      for(int j=0; j<N; j++) {
	jcolM[colind] = k*N+j;         M[colind] = SD[j]-Pmax[j]; colind++;
	jcolM[colind] = (k+1)*N+j;     M[colind] = RD[j]-SD[j];   colind++; 
	krowM[rowind+1] = colind; rowind++;
      }
    }
    int nnz; fnnzC(data, id, &nnz);
#ifdef DEBUG
    assert(rowind==2*(T-1) + 3*N*(T-1) + N*(T-1) + N*(T-1));
    assert(colind==2*N*(T-1) + 2*N*(T-1)+N + 2*N*(T-1)+N); 

    int mz = data->get_mz(2);
    assert( mz==rowind);
    assert(nnz==colind);

    //a whole new level of checking
    int nx = data->get_nx(1);
    for(int i=0; i<nnz; i++) {
      if(jcolM[i]>=nx) {
	//printf("fC-2- jcolM[%d]=%d > %d\n", i, jcolM[i], nx);
	assert(false);
      }
    }

#endif
    for(int i=0; i<nnz; i++) M[i] /= SCALE_C;
  }

  return 0; 
}

extern "C"
int fD(void* user_data, int id, int* krowM, int* jcolM, double* M)
{ 
  PbData* data = (PbData*)user_data;
  int N=data->N; int T=data->T; int NW=data->NW;

  if(id==0) {assert(false);}
  else {
    int rowind = 0, colind = 0;
    krowM[rowind] = colind;
    ////////////////////////////////////
    //   ------ 2nd stage ---------
    ////////////////////////////////////
    // 2nd stage variables
    // ------- 5 ------
    double NWT = data->NWT;
    for(int k=0; k<T-1; k++) {
      for(int i=0; i<N; i++) {
	jcolM[colind] = k*N+i; 
	M[colind] = 1.0; 
	colind++;
      }
      for(int i=0; i<NW; i++) {
	jcolM[colind] = N*(T-1)+k*NW+i; //!2
	M[colind] =-NWT; 
	colind++;
      }
      krowM[rowind+1] = colind; rowind++;
    }
#ifdef DEBUG
    assert(rowind==T-1);
    assert(colind==(N+NW)*(T-1));
#endif
    // ------- 6 ------
    for(int k=0; k<T-1; k++) {
      for(int i=0; i<NW; i++) {
	jcolM[colind] = N*(T-1)+k*NW+i;//!2
	M[colind] =-NWT; 
	colind++;
      }

     for(int i=0; i<N; i++) {
       jcolM[colind] = (N+NW)*(T-1)+k*N+i; //pcmax //!2
	M[colind] = 1.0; 
	colind++;
      }
     krowM[rowind+1] = colind; rowind++;
    }
#ifdef DEBUG
    assert(rowind==T-1 + T-1);
    assert(colind==(N+NW)*(T-1) + (N+NW)*(T-1));
#endif
    // ------- 8 ------
    for(int k=0; k<T-1; k++) {
      for(int j=0; j<N; j++) {
	jcolM[colind] = k*N+j;
	M[colind] = 1.0; 
	colind++;

	krowM[rowind+1] = colind; rowind++;
      }
    }
#ifdef DEBUG
    assert(rowind==T-1 + T-1  + N*(T-1));
    assert(colind==(N+NW)*(T-1) + (N+NW)*(T-1) + N*(T-1));
#endif

    // ------- 9 ------
    for(int k=0; k<T-1; k++) {
      for(int j=0; j<N; j++) {
	jcolM[colind] = k*N+j;                M[colind] =-1.0; colind++;
	jcolM[colind] = k*N+j+(N+NW)*(T-1); M[colind] = 1.0; colind++; //!2
	krowM[rowind+1] = colind; rowind++;
      }
    }
#ifdef DEBUG
    assert(rowind==T-1 + T-1  + N*(T-1) + N*(T-1));
    assert(colind==(N+NW)*(T-1) + (N+NW)*(T-1) + N*(T-1) + 2*N*(T-1));
#endif
    // ------- 10 ------
    for(int k=0; k<T-1; k++) {
      for(int j=0; j<N; j++) {
	jcolM[colind] = k*N+j+(N+NW)*(T-1); M[colind] =-1.0; colind++; //!2
	krowM[rowind+1] = colind; rowind++;
      }

    }
#ifdef DEBUG
    assert(rowind==T-1 + T-1  + N*(T-1) + N*(T-1) + N*(T-1));
    assert(colind==(N+NW)*(T-1) + (N+NW)*(T-1) + N*(T-1) + 2*N*(T-1) + N*(T-1));
#endif
    // ------ 11 -------
    // --- k=2
    for(int j=0; j<N; j++) {
      // no B
      // now the C -> -identity
      jcolM[colind] = N*(T-1)+NW*(T-1) + j; M[colind] = -1.0; colind++; //!2
      krowM[rowind+1] = colind; rowind++;
    }
    // --- k>2
    for(int k=1; k<T-1; k++) {
      for(int j=0; j<N; j++) {
	//B -> identity
	jcolM[colind] = (k-1)*N+j;                M[colind] = 1.0; colind++;
	//C-> -identity
	jcolM[colind] = N*(T-1)+NW*(T-1) + k*N+j; M[colind] =-1.0; colind++; //!2
	krowM[rowind+1] = colind; rowind++;
      }
    }
#ifdef DEBUG
    assert(rowind==T-1 + T-1  + N*(T-1) + N*(T-1) + N*(T-1) + N*(T-1));
    assert(colind==(N+NW)*(T-1) + (N+NW)*(T-1) + N*(T-1) + 2*N*(T-1) + N*(T-1) 
	   + N*(T-1)+N*(T-2));
#endif
    // ------ 12 -------
    for(int k=0; k<T-1; k++) {
      for(int j=0; j<N; j++) {
	jcolM[colind] = k*N+j; M[colind] = 1.0; colind++;

	if(k>0) {
	  jcolM[colind] = N*(T-1)+NW*(T-1) + (k-1)*N+j;//!2
	  M[colind] =-1.0;
	  colind++;
	}
	krowM[rowind+1] = colind; rowind++;
      }
    }
    int nnz; fnnzD(data, id, &nnz);
#ifdef DEBUG
    int mz = data->get_mz(2);
    assert(rowind==mz);

    assert(colind==nnz);
    //a whole new level of checking
    int nx = data->get_nx(2);
    for(int i=0; i<nnz; i++) {
      if(jcolM[i]>=nx) {
	//printf("fC-2- jcolM[%d]=%d > %d\n", i, jcolM[i], nx);
	assert(false);
      }
    }
#endif

    for(int i=0; i<nnz; i++) M[i] /= SCALE_C;
  }
  return 0; 
}


extern "C"
int fb(void* user_data, int id, double* vec, int len)
{
  PbData* data = (PbData*)user_data;
  int N=data->N; int T=data->T; int NW=data->NW; int S=data->S;
  int* G = data->G; int* L = data->L;

  if(id==0) {
    assert(len==0);

    for(int i=0; i<len; i++) vec[i] = 0.0;

    // 16 & 19 not implemented.
   for(int j=0; j<N; j++) 
     if(G[j]>0 || L[j]>0) cerr << " feature not implemented (16 & 19 blocks)" << endl;

  } else {
    assert(len==0);
    // 2nd stage
    //!2    assert(len==N*(T-1));
    // 1
    //!2 for(int i=0; i<len; i++) vec[i] = 0.0;
  }
  return 0;
}

extern "C"
int fclow (void* user_data, int id, double* vec, int len)
{  

  PbData* data = (PbData*)user_data;
  int N=data->N; int T=data->T; int NW=data->NW; int ND=data->ND;
  
  if(id==0) {
    int it=0;
    // ------- 2 ------
    double** K = data->K; 
    int* v0 = data->v0;
    for(int t=1; t<=ND; t++) {
      for(int k=0; k<t; k++) {
	for(int j=0; j<N; j++) vec[it++] = -K[j][t-1] * (t-k) * v0[j]; //!1
      }

      for(int k=t; k<T; k++) {
	for(int j=0; j<N; j++) vec[it++] = 0.0;
      }
    } // end 't' iteration
#ifdef DEBUG
    assert(it==N*T*ND);
#endif
    // ------ 3 ------
    double* C = data->C; 
    //first N
    for(int j=0; j<N; j++) vec[it++] = C[j]*v0[j];
    //next N*(T-1)
    for(int i=0; i<N*(T-1); i++) vec[it++] = 0.0;

#ifdef DEBUG
    assert(it==N*T*ND+N*T);
#endif

    // ------ 5 ------
    double* demand = data->demand; int NWT = data->NWT;
    double* wind_scen = new double[NW]; // wind scenario for k=0
    data->loadWindScenario(0, 0, wind_scen, NW);

    double sum_wind = 0.0; for(int j=0; j<NW; j++) sum_wind += wind_scen[j];
    delete[] wind_scen;

    vec[it++] = demand[0]-NWT*sum_wind;

    // ------ 6 ------
    double* R = data->R; //reserve
    vec[it++] = R[0] + demand[0] - NWT*sum_wind;

    // ------ 8 ------
    for(int j=0; j<N; j++) vec[it++] = 0.0;

    // ------ 9 ------
    for(int j=0; j<N; j++) vec[it++] = 0.0;

    // ------ 10 ------
    for(int j=0; j<N; j++) vec[it++] = 0.0;
#ifdef DEBUG
    assert(it==N*T*ND+N*T +1+1+ 1*(N+N+N));
#endif

    // ------ 17 ------
    int* UT = data->UT; int* G = data->G;
    for(int j=0; j<N; j++) {
      //might happen k=G[j]=0, in which case the rhs is not 0.0
      vec[it++] = (G[j]==0 ? 0.0-UT[j]*v0[j] : 0.0);
      for(int k=G[j]+1; k<T-UT[j]+1; k++) vec[it++] = 0.0;
    }
    // ------ 18 ------
    for(int j=0; j<N; j++) {
      for(int k=T-UT[j]+2; k<=T; k++) vec[it++] = 0.0;
    }
    // ------ 20 ------
    int* DT = data->DT; int* L = data->L;
    for(int j=0; j<N; j++) {
      vec[it++] = (L[j]==0 ? -DT[j]+DT[j]*v0[j] : -DT[j]);
      for(int k=L[j]+1; k<T-DT[j]+1; k++) vec[it++] = -DT[j];
    }
    // ------ 21 ------
    for(int j=0; j<N; j++) {
      for(int k=T-DT[j]+2; k<=T; k++) vec[it++] = -1.0-T+k;

    }
#ifdef DEBUG
    int mz = data->get_mz(1);
    assert(it==mz);

    //for(int i=start; i<end; i++)
    //  printf("[%d]=%g\n", i, vec[i]);
    //assert(false);
#endif
  } else {
    int it=0;
    // ------ 5 & 6 ------ in the same time to reduce the calls to load... fcn
    double* demand = data->demand; double* R = data->R;
    int NWT = data->NWT;
    double* wind_scen = new double[NW];
    for(int k=0; k<T-1; k++) {
      data->loadWindScenario(id-1, k+1, wind_scen, NW);
      double sum=0.0; 
      for(int j=0; j<NW; j++) sum += wind_scen[j];

      sum *= NWT;
      vec[it] = demand[k+1] - sum;
      vec[it+T-1] = vec[it] + R[k+1];
      it++;
    } 
    it += (T-1);
    delete[] wind_scen;
#ifdef DEBUG
    assert(it==2*(T-1));
#endif

    // ------ 8 9 10  ------
    int n=3*N*(T-1);
    for(int i=0; i<n; i++) vec[it++] = 0.0;
#ifdef DEBUG
    assert(it==2*(T-1) + 3*N*(T-1));
#endif

    // ------ 11-----
    double* Pmax = data->Pmax;
    for(int k=0; k<T-1; k++) 
      for(int j=0; j<N; j++)
	vec[it++] = -Pmax[j];

    // ------ 12 -----
    for(int k=0; k<T-1; k++) 
      for(int j=0; j<N; j++)
	vec[it++] = -Pmax[j];
#ifdef DEBUG
    int mz = data->get_mz(2);
    assert(it==mz);
#endif
  }

  for(int i=0; i<len; i++) vec[i] /= SCALE_C;
  return 0; 
}

extern "C"
int ficlow(void* user_data, int id, double* vec, int len)
{ 
  for(int i=0; i<len; i++) vec[i] = 1.0;
  return 0; 
}

extern "C"
int fcupp (void* user_data, int id, double* vec, int len)
{ 
  for(int i=0; i<len; i++) vec[i] = 0.0;
  //for(int i=0; i<len; i++) vec[i] /= SCALE_C;
  return 0; 
}


extern "C"
int ficupp(void* user_data, int id, double* vec, int len)
{ 
  for(int i=0; i<len; i++) vec[i] = 0.0;
  return 0; 
}


extern "C"
int fxlow (void* user_data, int id, double* vec, int len)
{
  PbData* data = (PbData*)user_data;
  int N=data->N; int ND=data->ND; int T=data->T; int NW=data->NW; int S=data->S;

  if(id==0) {
    assert(len==3*N*T+(2*N+NW));
    for(int i=0; i<len; i++) vec[i]= 0.0; //no limit or 0 as lower limit
    
  } else {
    assert((T-1)*(2*N+NW)==len);
    for(int i=0; i<len; i++) vec[i]=0.0; //no limit or 0 as lower limit
  }
  return 0;
}


extern "C"
int fixlow(void* user_data, int id, double* vec, int len)
{
  PbData* data = (PbData*)user_data;
  int N=data->N; int T=data->T; int NW=data->NW; int S=data->S;
  int NT = N*T;
  if(id==0) {
    for(int i=0; i<NT; i++) vec[i]=1.0; //lower limit on v
    for(int i=NT; i<3*NT; i++) vec[i] = 1.0; //lower limit on cu and cd;
    for(int i=3*NT; i<3*NT+N; i++) vec[i] = 1.0; //lower limit on pc[j,1]

    for(int i=3*NT+N; i<3*NT+2*N; i++) vec[i] = 1.0; //lower limit on pcmax[j,1]
    for(int i=3*NT+2*N; i<3*NT+2*N+NW; i++) vec[i] = 1.0; //lower limit on pwind[j,1,1]
  } else {
    int NTm1=N*(T-1); assert(len==2*NTm1+NW*(T-1));
    for(int i=0; i<NTm1; i++) vec[i]=1.0; //low limit on pc
    //for(int i=NTm1; i<2*NTm1; i++) vec[i]=0.0; //no low limit cp //!2
    for(int i=NTm1; i<NTm1+NW*(T-1); i++) vec[i]=1.0; //low on pcwind
    for(int i=NTm1+NW*(T-1); i<2*NTm1+NW*(T-1); i++) vec[i]=1.0; //low on pcmax
    assert(len == 2*NTm1+NW*(T-1));
  }
  return 0;
}

extern "C"
int fxupp (void* user_data, int id, double* vec, int len)
{
  PbData* data = (PbData*)user_data;
  int N=data->N; int T=data->T; int NW=data->NW; int S=data->S;
  int NT = N*T;

  if(id==0) {
#ifdef DEBUG
    int nx = data->get_nx(1);
    assert(nx==len);
    assert(nx==3*NT+2*N+NW);
#endif
    for(int i=0; i<N*T; i++) vec[i]=1.0; //upper on v
    for(int i=N*T; i<3*NT+2*N; i++) vec[i]=0.0; // no limit on cu, cd, pc[j,1,1] and pcmax
    
    data->loadWindScenario(0, 0, &vec[3*NT+2*N], NW); 
  } else { // 2nd stage
#ifdef DEBUG
    int nx = data->get_nx(2);
    assert(nx==len);
#endif
 
    for(int i=0; i<N*(T-1); i++) vec[i]=0.0; // no upper limit for pc //!2
    
    for(int k=1; k<T; k++)
      data->loadWindScenario(id-1, k, &vec[N*(T-1)+(k-1)*NW], NW);

    // no upper limit for pcmax
    for(int i=N*(T-1)+NW*(T-1); i<2*N*(T-1)+NW*(T-1); i++) vec[i] = 0.0;
#ifdef DEBUG
    assert(nx==2*N*(T-1)+NW*(T-1));
#endif
  }
  return 0;
}

extern "C"
int fixupp(void* user_data, int id, double* vec, int len)
{
  PbData* data = (PbData*)user_data;
  int N=data->N; int T=data->T; int NW=data->NW; int S=data->S;
  int NT = N*T;

  if(id==0) {
    for(int i=0; i<N*T; i++) vec[i]=1.0; //upper on v
    for(int i=N*T; i<3*NT+2*N; i++) vec[i]=0.0; // no limit //!2
    for(int i=3*NT+2*N; i<3*NT+2*N+NW; i++) vec[i]=1.0; //pwind upper limit
  } else {
    for(int i=0; i<N*(T-1); i++) vec[i]=0.0; // no upper limit for pc//! 2
    for(int i=N*(T-1); i<N*(T-1)+NW*(T-1); i++) vec[i] = 1.0; // upper limit for pwind//! 2
    // no upper limit for pcmax //! 2
    for(int i=N*(T-1)+NW*(T-1); i<2*N*(T-1)+NW*(T-1); i++) vec[i] = 0.0; 
  }
  //for(int i=0; i<len; i++) vec[i]=1.0;
  return 0;
}

int PbData::get_nx(int stage)
{ return stage==1?3*N*T+(2*N+NW) : (T-1)*(2*N+NW); }

int PbData::get_my(int stage)
{ 
  int my;
  if(stage==1) {
    //!2 my = N; //1
    my=0;
    for(int j=0; j<N; j++) {
      if(G[j]>=1) my++; //16
      if(L[j]>=1) my++; //19
    }
  } else {
    assert(stage==2);
    //!2 my = N*(T-1); //1
    my = 0;
  }
  return my;
}
int PbData::get_mz(int stage)
{
  int mz;
  if(stage==1) {
    //    mz = N*T*(1+ND) + N; //2 3
    mz = N*T*ND  +  N*T; //2 + 3
    mz += (1+1+N+N+N); //5 6 8 9 10
    
    for(int j=0; j<N; j++) {
      assert(T-UT[j]-G[j]+1 >= 0);
      assert(UT[j]-1        >= 0);
      assert(T-DT[j]-L[j]+1 >= 0);
      assert(DT[j]-1        >= 0);
      
      mz += (T-G[j]); //17 + 18
      mz += (T-L[j]); //20 + 21
    }
  } else {
    assert(stage==2);
    mz = 2*(T-1)+5*N*(T-1);//5,6 & 8,9,10,11,12
    
  }
  return mz;
}
void PbData::compute1stStageSizes(int& nx, int& my, int&mz)
{
  nx = get_nx(1);
  my = get_my(1);
  mz = get_mz(1);
}

void PbData::compute2ndStageSizes(int& nx, int& my, int&mz)
{
  nx = get_nx(2);
  my = get_my(2);
  mz = get_mz(2);
}

void PbData::loadWindScenario(int scenario, int tmStep, //time 0..T-1
			      double* data, int noWindUnits)
{

  assert(noWindUnits==NW); 
  assert(tmStep>=0   && tmStep<T);
  assert(scenario>=0 && scenario<S);


  if(!scen_loaded) {
    wind_scen = new double**[S];
    for(int s=0; s<S; s++) {
      wind_scen[s] = new double*[T];
      for(int k=0; k<T; k++) {
	wind_scen[s][k] = new double[NW];
	for(int j=0; j<NW; j++) wind_scen[s][k][j] = 0.0;
      }
    }


    char file[64]; sprintf(file, "wind_scen%03d.dat", S);
    string strFile = g_strDataDir+file;
    fstream file_op(strFile.c_str(), ios::in);
    if (!file_op.is_open()) {
      int mype; MPI_Comm_rank(MPI_COMM_WORLD, &mype);
      printf("ERROR: processor %d couldn't open scenario file, message: %s\n",mype,strerror(errno));
      assert(0);
    }
    
    double dummy; int e=0;
    for(int n=0; n<NW; n++)
      for(int k=0; k<T; k++)
	for(int s=0; s<S; s++) {
	  /*
	    file_op >> dummy;
	    if(s==scenario && k==tmStep)
	    data[n] = dummy;
	  */
	  file_op >> wind_scen[s][k][n];
	}
    
    file_op.close();
    scen_loaded=1;
  } 


  if(tmStep==0 && scenario==0) {
    
    for(int n=0; n<NW; n++) {
      double max=0.0; double mean=0.0; double min=1.0e10;
      for(int s=0; s<S; s++) {
      
	if(wind_scen[s][0][n]>max) {	
	  max = wind_scen[s][0][n];
	} 

	if(wind_scen[s][0][n]<min) {
	  min = wind_scen[s][0][n];
	}

	mean += wind_scen[s][0][n];
      }
      //printf("unit %d max=%g\n", n, max);
      //data[n] = (mean/S);
      data[n] = min;
      //data[n] = max;
    }
  } else {

    for(int n=0; n<NW; n++)
      data[n] = wind_scen[scenario][tmStep][n];
  }

  //for(int n=0; n<NW; n++) cout << wind_scen[29][14][n]  << endl;
  //assert(false);
}

#define divN 10
void PbData::loadFromFile(const string& file)
{
  string strFile = g_strDataDir+file;
  fstream file_op(strFile.c_str(), ios::in);
  if (!file_op.is_open()) {
    int mype; MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    printf("ERROR: processor %d couldn't open config file at %s , message: %s\n",mype,strFile.c_str(),strerror(errno));
    assert(0);
  }
  
  string strLine;istringstream sstrLine;
  getline(file_op, strLine); sstrLine.str(strLine);
  sstrLine >> N; //Number of units
  getline(file_op, strLine); sstrLine.str(strLine);
  sstrLine >> T; //Number of periods
  getline(file_op, strLine); sstrLine.str(strLine);
  sstrLine >> ND;//Number of intervals of the stairwise startup cost function
#ifdef TIMING
  if (OVERunits > 0) N = OVERunits;
#endif
  // Pmax - Maximum power output
  Pmax = new double[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, Pmax);

  // Pmin - Minimum power output
  Pmin = new double[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, Pmin);
  
  // SU - Startup ramp limit
  SU = new double[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, SU);
 
  // SD - Shutdown ramp limt
  SD = new double[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, SD);

  // RU - Ramp-up limit
  RU = new double[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, RU);

  //RD - Ramp-down limit
  RD = new double[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, RD);

  // UT - Minimum up time
  UT = new int[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, UT);

  // DT - Minimum down time
  DT = new int[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, DT);

  // v0 - Initial state
  v0 = new int[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, v0);

  // a - Coefficients of the quadratic production cost function(a)
  a = new double[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, a);

  // b - Coefficients of the quadratic production cost function(b)
  b = new double[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, b);

  // c - Coefficients of the quadratic production cost function(c)
  c = new double[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, c);

  // hc - Coefficients of the startup cost function
  hc = new double[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, hc);

  // cc - Coefficients of the startup cost function
  cc = new double[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, cc);

  // tcold - Coefficients of the startup cost function
  tcold = new double[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, tcold);

  // C - Shutdown cost
  C = new double[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, C);

  // demand - expected load demand (next T entries, may be on more than one line) 
  demand = new double[T];
  for(int i=0; i<T; i++) { file_op >> demand[i]; demand[i] *= (1.0*N/10); }//!/10
  getline(file_op, strLine); //skip comments and or empty lines

  // Rp - reserve requierement (next T entries, may be on more than one line)
  Rp = new double[T];
  for(int i=0; i<T; i++) file_op >> Rp[i];
  getline(file_op, strLine); //skip comments and or empty lines

  // R - reserve requierement (next T entries, may be on more than one line)
  R = new double[T];
  for(int i=0; i<T; i++) R[i] = Rp[i]*N/10;

  // initstate - 
  initstate = new int[N];
  getNextNonEmptyLine(file_op, strLine); 
  paramFromLine(strLine, initstate);

  // S - number of scenarios
  getNextNonEmptyLine(file_op, strLine); sstrLine.str(strLine); sstrLine >> S;
#ifdef TIMING
  if (OVERscen > 0 ) S = OVERscen;
#endif
  // NW - number of wind farms(units)
  getNextNonEmptyLine(file_op, strLine); sstrLine.str(strLine); sstrLine >> NW;

  // NWT - number of wind turbines in a wind farm
  getNextNonEmptyLine(file_op, strLine); sstrLine.str(strLine); sstrLine >> NWT; 
  NWT *= (N/10);

  //S0 -  number of offline periods prior to first
  S0 = new int[N];
  for(int i=0; i<N; i++) 
    if(initstate[i]<0.0) S0[i] = -initstate[i];
    else                 S0[i] = 0;

  //U0 -  number of online periods prior to first
  U0 = new int[N];
  for(int i=0; i<N; i++) 
    if(initstate[i]>0.0) U0[i] = initstate[i];
    else                 U0[i] = 0;
  
  //G - number of initial periods during which unit j must be online
  G = new int[N];
  for(int i=0; i<N; i++) 
    G[i] = min(T, (UT[i]-U0[i])*v0[i]);


  //L - number of initial periods during which unit j must be online
  L = new int[N];
  for(int i=0; i<N; i++) 
    L[i] = min(T,(DT[i]-S0[i])*(1-v0[i]));
  //for(int i=0; i<N; i++) cout << L[i] << "  ";
  //for(int i=0; i<T; i++) cout << R[i] << " "; cout << endl;

  K = new double*[N];
  for(int j=0; j<N; j++) {
    K[j] = new double[ND];
    for(int t=0; t<ND; t++)
      if(t<=DT[j]+tcold[j])
	K[j][t] = hc[j];
      else K[j][t] = cc[j];
  }
  file_op.close();

  //for(int j=0; j<N; j++) cout << v0[j] << endl;
  //for(int k=0; k<T; k++)  cout << [k]  << endl;
    //for(int k=0; k<ND; k++) 
    //  printf("[%d , %d ]=%g\n", j+1,k+1, K[j][k]);

}

int PbData::getNextNonEmptyLine(fstream& file, string& line)
{
  line=""; 
  while(line.size()==0) {
    getline(file, line);
  }
  return 1;
}

void PbData::paramFromLine(const string& strLine, double* param)
{
  istringstream sstrLine(strLine);
  for(int i=0; i<divN; i++)   sstrLine >> param[i];
  for(int j=1; j<N/divN; j++) 
    for(int i=0; i<divN; i++)
      param[j*divN+i] = param[i];
}
void PbData::paramFromLine(const string& strLine, int* param)
{
  istringstream sstrLine(strLine);
  for(int i=0; i<divN; i++)   sstrLine >> param[i];
  for(int j=1; j<N/divN; j++) 
    for(int i=0; i<divN; i++)
      param[j*divN+i] = param[i];
}

PbData::PbData(const string& strFile)
  : scen_loaded(0)
{
  loadFromFile(strFile);
}

PbData::~PbData()
{
  delete[] Pmax;
  delete[] Pmin;
  delete[] SU;
  delete[] SD;
  delete[] RU;
  delete[] RD;
  delete[] UT;
  delete[] DT;
  delete[] v0;
  delete[] a;
  delete[] b;
  delete[] c;
  delete[] hc;
  delete[] cc;
  delete[] tcold;
  delete[] C;
  delete[] demand;
  delete[] Rp;
  delete[] R;
  delete[] initstate;
  delete[] S0;
  delete[] U0;
  delete[] G;
  delete[] L;

  for(int j=0; j<N; j++) 
    delete[] K[j];
  delete[] K;
  if(scen_loaded) {
    for(int s=0; s<S; s++) {
      for(int k=0; k<T; k++) {
	delete[] wind_scen[s][k];
      }
      delete[] wind_scen[s];
    }
    delete[] wind_scen;
  }
}


void qpgenSolve(PbData& data, int printx)
{
  int NX, MY, MZ;
  int nnzA, nnzQ=0, nnzC;

  int S = data.S; int T = data.T; int N = data.N;

  NX = data.get_nx(1)+ S*data.get_nx(2);
  MY = data.get_my(1)+ S*data.get_my(2);
  MZ = data.get_mz(1)+ S*data.get_mz(2);
  printf("NX=%d MY=%d MZ=%d\n", NX, MY, MZ);

  int tmp;
  fnnzA(&data, 0, &nnzA); fnnzC(&data, 0, &nnzC);  
  fnnzB(&data, 0, &tmp); nnzA += tmp;
  fnnzD(&data, 0, &tmp); nnzC += tmp;
  for(int i=1; i<=S; i++) {

    fnnzA(&data, i, &tmp); nnzA+= tmp;
    fnnzB(&data, i, &tmp); nnzA+= tmp;

    fnnzC(&data, i, &tmp); nnzC+= tmp;
    fnnzD(&data, i, &tmp); nnzC+= tmp;
  }
  printf("nnzQ=%d nnzA=%d nnzC=%d\n", nnzQ, nnzA, nnzC);

  /////////////////////////////////////////
  //Q
  /////////////////////////////////////////
  int* krowQ=new int[NX+1]; int*jcolQ=new int[nnzQ]; double* dQ = new double[nnzQ];
  for(int i=0; i<=NX; i++) krowQ[i] = 0;

  /////////////////////////////////////////
  // A and b
  /////////////////////////////////////////
  int* krowA=new int[MY+1]; int*jcolA=new int[nnzA]; double* dA = new double[nnzA];
  double* b = new double[MY];
  // populate A

  int rowind=0, colind=0; 
  //root
  {
    int my0 = data.get_my(1);
    int NNZA; fnnzA(&data, 0, &NNZA);
    fA(&data, 0, krowA, jcolA, dA);
    rowind = my0; colind = NNZA;
    
    fb(&data, 0, &b[0], my0);
  }
  for(int s=1; s<=S; s++) {
    int NNZA; fnnzA(&data, s, &NNZA);
    int NNZB; fnnzB(&data, s, &NNZB);
    int my  = data.get_my(2); int my0= data.get_my(1);
    int nx0 = data.get_nx(1); int nx = data.get_nx(2);
    
    int* krowAs=new int[my+1]; int*jcolAs=new int[NNZA]; double* dAs = new double[NNZA];
    int* krowBs=new int[my+1]; int*jcolBs=new int[NNZB]; double* dBs = new double[NNZB];

    fA(&data, s, krowAs, jcolAs, dAs);
    fB(&data, s, krowBs, jcolBs, dBs);
    
    for(int i=0; i<my; i++) {
      krowA[rowind] = colind;
      
      for(int k=krowAs[i]; k<krowAs[i+1]; k++) {
	jcolA[colind] = jcolAs[k]; dA[colind] = dAs[k]; colind++;
      }
      for(int k=krowBs[i]; k<krowBs[i+1]; k++) {
	jcolA[colind] = jcolBs[k]+nx0+(s-1)*nx;
	dA[colind] = dBs[k]; colind++;
      }

      krowA[rowind+1]=colind; rowind++;
    }
#ifdef DEBUG
    assert(rowind==my0+s*my);
#endif
    
    delete[] krowAs; delete[] jcolAs; delete[] dAs;
    delete[] krowBs; delete[] jcolBs; delete[] dBs;

    //b
    fb(&data, s, &b[my0+(s-1)*my], my);
  }
  assert(rowind==MY); 
  assert(colind==nnzA);
  //~ done with A and b

  /////////////////////////////////////////
  // C and (i)clow, (i)cupp
  /////////////////////////////////////////
  int* krowC=new int[MZ+1]; int*jcolC=new int[nnzC]; double* dC = new double[nnzC];
  char* iclow=new char[MZ]; double* clow=new double[MZ];
  char* icupp=new char[MZ]; double* cupp=new double[MZ];

  rowind=0, colind=0; 
  //root
  {
    int mz0=data.get_mz(1);
    int NNZC; fnnzC(&data, 0, &NNZC);
    fC(&data, 0, krowC, jcolC, dC);
    rowind=mz0; colind=NNZC;
    double* iclowd=new double[mz0]; double* icuppd=new double[mz0];

    ficlow(&data, 0, iclowd, mz0); fclow(&data, 0, clow, mz0);
    ficupp(&data, 0, icuppd, mz0); fcupp(&data, 0, cupp, mz0);
    //convert from double to char
    for(int i=0; i<mz0; i++) {
      if(iclowd[i]==0.0) iclow[i]=0; else iclow[i]=1;
      if(icuppd[i]==0.0) icupp[i]=0; else icupp[i]=1;
    }
    delete[] iclowd; delete[] icuppd;
  }
  for(int s=1; s<=S; s++) {
    int NNZC; fnnzC(&data, s, &NNZC);
    int NNZD; fnnzD(&data, s, &NNZD);
    int mz  = data.get_mz(2); int mz0= data.get_mz(1);
    int nx0 = data.get_nx(1); int nx = data.get_nx(2);
    
    int* krowCs=new int[mz+1]; int*jcolCs=new int[NNZC]; double* dCs = new double[NNZC];
    int* krowDs=new int[mz+1]; int*jcolDs=new int[NNZD]; double* dDs = new double[NNZD];
    

    fC(&data, s, krowCs, jcolCs, dCs);
    fD(&data, s, krowDs, jcolDs, dDs);
    
    for(int i=0; i<mz; i++) {
      krowC[rowind] = colind;
      for(int k=krowCs[i]; k<krowCs[i+1]; k++) {
	jcolC[colind] = jcolCs[k];
	dC   [colind] = dCs[k];
	colind++;
      }
      for(int k=krowDs[i]; k<krowDs[i+1]; k++) {
	jcolC[colind] = jcolDs[k] + nx0+(s-1)*nx;;
	dC   [colind] = dDs[k];
	colind++;
      }
      krowC[rowind+1] = colind; rowind++;
    }
    delete[] krowCs; delete[] jcolCs; delete[] dCs;
    delete[] krowDs; delete[] jcolDs; delete[] dDs;

#ifdef DEBUG
    assert(rowind==mz0+s*mz);
#endif

    double* iclowd=new double[mz]; double* icuppd=new double[mz];
    ficlow(&data, s, iclowd, mz); fclow(&data, s, &clow[mz0+(s-1)*mz], mz);
    ficupp(&data, s, icuppd, mz); fcupp(&data, s, &cupp[mz0+(s-1)*mz], mz);

    for(int i=0; i<mz; i++) {
      int gi = i+mz0+(s-1)*mz;
      if(iclowd[i]==0.0) iclow[gi]=0; else iclow[gi]=1;
      if(icuppd[i]==0.0) icupp[gi]=0; else icupp[gi]=1;
    }
    delete[] iclowd; delete[] icuppd;
  }
  assert(rowind==MZ); 
  assert(colind==nnzC);
  //~ done with C and and (i)clow, (i)cupp

  ///////////////////////////////////////////
  // (i)xlow and (i)xupp and c
  //////////////////////////////////////////
  char* ixlow=new char[NX];double* xlow=new double[NX];
  char* ixupp=new char[NX];double* xupp=new double[NX];
  double* c = new double[NX];
  //root
  {
    int nx0 = data.get_nx(1);

    fxlow(&data, 0, xlow, nx0); fxupp(&data, 0, xupp, nx0);

    double* ixlowd=new double[nx0]; double* ixuppd=new double[nx0];
    fixlow(&data, 0, ixlowd, nx0); fixupp(&data, 0, ixuppd, nx0);
    for(int i=0; i<nx0; i++) {
      if(ixlowd[i]==0.0) ixlow[i]=0; else ixlow[i]=1;
      if(ixuppd[i]==0.0) ixupp[i]=0; else ixupp[i]=1;
    }
    delete[] ixlowd; delete[] ixuppd;

    fc(&data, 0, c, nx0);
  }
  for(int s=1; s<=S; s++) {
    int nx0 = data.get_nx(1); int nx = data.get_nx(2);

    fxlow(&data, s, &xlow[nx0+(s-1)*nx], nx); fxupp(&data, s, &xupp[nx0+(s-1)*nx], nx);

    double* ixlowd=new double[nx]; double* ixuppd=new double[nx];
    fixlow(&data, s, ixlowd, nx); fixupp(&data, s, ixuppd, nx);

    for(int i=0; i<nx; i++) {
      int gi = i+nx0+(s-1)*nx;
      if(ixlowd[i]==0.0) ixlow[gi]=0; else ixlow[gi]=1;
      if(ixuppd[i]==0.0) ixupp[gi]=0; else ixupp[gi]=1;
    }
    delete[] ixlowd; delete[] ixuppd;

    fc(&data, s, &c[nx0+(s-1)*nx], nx);
  }

  /******************************************
   * scaling related stuff
   *****************************************/
  /*{
    double maxC=0.0, maxcb=0.0, maxxbl=0.0, maxxbu=0.0, maxobj=0.0;
    for(int i=0; i<nnzC; i++) if(maxC<fabs(dC[i])) maxC = fabs(dC[i]);
    for(int i=0; i<MZ; i++) if(maxcb<fabs(clow[i])) maxcb = fabs(clow[i]);
    for(int i=0; i<NX; i++) {
      if(maxxbl<fabs(xlow[i])) maxxbl = fabs(xlow[i]);
      if(maxxbu<fabs(xupp[i])) maxxbu = fabs(xupp[i]);
      if(maxobj<fabs(   c[i])) maxobj = fabs(   c[i]);
    }
    printf("maxC=%g  maxcb=%g  maxxbl=%g  maxxbu=%g  maxobj=%g\n",
	   maxC, maxcb, maxxbl, maxxbu, maxobj);
    
    for(int i=0; i<nnzC; i++) dC[i] /= (maxC/2);
    for(int i=0; i<MZ; i++) clow[i] /= (maxC/2);
    for(int i=0; i<NX; i++)    c[i] /= (maxobj/2);
    }*/
  {
    double maxC=0.0, maxcb=0.0, maxxbl=0.0, maxxbu=0.0, maxobj=0.0;
    for(int i=0; i<nnzC; i++) if(maxC<fabs(dC[i])) maxC = fabs(dC[i]);
    for(int i=0; i<MZ; i++) if(maxcb<fabs(clow[i])) maxcb = fabs(clow[i]);
    for(int i=0; i<NX; i++) {
      if(maxxbl<fabs(xlow[i])) maxxbl = fabs(xlow[i]);
      if(maxxbu<fabs(xupp[i])) maxxbu = fabs(xupp[i]);
      if(maxobj<fabs(   c[i])) maxobj = fabs(   c[i]);
    }
    printf("maxC=%g  maxcb=%g  maxxbl=%g  maxxbu=%g  maxobj=%g\n",
	   maxC, maxcb, maxxbl, maxxbu, maxobj);
  }
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
  
  //rusage before_solve;
  //getrusage( RUSAGE_SELF, &before_solve );

  int result = s->solve(prob,vars, resid);

  //rusage  after_solve;
  //getrusage( RUSAGE_SELF, &after_solve );  

  //double solve_time =
  //  (after_solve.ru_utime.tv_sec - before_solve.ru_utime.tv_sec)
  //  + (after_solve.ru_utime.tv_usec - before_solve.ru_utime.tv_usec)
  //  / 1000000.0;
  
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

  //cout  << "QP solved in " << solve_time << " seconds.\n";

  delete s;
  delete vars;  
  delete resid;
  delete prob;
  delete qp;
}
