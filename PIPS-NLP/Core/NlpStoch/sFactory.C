/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "mpi.h"
#include "sFactory.h"

#include "sData.h"
#include "stochasticInput.hpp"

#include "sTreeImpl.h"

#include "sTreeCallbacks.h"
#include "StochInputTree.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"


#include "sVars.h"
#include "sResiduals.h"

#include "sLinsysRoot.h"
#include "sLinsysLeaf.h"

#ifdef WITH_MA27
#include "Ma27Solver.h"
#endif
#ifdef WITH_MA57
#include "Ma57Solver.h"
#endif
#ifdef WITH_PARDISO
#include "PardisoSolver.h"
#endif

#include "DeSymIndefSolver.h"

#ifdef WITH_UMFPACK
#include "UmfPackSolver.h"
#endif 


#include <stdio.h>
#include <stdlib.h>

#include "NlpInfo.h"
#include "NlpGen.h"

#include "sInfo.h"
#include "sData.h"

using namespace std;

extern int gSymLinearSolver;


sFactory::sFactory( stochasticInput& in, MPI_Comm comm)
  : NlpGen(0,0,0), data(NULL), m_tmTotal(0.0)
{
  tree = new sTreeImpl(in, comm);

  //decide how the CPUs are assigned 
  tree->assignProcesses(comm);
  tree->loadLocalSizes();

  tree->computeGlobalSizes();
  tree->GetGlobalSizes(nx, my, mz);

  datarootname = in.datarootname;
}

 
sFactory::sFactory( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : NlpGen( nx_, my_, mz_ ),
    nnzQ(nnzQ_), nnzA(nnzA_), nnzC(nnzC_),
    tree(NULL), data(NULL), resid(NULL), linsys(NULL),
    m_tmTotal(0.0)
{ };

sFactory::sFactory()
  : NlpGen( 0,0,0 ), data(NULL), m_tmTotal(0.0)
{ };

sFactory::~sFactory()
{
  if(tree) delete tree;
}

sLinsysLeaf* 
sFactory::newLinsysLeaf()
{
  assert(false && "not supported");
  return NULL;
}
sLinsysLeaf* 
sFactory::newLinsysLeaf(sData* prob,
			OoqpVector* dd, OoqpVector* dq,
			OoqpVector* nomegaInv, OoqpVector* rhs, OoqpVector* additiveDiag_)
{
  sLinsysLeaf *resultSLin=NULL;
#ifdef TIMING
  cout<< "newLinsysLeaf:" <<gSymLinearSolver << endl;
#endif
  if(0==gSymLinearSolver){
#ifdef WITH_MA27
	Ma27Solver* sMA27=NULL; 
	resultSLin = new sLinsysLeaf(this, prob, dd, dq, nomegaInv, rhs, additiveDiag_, sMA27);
#endif    
  }else if(1==gSymLinearSolver){
#ifdef WITH_MA57
	Ma57Solver* sMA57=NULL; 
	resultSLin = new sLinsysLeaf(this, prob, dd, dq, nomegaInv, rhs, additiveDiag_, sMA57);
#endif   
  }else if(2==gSymLinearSolver){
#ifdef WITH_PARDISO
	PardisoSolver* sPardiso=NULL; 
	resultSLin = new sLinsysLeaf(this, prob, dd, dq, nomegaInv, rhs, additiveDiag_, sPardiso);
#endif
  }else if(3==gSymLinearSolver){
#ifdef WITH_UMFPACK
	UmfPackSolver* sUmfPack=NULL; 
    resultSLin = new sLinsysLeaf(this, prob, dd, dq, nomegaInv, rhs, additiveDiag_, sUmfPack);
#endif
  }else if(4==gSymLinearSolver){ //sMA57 as dummy
#ifdef WITH_MA57
//	Ma57Solver* sMA57=NULL; 
//	resultSLin = new sLinsysLeaf(this, prob, dd, dq, nomegaInv, rhs, additiveDiag_, sMA57);
#endif 
  }
  if(NULL==resultSLin) assert(0 && "PIPS was not compiled with the linear solver required");

  return resultSLin;
}


#ifdef TIMING
#define TIM t = MPI_Wtime();
#define REP(s) t = MPI_Wtime()-t;if (p) printf("makeData: %s took %f sec\n",s,t);
#else
#define TIM
#define REP(s)
#endif


void dumpaug(int nx, SparseGenMatrix &A, SparseGenMatrix &C) {

    long long my, mz, nx_1, nx_2;
    A.getSize(my,nx_1);
    C.getSize(mz,nx_2);
    assert(nx_1 == nx_2);
	
    int nnzA = A.numberOfNonZeros();
    int nnzC = C.numberOfNonZeros();
    cout << "augdump  nx=" << nx << endl;
    cout << "A: " << my << "x" << nx_1 << "   nnz=" << nnzA << endl
	 << "C: " << mz << "x" << nx_1 << "   nnz=" << nnzC << endl;
	
	vector<double> eltsA(nnzA), eltsC(nnzC), elts(nnzA+nnzC);
	vector<int> colptrA(nx_1+1),colptrC(nx_1+1), colptr(nx_1+1), rowidxA(nnzA), rowidxC(nnzC), rowidx(nnzA+nnzC);
	A.getStorageRef().transpose(&colptrA[0],&rowidxA[0],&eltsA[0]);
	C.getStorageRef().transpose(&colptrC[0],&rowidxC[0],&eltsC[0]);
	
	int nnz = 0;
	for (int col = 0; col < nx_1; col++) {
		colptr[col] = nnz;
		for (int r = colptrA[col]; r < colptrA[col+1]; r++) {
			int row = rowidxA[r]+nx+1; // +1 for fortran
			rowidx[nnz] = row;
			elts[nnz++] = eltsA[r];
		}
		for (int r = colptrC[col]; r < colptrC[col+1]; r++) {
			int row = rowidxC[r]+nx+my+1;
			rowidx[nnz] = row;
			elts[nnz++] = eltsC[r];
		}
	}
	colptr[nx_1] = nnz;
	assert(nnz==nnzA+nnzC);

	ofstream fd("augdump.dat");
	fd << scientific;
	fd.precision(16);
	fd << (nx + my + mz) << endl;
	fd << nx_1 << endl;
	fd << nnzA+nnzC << endl;
	int i;
	for (i = 0; i <= nx_1; i++)
	fd << colptr[i] << " ";
	fd << endl;
	for (i = 0; i < nnz; i++)
	fd << rowidx[i] << " ";
	fd << endl;
	for (i = 0; i < nnz; i++)
	fd << elts[i] << " ";
	fd << endl;
	printf("finished dumping aug\n");


}

Data * sFactory::makeData()
{
  double t,t2=MPI_Wtime();
  int mype; MPI_Comm_rank(tree->commWrkrs,&mype);
  bool p = (mype == 0);

  // use the tree to get the data from user and to create OOQP objects

  TIM;
  StochGenMatrixHandle     Jeq( tree->createA() );
  REP("A");
  TIM;
  StochVectorHandle        b( tree->createb() );
  REP("b");
  TIM;
  StochGenMatrixHandle     Jineq( tree->createC() );  
  REP("C");
  TIM;
  StochVectorHandle     clow( tree->createclow()  );
  REP("clow");
  TIM;
  StochVectorHandle    iclow( tree->createiclow() );
  REP("iclow");
  TIM;
  StochVectorHandle     cupp( tree->createcupp()  );
  REP("cupp");
  TIM;
  StochVectorHandle    icupp( tree->createicupp() );
  REP("icupp");


  TIM;
  StochSymMatrixHandle     Q( tree->createQ() );
  REP("Q");
  TIM;
  StochVectorHandle        grad( tree->createc() );
  REP("c");
  TIM;
  StochVectorHandle     xlow( tree->createxlow()  );
  REP("xlow");
  TIM;
  StochVectorHandle    ixlow( tree->createixlow() );
  REP("ixlow");
  TIM;
  StochVectorHandle     xupp( tree->createxupp()  );
  REP("xupp");
  TIM;
  StochVectorHandle    ixupp( tree->createixupp() );
  REP("ixupp");


  TIM;
  StochVectorHandle	  CeqBody( tree->createCeqBody()	);
  REP("CeqBody");
  TIM;
  StochVectorHandle	 CIneqBody( tree->createCineqBody() );
  REP("CIneqBody");
  TIM;
  StochVectorHandle	  trialBarrGrad_x( tree->newPrimalVector()	);
  REP("trialBarrGrad_x");
  TIM;
  StochVectorHandle	 trialBarrGrad_s( tree->newDualZVector() );
  REP("trialBarrGrad_s");
  TIM;
  StochVectorHandle	  trialCeqBody( tree->newDualYVector()	);
  REP("trialCeqBody");
  TIM;
  StochVectorHandle	 trialCIneqBody( tree->newDualZVector() );
  REP("trialCIneqBody");
  
  TIM;
  StochVectorHandle dampind_xL_v  ( tree->newPrimalVector()  ); 
  REP("dampind_xL_v");
  TIM;
  StochVectorHandle dampind_xU_w  ( tree->newPrimalVector()  );  
  REP("dampind_xU_w");
  TIM;
  StochVectorHandle dampind_sL_t  ( tree->newDualZVector()   );
  REP("dampind_sL_t");
  TIM;
  StochVectorHandle dampind_sU_u  ( tree->newDualZVector()   );
  REP("dampind_sU_u");
  

  

#ifdef TIMING
  MPI_Barrier(tree->commWrkrs);
  t2 = MPI_Wtime() - t2;
  if (mype == 0) {
    cout << "IO second part took " << t2 << " sec\n";
  }
#endif

  
  data = new sData( 0, tree, 
		    grad, Q, 
		    xlow, ixlow, ixlow->numberOfNonzeros(),
		    xupp, ixupp, ixupp->numberOfNonzeros(),
		    Jeq, b,
		    Jineq, clow, iclow, iclow->numberOfNonzeros(),
		    cupp, icupp, icupp->numberOfNonzeros(),
		    CeqBody, CIneqBody,
		    trialBarrGrad_x, trialBarrGrad_s,
	        trialCeqBody, trialCIneqBody,
	        dampind_xL_v, dampind_xU_w,
	        dampind_sL_t, dampind_sU_u);

  data->datarootname = datarootname;
  	
  return data;
}



Data * sFactory::makeDataMulti()
{
  double t,t2=MPI_Wtime();
  int mype; MPI_Comm_rank(tree->commWrkrs,&mype);
  bool p = (mype == 0);

  // use the tree to get the data from user and to create OOQP objects

  TIM;
  StochGenMatrixHandle     Jeq( tree->createA() );
  REP("A");
  TIM;
  StochVectorHandle        b( tree->createb() );
  REP("b");
  TIM;
  StochGenMatrixHandle     Jineq( tree->createC() );  
  REP("C");
  TIM;
  StochVectorHandle     clow( tree->createclow()  );
  REP("clow");
  TIM;
  StochVectorHandle    iclow( tree->createiclow() );
  REP("iclow");
  TIM;
  StochVectorHandle     cupp( tree->createcupp()  );
  REP("cupp");
  TIM;
  StochVectorHandle    icupp( tree->createicupp() );
  REP("icupp");


  TIM;
  StochSymMatrixHandle     Q( tree->createQ() );
  REP("Q");
  TIM;
  StochVectorHandle        grad( tree->createc() );
  REP("c");
  TIM;
  StochVectorHandle     xlow( tree->createxlow()  );
  REP("xlow");
  TIM;
  StochVectorHandle    ixlow( tree->createixlow() );
  REP("ixlow");
  TIM;
  StochVectorHandle     xupp( tree->createxupp()  );
  REP("xupp");
  TIM;
  StochVectorHandle    ixupp( tree->createixupp() );
  REP("ixupp");


  TIM;
  StochVectorHandle	  CeqBody( tree->createCeqBody()	);
  REP("CeqBody");
  TIM;
  StochVectorHandle	 CIneqBody( tree->createCineqBody() );
  REP("CIneqBody");
  TIM;
  StochVectorHandle	  trialBarrGrad_x( tree->newPrimalVector()	);
  REP("trialBarrGrad_x");
  TIM;
  StochVectorHandle	 trialBarrGrad_s( tree->newDualZVector() );
  REP("trialBarrGrad_s");
  TIM;
  StochVectorHandle	  trialCeqBody( tree->newDualYVector()	);
  REP("trialCeqBody");
  TIM;
  StochVectorHandle	 trialCIneqBody( tree->newDualZVector() );
  REP("trialCIneqBody");
  
  TIM;
  StochVectorHandle dampind_xL_v  ( tree->newPrimalVector()  ); 
  REP("dampind_xL_v");
  TIM;
  StochVectorHandle dampind_xU_w  ( tree->newPrimalVector()  );  
  REP("dampind_xU_w");
  TIM;
  StochVectorHandle dampind_sL_t  ( tree->newDualZVector()   );
  REP("dampind_sL_t");
  TIM;
  StochVectorHandle dampind_sU_u  ( tree->newDualZVector()   );
  REP("dampind_sU_u");
  

  

#ifdef TIMING
  MPI_Barrier(tree->commWrkrs);
  t2 = MPI_Wtime() - t2;
  if (mype == 0) {
    cout << "IO second part took " << t2 << " sec\n";
  }
#endif

  
  data = new sData( 1, tree, 
		    grad, Q, 
		    xlow, ixlow, ixlow->numberOfNonzeros(),
		    xupp, ixupp, ixupp->numberOfNonzeros(),
		    Jeq, b,
		    Jineq, clow, iclow, iclow->numberOfNonzeros(),
		    cupp, icupp, icupp->numberOfNonzeros(),
		    CeqBody, CIneqBody,
		    trialBarrGrad_x, trialBarrGrad_s,
	        trialCeqBody, trialCIneqBody,
	        dampind_xL_v, dampind_xU_w,
	        dampind_sL_t, dampind_sU_u);

  data->datarootname = datarootname;
  	
  return data;
}


Data * sFactory::makeData(NlpInfo *updateNlp)
{
  data = dynamic_cast<sData*> (makeData());

  data->SetInputNlpPara(updateNlp);
 
  data->inputNlp = updateNlp; 

  return data;
}


Variables* sFactory::makeVariables( Data * prob_in )
{
  sData* prob = dynamic_cast<sData*>(prob_in);

  OoqpVectorHandle x      = OoqpVectorHandle( tree->newPrimalVector() );
  OoqpVectorHandle s      = OoqpVectorHandle( tree->newDualZVector() );
  OoqpVectorHandle y      = OoqpVectorHandle( tree->newDualYVector() );
  OoqpVectorHandle z      = OoqpVectorHandle( tree->newDualZVector() );

  OoqpVectorHandle v      ; 
  OoqpVectorHandle gamma  ;
  OoqpVectorHandle w      ;
  OoqpVectorHandle phi    ;
  OoqpVectorHandle t      ;
  OoqpVectorHandle lambda ;
  OoqpVectorHandle u      ;
  OoqpVectorHandle pi     ;



  if ( prob->nxlow > 0 ) { 
    v     = OoqpVectorHandle( tree->newPrimalVector() ); 
    gamma  = OoqpVectorHandle( tree->newPrimalVector() );
  }else{
	v 	= OoqpVectorHandle( tree->newPrimalVectorEmpty() ); 
	gamma  = OoqpVectorHandle( tree->newPrimalVectorEmpty() );
  }
  if ( prob->nxupp > 0 ) { 
    w     = OoqpVectorHandle( tree->newPrimalVector() ); 
    phi  = OoqpVectorHandle( tree->newPrimalVector() );
  }else{
	w 	= OoqpVectorHandle( tree->newPrimalVectorEmpty() ); 
	phi  = OoqpVectorHandle( tree->newPrimalVectorEmpty() );
  }
  if ( prob->mclow > 0 ) { 
    t     = OoqpVectorHandle( tree->newDualZVector() ); 
    lambda  = OoqpVectorHandle( tree->newDualZVector() );
  }else{
	t	  = OoqpVectorHandle( tree->newDualZVectorEmpty() ); 
	lambda  = OoqpVectorHandle( tree->newDualZVectorEmpty() );
  }  
  if ( prob->mcupp > 0 ) { 
    u     = OoqpVectorHandle( tree->newDualZVector() ); 
    pi  = OoqpVectorHandle( tree->newDualZVector() );
  }else{
	u 	= OoqpVectorHandle( tree->newDualZVectorEmpty() ); 
	pi  = OoqpVectorHandle( tree->newDualZVectorEmpty() );
  }
  
  // OoqpVector * s      = tree->newDualZVector();
  // OoqpVector * y      = tree->newDualYVector();
  // OoqpVector * z      = tree->newDualZVector();
  // OoqpVector * v      = tree->newPrimalVector(); 
  // OoqpVector * gamma  = tree->newPrimalVector();
  // OoqpVector * w      = tree->newPrimalVector(); 
  // OoqpVector * phi    = tree->newPrimalVector();
  // OoqpVector * t      = tree->newDualZVector();
  // OoqpVector * lambda = tree->newDualZVector();
  // OoqpVector * u      = tree->newDualZVector(); 
  // OoqpVector * pi     = tree->newDualZVector();

  sVars* vars = new sVars( tree, x, s, y, z,
			   v, gamma, w, phi,
			   t, lambda, u, pi, 
			   prob->ixlow, prob->ixlow->numberOfNonzeros(),
			   prob->ixupp, prob->ixupp->numberOfNonzeros(),
			   prob->iclow, prob->iclow->numberOfNonzeros(),
			   prob->icupp, prob->icupp->numberOfNonzeros());
  registeredVars.push_back(vars);
  return vars;
}


Residuals* sFactory::makeResiduals( Data * prob_in )
{
  sData* prob = dynamic_cast<sData*>(prob_in);
  resid =  new sResiduals(tree, 
			  prob->ixlow, prob->ixupp,
			  prob->iclow, prob->icupp);
  return resid; 
}



LinearSystem* sFactory::makeLinsys( Data * prob_in )
{  
  linsys = NULL;
  linsys = newLinsysRoot();

  assert(linsys);
  return linsys; 
}

void sFactory::joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			  OoqpVector& rhs2_in, OoqpVector& rhs3_in )
{
  assert(0 && "not implemented here");
}

void
sFactory::separateVars( OoqpVector& x_in, OoqpVector& y_in,
			OoqpVector& z_in, OoqpVector& vars_in )
{
  assert(0 && "not implemented here");
}

void sFactory::joinRHSXSYZ( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
		  OoqpVector& rhs2_in, OoqpVector& rhs3_in, OoqpVector& rhs4_in )
{
  assert(0 && "not implemented here");
}		  

void sFactory::separateVarsXSYZ( OoqpVector& x_in, OoqpVector& s_in, 
		  OoqpVector& y_in, OoqpVector& z_in, OoqpVector& vars_in)
{
  assert(0 && "not implemented here");
}		  


void sFactory::iterateStarted()
{
  iterTmMonitor.recIterateTm_start();
  tree->startMonitors();
}


void sFactory::iterateEnded()
{
  tree->stopMonitors();

  if(tree->balanceLoad()) {
    // balance needed
    data->sync();
      
    for(size_t i=0; i<registeredVars.size(); i++)
      registeredVars[i]->sync();
    
    resid->sync();
    
    linsys->sync();

    printf("Should not get here! OMG OMG OMG\n");
  }
  //logging and monitoring
  iterTmMonitor.recIterateTm_stop();
  m_tmTotal += iterTmMonitor.tmIterate;

  if(tree->rankMe==tree->rankZeroW) {
#ifdef TIMING
    extern double g_iterNumber;
    printf("TIME %g SOFAR %g ITER %d\n", iterTmMonitor.tmIterate, m_tmTotal, (int)g_iterNumber);
    //#elseif STOCH_TESTING
    //printf("ITERATION WALLTIME: iter=%g  Total=%g\n", iterTmMonitor.tmIterate, m_tmTotal);
#endif
  }
}

