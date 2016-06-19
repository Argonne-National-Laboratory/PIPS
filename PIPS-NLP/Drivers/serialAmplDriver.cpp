/* PIPS-NLP                                                             
 * Author: Nai-Yuan Chiang
 * (C) 2015 Argonne National Laboratory
 */

#include "NlpGenData.h"
#include "NlpGenVars.h"
#include "NlpGenResiduals.h"
#include "NlpInfoAMPL.h"


#include "FilterIPMSolver.h"

#include <iostream>
#include <string.h>
#include <stdio.h>


#include <cstdio>
#include <cassert>

#include <asl_pfgh.h>
#include "getstub.h"

#include "NlpGenSparseWithSolver.h"
#include "cNlpGenSparseNLP.h"


#include "Status.h"

#include "getAmplFunction.h"
#include "pipsOptions.h"


#include "AmplData_NL.hpp"


extern int gSymLinearSolver;
extern int gBuildSchurComp;
extern int gUseReducedSpace;

#ifdef WITH_PETSC
#include "petscksp.h"
#endif


#ifdef TIMING
  #include "mpi.h"
  double timeFromAMPL;
  double probGenTime;
  int call_sol_Times_MA57;
  int call_fact_Times_MA57;
#endif


void LocalNegate(double * x, int nx)
{
  for (int i = 0;  i < nx;  i++)
    x[i] = -x[i];
}


int main( int argc, char * argv[] )
{

#ifdef WITH_PETSC
PetscInitialize(NULL,NULL,"petsc_option.Opt",NULL);
MPI_Init(&argc, &argv);
#endif


#ifdef TIMING
	MPI_Init(&argc, &argv);
	double tTot=MPI_Wtime();
	timeFromAMPL = 0.0;
	probGenTime  = 0.0;
#endif
  
  pipsOptions *pipsOpt = new pipsOptions();
  pipsOpt->readFile();
  pipsOpt->defGloOpt();
  if(gUseReducedSpace!=0) 
  	gUseReducedSpace=1;
  				  
  int NlpPrintLevel = 0;

  double  *Objgrad     = NULL,  *xlow   = NULL,  *xupp = NULL;
  int      nx;
  int     *irowQ = NULL,  *jcolQ  = NULL;
  int      nnzQ;
  int     *irowA = NULL,  *jcolA  = NULL;
  int      nnzA, my;
  int     *irowC = NULL,  *jcolC  = NULL;
  int      nnzC, mz;
  int 	   result;

  int nxL=0, nxU=0, nsL=0, nsU=0;
  int nnzCL, nnzCU;
  int nnzCFull, nnzCFull_Jac;
  int AMPLnnzQ,AMPLnnzA,AMPLnnzC;

  double  *dQ    = NULL,  *dA     = NULL,  *dC   = NULL;
  double  *bA     = NULL,  *clow   = NULL,  *cupp = NULL;

  char    *ixlow = NULL,  *ixupp  = NULL;
  char    *iclow = NULL,  *icupp  = NULL;

  double  *x     = NULL,  *gamma  = NULL,  *phi  = NULL;
  double  *y     = NULL;
  double  *z     = NULL,  *lambda = NULL,  *pi   = NULL;
  int ierr;
  int full_size_Hessian = 0;
  double * HessianElts = NULL;

  // Allocate the ampl solver library (ASL) context.
  AmplData_NL *localNL = new AmplData_NL();

  // Add the suffix   
  AmplSuffix *amplSuffix = new AmplSuffix();
  amplSuffix->DefineSuffix("pipsNLP_DecisionVar_in", Suffix_Var, Suffix_Int);  

  // Allocate the ampl solver library (ASL) context.
  localNL->initASL(argv, amplSuffix);	
  ASL_pfgh * asl = localNL->locASL;
  localNL->initialize(n_var,n_con);


  int *decisionVarIDX 	 = NULL;
  int decisionVarDim	 = 0;
  int *schurVarConIDinNL = NULL;

  if(0!=gBuildSchurComp && 0!=gUseReducedSpace){
  	decisionVarIDX = amplSuffix->GetSuffixVal_Int(asl, "pipsNLP_DecisionVar_in", Suffix_Var);
	//find correct index in C:	in ampl model, index of decision var starts from 1, here we correct it as zero;
	for(int j=0; j<n_var; j++){
	  decisionVarIDX[j] -=1;
	  if(decisionVarIDX[j]>=0) 
		decisionVarDim++;
	}
	schurVarConIDinNL = (int*)malloc(decisionVarDim*sizeof(int));

	for(int findSCVar=0; findSCVar<decisionVarDim;findSCVar++){
	  for(int j=0; j<n_var; j++){
	    if(decisionVarIDX[j] == findSCVar){
		  schurVarConIDinNL[findSCVar]=j;
		  break;
	    }
	  }
	}
  }

  int *irow, *jcol;
  int *rowMap;

  // Now initialise the Hessian. 
  assert(0==full_size_Hessian);
  nnzQ = ampl_get_nnz_Hessian_Tri();
  if (0==nnzQ){
	printf("This is LP problem, please use OOQP in LP/QP version.");
  }
	
  irow = (int*)sputinfo->hrownos;
  jcol = (int*)sputinfo->hcolstarts;
  
  // Count the problem sizes, and create a map between the rows of the
  // Jacobian and A and C (A is the Jacobian of the equality constraints,
  // C is the Jacobian of the inequality constraints.)
  // rowmap is moved into getAmplFunction as global variable

  ampl_count_sizes_SplitSlack( irow, jcol, nx, nnzQ, my, nnzA, mz, nnzC,full_size_Hessian,
   			nnzCL, nnzCU, nxL, nxU, nsL, nsU);

  // Note that n_con and n_var are MACROS the translate to
  // asl->n_con and asl->n_var. Arghh....
  assert( my + mz == n_con );
  assert( nnzA + nnzC == nzc);


  HessianElts = (double*) malloc(nnzQ*sizeof(double));


  // Allocate the structure of the NLP
  newNlpGenSparse( &Objgrad, 	nx,
		&irowQ, nnzQ,  &jcolQ,	&dQ,
		&xlow,		   &ixlow,	&xupp, &ixupp,
		&irowA, nnzA,  &jcolA,	&dA,
		&bA, 	my,
		&irowC, nnzC,  &jcolC,	&dC,
		&clow,	mz,    &iclow,	&cupp, &icupp,
		&ierr );

  if( ierr != 0 ) {
    fprintf( Stderr, "Couldn't allocate enough memory\n" );
    exit( 1 );
  }

    
  double *dxWrk = (double*) malloc (n_var*sizeof(double));
  ampl_get_InitX0(dxWrk);
  double init_Obj = ampl_get_Obj(dxWrk);
	
  ampl_get_ObjGrad(dxWrk, Objgrad);
  ampl_get_bounds(
		   xlow, nx, ixlow, xupp, ixupp,
		   bA, my,
		   clow, mz, iclow, cupp, icupp );

  ampl_get_matrices( irow, jcol, HessianElts,
		   nx, nnzQ, my, nnzA, mz, nnzC,
		   irowQ, jcolQ, dQ,
		   irowA, jcolA, dA,
		   irowC, jcolC, dC, 
		   dxWrk, 
		   full_size_Hessian );

  free(dxWrk);
  free(HessianElts);

  double objectiveConstant = objconst(0);


  if (objtype[0] == kMaximize) {
    LocalNegate(Objgrad, nx);
    LocalNegate(dQ, nnzQ);
  }


  double objectiveValue=0;

  printf("  \n  -----------------------------------------------\n");
  printf("  NLP Solver \n");
  printf("  Nai-Yuan Chiang & V.M. Zavala, Argonne National Laboratory, 2013\n");
  printf("  -----------------------------------------------\n");
  { 
	printf("  Variables ================ %-5d \n",nx);
	printf("  Equality Constraints ===== %-5d \n",my);
	printf("  Inequality Constraints === %-5d \n",mz);  
	if(xlow)
	  printf("  Variable Lower Bounds ==== %-5d \n",nxL);
	if(xupp)
	  printf("  Variable Upper Bounds ==== %-5d \n",nxU);
	if(clow)
	  printf("  Inequality Lower Bounds == %-5d \n",nsL);
	if(cupp)
	  printf("  Inequality Upper Bounds == %-5d \n\n",nsU);
  }


  //create the NLP formulation factory
  NlpGenSparseWithSolver* nlp  = NULL;
  NlpGenData	  * prob = NULL;
  NlpGenVars     	* vars;
  FilterIPMSolver  * s;
  NlpGenResiduals  	* resid;

  NlpInfoAMPL *updateNlpAmpl = new NlpInfoAMPL(nx,my,mz,nnzQ,nnzA,nnzC,nxL, nxU, nsL, nsU);

  if(0==gSymLinearSolver)
	std::cout << "\n  Linear system solver ------	 Ma27.\n\n";
  else if(1==gSymLinearSolver)
	std::cout << "\n  Linear system solver ------	 Ma57.\n\n";
  else if (2==gSymLinearSolver)
	std::cout << "\n  Linear system solver ------	 Pardiso.\n\n";
  else if (3==gSymLinearSolver)
	std::cout << "\n  Linear system solver ------	 UMFPACK.\n\n";
  else{
	std::cout << "\n  Linear system solver ------	 unknown linear solver! \n\n";  
	assert("No linear solver!" && 0);
  }

  nlp = new NlpGenSparseWithSolver( nx, my, mz, nnzQ, nnzA, nnzC);
  prob = (NlpGenData * ) nlp->copyDataFromSparseTriple(
								  Objgrad,		
								  irowQ,  nnzQ,   jcolQ,  dQ,
								  xlow,   ixlow,  xupp,   ixupp,
								  irowA,  nnzA,   jcolA,  dA,	  bA,
								  irowC,  nnzC,   jcolC,  dC,
								  clow,   iclow,  cupp,   icupp , rowMap,
								  nxL, nxU, nsL, nsU, updateNlpAmpl);

  if(0!=gBuildSchurComp|| 0!=gUseReducedSpace){
  	((NlpGenData *) prob)->schurVarConID = schurVarConIDinNL;
	((NlpGenData *) prob)->schurSize 	 = decisionVarDim;
  }
  
  vars	= (NlpGenVars * )  nlp->makeVariables( prob);
  resid = (NlpGenResiduals* ) nlp->makeResiduals( prob );
  s 	= new FilterIPMSolver( nlp, prob );
	

  prob->nxLOri = nxL;
  prob->nxUOri = nxU;
  prob->nsLOri = nsL;
  prob->nsUOri = nsU;

  s->monitorSelf();


#ifdef TIMING
	double tGenTime = MPI_Wtime()-tTot;
#endif

#ifdef TIMING
	double tSolTime = MPI_Wtime();
#endif

  result = s->solve(prob,vars, resid);

#ifdef TIMING
	tSolTime = MPI_Wtime()-tSolTime;
#endif


  if(result !=0 ){
	std::cout << "Could not solve the problem.\n";
  } else {
    if( NlpPrintLevel >= 2) {

      objectiveValue = prob->objectiveValue(vars);
	 
	  if (objtype[0] == kMaximize)
	    objectiveValue = -objectiveValue;

      std::cout.precision(4);
	  std::cout << "Optimal Objective: " << objectiveValue << std::endl;

	  std::cout << "Solution: \n"; 
      vars->x->writefToStream( std::cout, "x[%{index}] = %{value}" );
    } 
	
  }

  double *solx = new double[nx];
  double *soly = new double[my];
  double *solz = new double[mz];
  vars->x->copyIntoArray(solx);
  vars->y->copyIntoArray(soly);
  vars->z->copyIntoArray(solz);
  ampl_write_solution( solx, soly, solz);
  delete []solx;
  delete []soly;
  delete []solz;

#ifdef TIMING
		tTot = MPI_Wtime()-tTot;
		std::cout << "\n" << std::endl;
		std::cout << "Problem Generation Time " << tGenTime << std::endl;
		if(gUseReducedSpace!=0) std::cout << "Total Problem Generation Time (With Saddle-Point Solver))" << tGenTime + probGenTime<< std::endl;
		std::cout << "Total Solve Time " << tSolTime << std::endl;
		std::cout << "Total Running Time " << tTot << std::endl;
		std::cout << "Total AMPL Time " << timeFromAMPL << std::endl;
		std::cout << "\n" << std::endl;
#endif


  freeNlpGenSparse(&Objgrad,
  	&irowQ,&jcolQ,&dQ,
  	&xlow,&ixlow,&xupp,&ixupp,
  	&irowA,&jcolA,&dA,
  	&bA,
    &irowC,&jcolC,&dC,
	&clow,&iclow,&cupp,&icupp);

  ampl_free_mapinfo();

  if(updateNlpAmpl) delete updateNlpAmpl;
  delete pipsOpt;
  delete amplSuffix;
  delete localNL;


  delete s;
  delete vars;  
  delete resid;
  delete prob;
  delete nlp;
  
};






