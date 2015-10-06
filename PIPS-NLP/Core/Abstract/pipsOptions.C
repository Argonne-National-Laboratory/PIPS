/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cassert>

#include "pipsOptions.h"
#include <mpi.h>


int gDoIR_Aug;
int gDoIR_Full;
int gMaxIR;
double gIRtol;

double gconv_tol;
int gmax_iter;

int separateHandDiag;
int gSymLinearSolver;
int gSymLinearAlgSolverForDense;


int gBuildSchurComp;
int gSolveSchurScheme;
int gUseReducedSpace;
int gRS_SchurSolver;
int gRS_MaxIR;
int gPipsPrtLV;
int gdWd_test;
int gUseFilter;
int gDoTinyStepTest;
int gAssumeMatSingular;
int gFilterResetStep;
int gLineSearchMatStep;
int gdWd_test_soc;

int gNP_Alg;

double gRS_LU_PivotLV;

double gHSL_PivotLV;

double gkappa_tWt;

int gCheckSmallConstVio;

int gDoSOC;
int gkappaWithMu;

extern int gOuterSolve;

int gUsePetsc;
int gUser_Defined_PC;
int gUser_Defined_SymMat;
int gUsePetscOuter;
int gSCOPF_precond;

int gAddSlackParallelSetting;

int gUseDualRegAlg;
int gMA57_Ordering;
int gisNLP;

#include "constants.h"
PreCondInfo *preCond;


#ifndef FindMPI_ID
#define FindMPI_ID(FLAG,MYID) MPI_Initialized(&FLAG); \
		if(FLAG) MPI_Comm_rank(MPI_COMM_WORLD,&MYID); \
		else MYID=0;
#endif

#ifndef FindMPI_Size
#define FindMPI_Size(FLAG,SIZE) MPI_Initialized(&FLAG); \
		if(FLAG) MPI_Comm_size(MPI_COMM_WORLD,&SIZE); \
		else SIZE=1;
#endif


//pipsOptions *glopt;
pipsOptions* pipsOptions::defOpt = NULL;

/* ----------------------------------------------------------------------------
 pipsOptions::pipsOptions - Constructor
---------------------------------------------------------------------------- */
pipsOptions::pipsOptions()
  :	prtLvl(1),
	outerSolve(3),
	splitHesDiag(0),
    conv_tol(1.e-6),
	max_iter(500),	
	AddSlackParallelSetting(0),
    DoIR_Aug(1),
    DoIR_Full(0),
	MaxIR(10),
    IRtol(1e-12),
    SymLinearSolver(1),
    HSL_PivotLV(1e-4),
    MA57_Ordering(5), 
	dWd_test(0),
	dWd_test_soc(0),
	kappa_tWt(1e-12),
	kappaWithMu(1),
	DoSOC(1),
	UseFilter(1),
	FilterResetStep(5),
	DoTinyStepTest(1),
	AssumeMatSingular(0),
	CheckSmallConstVio(1),
	LineSearchMatStep(50),
	UsePetsc(0),
	User_Defined_PC(2),
	User_Defined_SymMat(1),
	UsePetscOuter(1),
	UseReducedSpace(0),
	RS_SchurSolver(1),
	RS_MaxIR(5),
	RS_LU_PivotLV(0.0001),
   BuildSchurComp(1),
   SolveSchurScheme(0),
   NP_Alg(0),
   SCOPF_precond(0),
   UseDualRegAlg(0),
   isNLP(1)
{}

/* ----------------------------------------------------------------------------
 pipsOptions::~pipsOptions - Destructor
---------------------------------------------------------------------------- */
pipsOptions::~pipsOptions() {}

/* ----------------------------------------------------------------------------
 pipsOptions::copyFrom
---------------------------------------------------------------------------- */
void
pipsOptions::defGloOpt()
{

  gDoIR_Aug	 = DoIR_Aug; 
  gDoIR_Full = DoIR_Full; 
  
  if(gDoIR_Aug==1 || gDoIR_Full==1){
    gMaxIR =  MaxIR;  
    gIRtol =  IRtol; 
  }else{
	gMaxIR =	0;  
	gIRtol =	0; 
  }	

  gconv_tol			= conv_tol;
  gmax_iter			= max_iter;

  gSymLinearSolver 	= SymLinearSolver;
  if(gSymLinearSolver<=2){
  	gSymLinearAlgSolverForDense = gSymLinearSolver;
  }else 
  	gSymLinearAlgSolverForDense=1;
  
  separateHandDiag	= splitHesDiag;  			// 0: default solve - add diag part (X^{-1}Z) to Q
												// 1: separate them : FIXME_NY: now only works if outerSolve =3

  gOuterSolve = outerSolve; //  0: Default solve - Schur complement based decomposition  
  					  // 1: Iterative refinement 
  					  // 2: BiCGStab
  					  // 3: Default solve - Schur complement based decomposition, do not compress!

  gBuildSchurComp 	= BuildSchurComp;
  gSolveSchurScheme = SolveSchurScheme;			
  gUseReducedSpace  = UseReducedSpace;
  gRS_SchurSolver	= RS_SchurSolver;
  gRS_MaxIR			= RS_MaxIR;
  gRS_LU_PivotLV	= RS_LU_PivotLV;

  gPipsPrtLV 		= prtLvl;
  gdWd_test			= dWd_test;
  gdWd_test_soc		= dWd_test_soc;

  gUseFilter		= UseFilter;	
  gDoTinyStepTest		= DoTinyStepTest;
  gAssumeMatSingular	= AssumeMatSingular;
  gFilterResetStep		= FilterResetStep;
  gLineSearchMatStep	= LineSearchMatStep;

  gHSL_PivotLV			= HSL_PivotLV;

  gNP_Alg			= NP_Alg;

  gkappa_tWt = kappa_tWt;

  gCheckSmallConstVio	= CheckSmallConstVio;

  gDoSOC = DoSOC;

  gkappaWithMu = kappaWithMu;

  gUsePetsc = UsePetsc;
  gUser_Defined_PC = User_Defined_PC;
  gUser_Defined_SymMat = User_Defined_SymMat;

  gUsePetscOuter = UsePetscOuter;

  gSCOPF_precond = SCOPF_precond;

  gAddSlackParallelSetting = AddSlackParallelSetting;


  gUseDualRegAlg = UseDualRegAlg;

  gisNLP = isNLP;
}


/* ----------------------------------------------------------------------------
 pipsOptions::copyFrom
---------------------------------------------------------------------------- */
void
pipsOptions::copyFrom(pipsOptions &os)
{
  this->prtLvl = os.prtLvl;
  this->max_iter = os.max_iter;
  this->conv_tol = os.conv_tol;
  this->MaxIR = os.MaxIR;

} 

/**
 *  Read the global options file.
 *
 *  All options are described by an identifier and a value.
 *  Characters after the value are ignored.
 *
 *  WSFSTWA  1  !Do Forward Sensitivity analysis of ten worst components
 *  WSFSTWTS 1  !And take step suggested by analysis
 */
void pipsOptions::readFile()
{
  FILE *optfile;
  char buffer[1000];
  
  int mype,doFlag;  FindMPI_ID(doFlag,mype);
  int prosize;	FindMPI_Size(doFlag,prosize);

  /** The first time that the global option file is read, these options
   *  are remembered in pipsOptions::defOpt */
  if (pipsOptions::defOpt==NULL){
    defOpt = this;
  }

  optfile = fopen("pipsnlp.parameter","r");
  if (optfile==NULL) {
    return;
  }
  
  /* Read one line of the options file */
  while(fgets(buffer, 999, optfile)!=NULL){
    parseLine(buffer);
  }
  fclose(optfile);

  if(splitHesDiag==1)
    assert(outerSolve==3);

  if(UseReducedSpace==1 && NP_Alg==0){
    BuildSchurComp = 1;
  }



  if (mype == 0){
  	printf("OPTION: Set printing level to %d\n",prtLvl);
	printf("OPTION: Set Iteration Limit to %d\n",max_iter);
	printf("OPTION: Set Convergence Tolerance to %.2e\n", conv_tol);

	printf("OPTION: Max Line search step:   %d\n",LineSearchMatStep);
  }

  if(AddSlackParallelSetting==1){
	if (mype == 0)printf("OPTION: Adding slacks in the parallel setting.\n");
  }

  if(UseReducedSpace==1){
//  	BuildSchurComp = 0;
	if(RS_MaxIR>10) RS_MaxIR=10;
	if (mype == 0)printf("OPTION: Use Reduced Space Solver.\n");
	if (mype == 0)printf("OPTION: Set Gen linear solver as UMFPACK.\n");
  }	

  if(SCOPF_precond==1){
  	if (mype == 0)printf("OPTION: Use SCOPF_preconditioner, only support serial run.\n");
  	assert(prosize==1);
//	UsePetscOuter=0;
//	User_Defined_PC=0;
//	User_Defined_SymMat=0;
  }
  
  if(UsePetsc==1){
  	assert(prosize==1);
  	if (mype == 0)printf("OPTION: Use Petsc iterative solver.\n");
	if(UsePetscOuter==1){
	  if (mype == 0)printf("OPTION: Use Petsc as outer linear system solver. Use user difined PC and Mat in petsc.. \n");
	  if(User_Defined_PC==0){
	  	if (mype == 0)printf("OPTION: Cannot use Petsc difined PC. Set User_Defined_PC=1. \n");
		User_Defined_PC=1;
	  }
	  User_Defined_SymMat=1;
	}

	if (User_Defined_PC==1 && mype == 0)printf("OPTION: Use user defined PC: use the full mat. \n");
	if (User_Defined_PC==2 && mype == 0)printf("OPTION: Use user defined PC: use the diagonals of H. \n");	
	if (User_Defined_SymMat==1 && mype == 0) printf("OPTION: Use user defined symmetric form.\n");
	if (User_Defined_PC==0 && User_Defined_SymMat!=0 && mype == 0){
	  printf("OPTION: Use petsc default form! set User_Defined_PC=0 and User_Defined_SymMat=0 \n");
	  assert("exit" && 0);
	}

	DoIR_Aug  = 0;
	DoIR_Full = 0;
	if (mype == 0)printf("OPTION: Use Petsc, switch off IR. \n");
  }




  if(HSL_PivotLV<=0)
    HSL_PivotLV=1e-8;
  else if(HSL_PivotLV>0.5)
  	HSL_PivotLV=0.5;
	
  if(this->SymLinearSolver==0){
	if (mype == 0)printf("OPTION: Set Sym linear solver as MA27.\n");
	if (mype == 0)printf("OPTION: MA27 pivot level = %.2e\n", HSL_PivotLV);	
  }else if(this->SymLinearSolver==1){
	if (mype == 0)printf("OPTION: Set Sym linear solver as MA57.\n");
	if (mype == 0)printf("OPTION: MA57 pivot level = %.2e\n", HSL_PivotLV);
    if(MA57_Ordering<=1 || MA57_Ordering>5)
  	  MA57_Ordering=5;	
	if (mype == 0)printf("OPTION: MA57 ordering method = %d, (see doc)\n", MA57_Ordering);
  }else if(this->SymLinearSolver==2){
	if (mype == 0)printf("OPTION: Set Sym linear solver as PARDISO.\n");
  }else if(this->SymLinearSolver==3){
	if (mype == 0)printf("OPTION: Set Sym linear solver as SaddlePoint.\n");
	if (mype == 0)printf("OPTION: LU solver max IR:  %d \n", RS_MaxIR);
	if (mype == 0)printf("OPTION: LU solver pivot tol:  %.2e \n", RS_LU_PivotLV);
  }else if(this->SymLinearSolver==4){
	if (mype == 0)printf("OPTION: Set Sym linear solver as PARDISO_Iter.\n");
  }else if(this->SymLinearSolver==5){
	if (mype == 0)printf("OPTION: Set Sym linear solver as PARDISO_Schur.\n");
  }else if(this->SymLinearSolver==6){
	if (mype == 0)printf("OPTION: Set Sym linear solver as Umfpack_symm.\n");
  }else {
	assert("need one sym solver"&&0);
  }

  if (DoIR_Aug != 0 ){
	if (mype == 0)printf("OPTION: Do IR on Aug sys.\n");
  }  
  if (DoIR_Full != 0 ){
	if (mype == 0)printf("OPTION: Do IR on Full sys.\n");
  }  
  if (DoIR_Aug != 0 || DoIR_Full !=0){
	if (mype == 0)printf("OPTION: Set MaxIR to %i\n", this->MaxIR);
	if (mype == 0)printf("OPTION: Set IRtol to %.2e\n", this->IRtol);
  }

  
  if(dWd_test==1){
  	if (mype == 0)printf("OPTION: Do dWd test. \n");
  }else if(dWd_test==2){
	if (mype == 0)printf("OPTION: Do dWd test, with checking switching condition in advance.\n");
  }else if(dWd_test==3){
	if (mype == 0)printf("OPTION: Do tWt test. Use d as direction\n");
  }else if(dWd_test==4){
	if (mype == 0)printf("OPTION: Do tWt test. Use n+t as direction\n");
  }else if(dWd_test==0){
	dWd_test_soc = 0;
  }

  if(dWd_test_soc==1 && dWd_test!=0 ){
  	if (mype == 0)printf("OPTION: Do dWd/tWt test (based on option dWd_test) for SOC step. \n");
  }

  if(dWd_test>=1){
  	if (mype == 0)printf("OPTION: set kappa_tWt (used in test dwd >= kappa_tWt d'd) : %.2e \n", kappa_tWt);
	if(kappaWithMu==1)
  	  if (mype == 0)printf("OPTION: add mu in the tWt test (used in test dwd >= kappa_tWt * mu * d'd) \n");		
  }

  if(CheckSmallConstVio!=0){
  	if (mype == 0)printf("OPTION: Require small constraint violation in SWC \n");
  }

  if(NP_Alg!=0){
  	if (mype == 0)printf("OPTION: Do Network partitioning. \n");
  }
  
  if(DoSOC == 0){
  	if (mype == 0)printf("OPTION: Skip Second Order Correction! \n");
  }

  if(UseDualRegAlg==1){
    if (mype == 0)printf("OPTION: Compute dual regularization with modified filter test.\n");
  }  

}

bool pipsOptions::parseLine(char *buffer)
{
  char label[100];
  double dval;
  
  bool found = false;
  sscanf(buffer, "%s%lf\n",label, &dval);  
  
  int mype,doFlag; FindMPI_ID(doFlag,mype);

  /* -----------------------------------------------------------------------
     General Options
     ----------------------------------------------------------------------- */
  /* Printing level  */
  if (strcmp(label, "prtLvl")==0){
//    if (mype == 0)printf("OPTION: Set printing level to %d\n",(int)dval);
    this->prtLvl = (int)dval;
    found = true;
  }    
  
  /* Set Iteration Limit */
  if (strcmp(label, "max_iter") == 0 && dval > 0.0) {
//    if (mype == 0)printf("OPTION: Set Iteration Limit to %d\n",(int)dval);
    this->max_iter = (int)dval;
    found = true;
  }    
  
  /* Set Convergence Tolerance */
  if (strcmp(label, "conv_tol") == 0 && dval > 0.0) {
//    if (mype == 0)printf("OPTION: Set Convergence Tolerance to %.2e\n", dval);
    this->conv_tol = dval;
    found = true;
  }    


  /* The value to use for ro_reg in num_factAS (pdreg in qnmfct) */
//  if (strcmp(label, "RO_REG") == 0 && dval >= 0.0) {
//    printf("OPTION: Use RO_REG = %g\n", dval);
//    this->ro_reg = dval;
//    found = true;
//  }    

  /* -----------------------------------------------------------------------
     Options About iterative refinement
    ----------------------------------------------------------------------- */

  if (strcmp(label, "DoIR_Aug") == 0 ) {
	if (mype == 0)printf("OPTION: Set DoIR_Aug to %i\n", int(dval));
	this->DoIR_Aug = int(dval);
	found = true;
  }

  if (strcmp(label, "DoIR_Full") == 0 ) {
	if (mype == 0)printf("OPTION: Set DoIR_Full to %i\n", int(dval));
	this->DoIR_Full = int(dval);
	found = true;
  }

  if (strcmp(label, "MaxIR") == 0 ) {
//	if (mype == 0)printf("OPTION: Set MaxIR to %i\n", int(dval));
	this->MaxIR = int(dval);
	found = true;
  } 

  if (strcmp(label, "IRtol") == 0 ) {
//	if (mype == 0)printf("OPTION: Set IRtol to %.2e\n", dval);
	this->IRtol = dval;
	found = true;
  }


  /* -----------------------------------------------------------------------
     Options About which linear solver to use
    ----------------------------------------------------------------------- */
  if (strcmp(label, "SymLinearSolver") == 0 ) {
	this->SymLinearSolver = int(dval);
	found = true;
  } 

  if (strcmp(label, "UsePetsc") == 0 ) {
	this->UsePetsc = int(dval);
	found = true;
  } 

  if (strcmp(label, "User_Defined_PC") == 0 ) {
	this->User_Defined_PC = int(dval);
	found = true;
  } 

  if (strcmp(label, "User_Defined_SymMat") == 0 ) {
	this->User_Defined_SymMat = int(dval);
	found = true;
  } 

  if (strcmp(label, "UsePetscOuter") == 0 ) {
	this->UsePetscOuter = int(dval);
	found = true;
  }

  if (strcmp(label, "SCOPF_precond") == 0 ) {
	this->SCOPF_precond = int(dval);
	found = true;
  }

  /* -----------------------------------------------------------------------
     Options About LP or NLP
    ----------------------------------------------------------------------- */
  if (strcmp(label, "outerSolve") == 0 ) {
	this->outerSolve = int(dval);
	if(this->outerSolve==1)
	  if (mype == 0) printf("OPTION: Split Hessian and Diag X^{-1}Z \n");	
	found = true;
  } 

  if (strcmp(label, "splitHesDiag") == 0 ) {
	this->splitHesDiag = int(dval);
	if(this->splitHesDiag==1)
	 if (mype == 0)  printf("OPTION: Split Hessian and Diag X^{-1}Z \n");
	found = true;
  }   

  /* -----------------------------------------------------------------------
     Options schur complement solver
    ----------------------------------------------------------------------- */
  if (strcmp(label, "BuildSchurComp") == 0 ) {
	this->BuildSchurComp = int(dval);
	found = true;
  } 
  
  if (strcmp(label, "SolveSchurScheme") == 0 ) {
	this->SolveSchurScheme = int(dval);
	found = true;
  } 


  /* -----------------------------------------------------------------------
     Reduced space solver
    ----------------------------------------------------------------------- */

  if (strcmp(label, "UseReducedSpace") == 0 ) {
	this->UseReducedSpace = int(dval);
	found = true;
  } 
  
  if (strcmp(label, "RS_SchurSolver") == 0 ) {
	this->RS_SchurSolver = int(dval);
	found = true;
  } 
  
  if (strcmp(label, "RS_MaxIR") == 0 ) {
	this->RS_MaxIR = int(dval);
	found = true;
  } 

  if (strcmp(label, "RS_LU_PivotLV") == 0 ) {
	this->RS_LU_PivotLV = (dval);
	found = true;
  } 

  /* -----------------------------------------------------------------------
     dWd test
    ----------------------------------------------------------------------- */
  if (strcmp(label, "dWd_test") == 0 ) {
	this->dWd_test = int(dval);
	found = true;
  } 
  if (strcmp(label, "dWd_test_soc") == 0 ) {
	this->dWd_test_soc = int(dval);
	found = true;
  } 


  /* -----------------------------------------------------------------------
     use tiny step check or not
    ----------------------------------------------------------------------- */
  if (strcmp(label, "UseFilter") == 0 ) {
	this->UseFilter = int(dval);
	found = true;
  }   

  /* -----------------------------------------------------------------------
     use filter or not
    ----------------------------------------------------------------------- */
  if (strcmp(label, "DoTinyStepTest") == 0 ) {
	this->DoTinyStepTest = int(dval);
	found = true;
  } 

  /* -----------------------------------------------------------------------
      assume Mat is singular
    ----------------------------------------------------------------------- */
  if (strcmp(label, "AssumeMatSingular") == 0 ) {
	this->AssumeMatSingular = int(dval);
	found = true;
  } 

  /* -----------------------------------------------------------------------
     reset Filter, this is the max number of previous iter rejected by filter 
    ----------------------------------------------------------------------- */
  if (strcmp(label, "FilterResetStep") == 0 ) {
	this->FilterResetStep = int(dval);
	found = true;
  } 

  

  /* -----------------------------------------------------------------------
     Options About MA57
    ----------------------------------------------------------------------- */
   if (strcmp(label, "HSL_PivotLV") == 0 ) {
	this->HSL_PivotLV = dval;
	found = true;
  }  
   if (strcmp(label, "MA57_Ordering") == 0 ) {
	this->MA57_Ordering = (int)dval;
	found = true;
  }     

  /* -----------------------------------------------------------------------
     Max step of line search
    ----------------------------------------------------------------------- */
  if (strcmp(label, "LineSearchMatStep")==0){
//    if (mype == 0)printf("OPTION: Max Line search step:   %d\n",(int)dval);
    this->LineSearchMatStep = (int)dval;
    found = true;
  } 


  /* -----------------------------------------------------------------------
     Use Partitioning Algorithm
    ----------------------------------------------------------------------- */
  if (strcmp(label, "NP_Alg")==0){
    this->NP_Alg = (int)dval;
    found = true;
  } 


  /* -----------------------------------------------------------------------
     this is the constant used in test dwd >= kappa_tWt d'd
    ----------------------------------------------------------------------- */
  if (strcmp(label, "kappa_tWt")==0){
    this->kappa_tWt = dval;
    found = true;
  }  
  
  /* -----------------------------------------------------------------------
     use mu  in the test dwd >= kappa_tWt * mu * d'd or not
    ----------------------------------------------------------------------- */
  if (strcmp(label, "kappaWithMu")==0){
    this->kappaWithMu = (int)dval;
    found = true;
  }   

  


  /* -----------------------------------------------------------------------
      check constraint violation in switching condition
    ----------------------------------------------------------------------- */
  if (strcmp(label, "CheckSmallConstVio")==0){
    this->CheckSmallConstVio = (int)dval;
    found = true;
  } 


  /* -----------------------------------------------------------------------
      do second order correction or not
    ----------------------------------------------------------------------- */
  if (strcmp(label, "DoSOC")==0){
    this->DoSOC = (int)dval;
    found = true;
  } 

  /* -----------------------------------------------------------------------
     use Carl's setting or not (adding slacks in the ampl model)
    ----------------------------------------------------------------------- */
  if (strcmp(label, "AddSlackParallelSetting")==0){
    this->AddSlackParallelSetting = (int)dval;
    found = true;
  }   


  /* -----------------------------------------------------------------------
     about regularization 
    ----------------------------------------------------------------------- */
  if (strcmp(label, "UseDualRegAlg")==0){
    this->UseDualRegAlg = (int)dval;
    found = true;
  }   


  
  return found;
}

void pipsOptions::print(){
  printf("iter_limit = %d\n",max_iter);
  printf("conv_tol = %f\n",conv_tol);

}

