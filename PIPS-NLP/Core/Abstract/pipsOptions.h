/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef PIPSOPTIONS_H
#define PIPSOPTIONS_H

/** Structure for the options that can be set through the OOPS control file. */
class pipsOptions
{
 public:

  static pipsOptions	*defOpt;

  /* --- adding slacks to get parallelism setting:  */
  int AddSlackParallelSetting;

  /* ----------------------- General options ---------------------------- */
  /** Printing level  */
  int prtLvl;

  /** Iteration limit */
  int max_iter;

  /** Convergence tolerance */
  double conv_tol;

  /* -------- Options max no of iterative refinement ---------- 
    SymLinearSolver   	= 	(1)	MA57 	
    				  	=	2	PARDISO
    					=	6	Umfpack	*/

  int SymLinearSolver;

  /* -------- Options about using PETSC ---------- */
  int UsePetsc;
  int User_Defined_PC;
  int User_Defined_SymMat;
  int UsePetscOuter;
  int SCOPF_precond;


  /* -----------------------  options for iterative refinement ---------------------------- */

  /* -------- do iterative refinement ---------- 
    DoIR_Aug    = 	0	default 	*/
  int DoIR_Aug;

  /* -------- do iterative refinement ---------- 
    DoIR_Full    = 	0	default 	*/
  int DoIR_Full;

  /* -------- Options max no of iterative refinement ---------- 
    MaxIR    = 	8	default 	*/
  int MaxIR;

  /* -------- tol of iterative refinement ---------- 
    IRtol    = 	1e-8	default 	*/
  double IRtol;




  /* -------- do tiny step test  -------- 
		  DoTinyStepTest 	=    	(0)	  
					  		(1) do it*/
  int DoTinyStepTest;

  /* -------- assume Mat is always singular  once detected-------- 
		   AssuneMatSingular	 =		(0)   
							 		(1) do it*/
  int AssumeMatSingular;




  int outerSolve; //  0: Default solve - Schur complement based decomposition  
  					  // 1: Iterative refinement 
  					  // 2: BiCGStab
  					  // 3: Default solve - Schur complement based decomposition, do not compress!

  // 0: default solve - add diag part (X^{-1}Z) to Q
  // 1: separate them : FIXME_NY: now only works if outerSolve =3
  int splitHesDiag; 				
						
  /* -------- about schur complement solver --------  */
  //(0): do not build SC   1: use dense Schur	2: use sparse Schur
  // 3: compute SC in sparse triplet format (and use MUMPS as a parallel solver); ignores 'SolveSchurScheme' below
  int BuildSchurComp; 
  // Number of MUMPS ranks
  int MUMPSranks; 
  double AbsTolForZero;

  /* -------- how to compute Schur --------  */
  int SolveSchurScheme; //   (0): LDLt  2: BICG

  /* -------- about Reduced space solver --------  */
  int UseReducedSpace;		//   (0): full sapce  1: reduced space
  int RS_SchurSolver;	 	//   (0): build schur from full sapce	1: rfrom educed space
  int RS_MaxIR; 				// do IR from LU solver	(0): not do IR  
  double RS_LU_PivotLV;

  /* -------- about dWd test --------  */
  int dWd_test; 				//   (0): not applied  1: applied dwd test 	2: do check in advance
  int dWd_test_soc; 			//   (0): not applied  1: applied dwd test


  /* -------- use filter --------  */
  int UseFilter; 			//   0: not applied  (1): applied filter

  int FilterResetStep;		//   0: not applied  (5): 5 rejection

  /* -------- ma57 parameter--------  */
  double HSL_PivotLV;		//  close to zero=fast, close to 0.5=stable   (1e-4)    1e-8 cannot solve pdegas!
  int MA57_Ordering;		// 5 automatic choice(MA47 or Metis); 4 use Metis (ND); 3 min degree ordering as in MA27; 2 use MC47; 

  /* -------- max number of Line search --------  */
  int LineSearchMatStep;

  /* Use Partitioning Algorithm*/
  int NP_Alg;


  /* this is the constant used in test dwd >= kappa_tWt d'd*/
  double kappa_tWt;

  /*   use mu  in the test dwd >= kappa_tWt * mu * d'd or not */
  int kappaWithMu;


  /* check constraint violation in switching condition*/
  int CheckSmallConstVio;

  /* do second order correction or not */
  int DoSOC;



  /* -------- about regularization --------  */

  /* use different method to compute dual regularization
	* 0: compute dual regularization from costant*\mu
	* 1: use VZ's method: set dual regualrization from constraint violation. here we need to change the filter tests
  */
  int UseDualRegAlg;


  /* -------- is NLP or QP/LP--------  */
  int isNLP; 			//   0: is QP/LP  (1): is NLP

  /* ============================ methods =============================== */

  /** base constructor */
  pipsOptions();

  /** trivial deconstructor */
  ~pipsOptions();

  /** Read the PIPS control file. */
  void readFile(void);

  /** parse a single line of the file. Return if option found or not */
  bool parseLine(char *line);

  void copyFrom(pipsOptions &pipsOpt);

  /** print Option settings to screen */
  void print();


  /** define global option */
  void defGloOpt();

  
};

#endif

