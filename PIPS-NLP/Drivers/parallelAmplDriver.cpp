/* PIPS-NLP                                                         	*
 * Author: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "amplGenStochInput.hpp"
#include "amplGenStochInput_AddSlack.hpp"

#include "NlpPIPSIpmInterface.h"

#include "sFactoryAug.h"

#include "FilterIPMStochSolver.h"

#include "sNlpInfoFromNL.h"

#include "pipsOptions.h"

#ifdef WITH_PETSC
#include "petscksp.h"
#endif


#include "mpi.h"


#ifdef TIMING
  double timeFromAMPL;
  double probGenTime;
  double PartSolver_GenTime;
  double PartSolver_SolTime;  
  double PartSolver_FactTime;
  int call_sol_Times;
  int call_fact_Times;

  int call_sol_Times_MA57;
  int call_fact_Times_MA57;  
  double genTime_localAmpl;
#endif

using namespace std;

extern int gInnerSCsolve;
extern int gNP_Alg;
extern int gAddSlackParallelSetting;
extern int gSymLinearSolver;
extern int gUseReducedSpace;

// using the latest interface PIPSIpmInterface 
int solve               (const string& datarootname, int nscen);


int main(int argc, char ** argv) {
#ifdef WITH_PETSC
	PetscInitialize(NULL,NULL,"petsc_option.Opt",NULL);
#endif
	
  MPI_Init(&argc, &argv);
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);

  if(argc<3) {
    if (mype == 0) printf("Usage: %s [rawdump root name] [num scenarios] [solution output root name]\n",argv[0]);
    return 1;
  }
  
  if (mype == 0) cout << argv[0] << " starting..." << endl;  
  string datarootname(argv[1]);
  int nscen = atoi(argv[2]);

#ifdef TIMING
	double tTot=MPI_Wtime();
	timeFromAMPL = 0.0;
	probGenTime  = 0.0;
	PartSolver_GenTime  = 0.0;
	PartSolver_SolTime  = 0.0;
	PartSolver_FactTime = 0.0;	
	call_sol_Times=0;
	call_fact_Times=0;	
	call_sol_Times_MA57=0;
	call_fact_Times_MA57=0;
	genTime_localAmpl=0.0;
#endif

  solve(datarootname, nscen);

#ifdef TIMING
  tTot = MPI_Wtime()-tTot;
  MPI_Barrier(MPI_COMM_WORLD);
  if(0==mype){
	  cout << " " << endl;  	
      cout << "Total Running Time " << tTot << endl;
  }
#endif


  MPI_Finalize();
  return 0;
}


int solve(const string& datarootname, int nscen)
{
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  int nprocs; MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  if(0==mype) cout << "Using a total of " << nprocs << " MPI processes." << endl;

#ifdef TIMING
	  double tGenTime=MPI_Wtime();
	  double tGenTime2,tGenTime3;
#endif

  pipsOptions *pipsOpt = new pipsOptions();
  pipsOpt->readFile();
  pipsOpt->defGloOpt();
  				  
  gInnerSCsolve=0;

  amplGenStochInput *s;
  if(gAddSlackParallelSetting==0)
    s = new amplGenStochInput(datarootname,nscen,MPI_COMM_WORLD);
  else if(gAddSlackParallelSetting==1)
    s = new amplGenStochInput_AddSlack(datarootname,nscen,MPI_COMM_WORLD);


  if (mype == 0) cout << "AMPL NL Input created .." << endl;
  
#ifdef TIMING
	
	MPI_Barrier(MPI_COMM_WORLD);
    tGenTime2 = MPI_Wtime();
	
    if(0==mype){
			cout << " " << endl;
			cout << "Problem Generation Time 1 " << tGenTime2-tGenTime << endl;
			cout << " " << endl;
    }
#endif

  NlpPIPSIpmInterface<sFactoryAug, FilterIPMStochSolver,  sNlpInfoFromNL> pipsIpm(*s);
	
#ifdef TIMING
    MPI_Barrier(MPI_COMM_WORLD);
    tGenTime3 = MPI_Wtime();
    if(0==mype){
		cout << " " << endl;
		cout << "Problem Generation Time 2 " << tGenTime3-tGenTime2 << endl;
		cout << "    in which, GenTime on local NL file:" << genTime_localAmpl << endl;
		cout << " " << endl;
    }
#endif


  if (mype == 0) cout << "PIPSIpmInterface created .." << endl;
//  delete s;
  if (mype == 0) cout << "AMPL NL  deleted ... solving" << endl;
 
  if (mype == 0){
    cout 	<< "  \n  -----------------------------------------------\n" 
     		<< "  NLP Solver \n" 
    		<< "  Nai-Yuan Chiang & V.M. Zavala, Argonne National Laboratory, 2013\n" 
    		<< "  -----------------------------------------------\n" << endl;

    if(gUseReducedSpace>0)
  	  cout << "\n  Reduced Space Solver ------	 Reduced Space Solver with Umfpack and following linear solver.\n";

    if(0==gSymLinearSolver)
	  cout << "\n  Linear system solver ------	 Ma27.\n\n";
    else if(1==gSymLinearSolver)
	  cout << "\n  Linear system solver ------	 Ma57.\n\n";
    else if(2==gSymLinearSolver)
	  cout << "\n  Linear system solver ------	 Pardiso.\n\n";
    else if(3==gSymLinearSolver)
	  cout << "\n  Linear system solver ------	 Umfpack.\n\n";
  } 

#ifdef TIMING
		double tSolTime=MPI_Wtime();
		MPI_Barrier(MPI_COMM_WORLD);
		tGenTime  = MPI_Wtime()-tGenTime;	
		MPI_Barrier(MPI_COMM_WORLD);	
#endif
	
  pipsIpm.go();

#ifdef TIMING
		tSolTime = MPI_Wtime()-tSolTime;
		MPI_Barrier(MPI_COMM_WORLD);
		if(gNP_Alg!=0)probGenTime += PartSolver_GenTime;
		
		if(0==mype){
			cout << " " << endl;
			cout << "Problem Generation Time " << tGenTime << endl;
			cout << "Total Gen Time (with Reduced Solver)" << tGenTime+probGenTime << endl;
			cout << "Total Solve Time " << tSolTime << endl;
	        cout << "Total AMPL Time " << timeFromAMPL << endl;		
			cout << "Call Fact Times _MA57 " << call_fact_Times_MA57<< endl;				
			cout << "Call Sol Times _MA57 " << call_sol_Times_MA57 << endl;		
		}
		if(gNP_Alg!=0){
			cout << " " << endl;
			cout << "PART_ALG Generation Time: " << PartSolver_GenTime << endl;
			cout << "PART_ALG Factorization Time:" << PartSolver_FactTime << endl;
			cout << "PART_ALG Solve Time:" << PartSolver_SolTime << endl;
			cout << "Call Fact Times  " << call_fact_Times<< endl;				
			cout << "Call Sol Times  " << call_sol_Times << endl;					
		}
#endif
  pipsIpm.writeSolution();
  delete pipsOpt;
  delete s;

  return 0;
}




