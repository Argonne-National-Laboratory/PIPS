/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */
#include <stdio.h>
#include <stdlib.h>

#include "rawInput.hpp"
#include "PIPSIpmInterface.h"

#include "sFactoryAugSchurLeaf.h"
#include "MehrotraStochSolver.h"

#include <string>
#include <sstream>

#ifdef TIMING_FLOPS
extern "C" {
    void HPM_Init(void);
    void HPM_Start(char *);
    void HPM_Stop(char *);
    void HPM_Print(void);
    void HPM_Print_Flops(void);
    void HPM_Print_Flops_Agg(void);
    void HPM_Terminate(char*);
}
#endif

int main(int argc, char ** argv) {
  MPI_Init(&argc, &argv);
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);

  if(argc<3) {
    if (mype == 0) printf("\nUsage:\n%s   [rawdump root name]   [num scenarios]   [outer solve (optional): 0 vanilla direct (default), 1 with iter.refin, 2 with BICGStab]   [inner solve (optional): 0 vanila direct (default), EXPERIMENTAL-> 1 iter.refin, 2. BiCGStab]\n\n",argv[0]);
    return 1;
  }
  
  string datarootname(argv[1]);
  int nscen = atoi(argv[2]);

  int outerSolve=2;
  if(argc>=4) {
    outerSolve = atoi(argv[3]);
    if(mype==0) cout << "Using option [" << outerSolve << "] for outer solve" << endl;
  }

  int innerSolve=0;
  if(argc>=5) {
    innerSolve = atoi(argv[4]);
     if(mype==0) cout << "Using option [" << innerSolve << "] for inner solve" << endl;
  }

  if(mype==0) cout << argv[0] << " starting ..." << endl;
  int nprocs; MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if(0==mype) cout << "Using a total of " << nprocs << " MPI processes." << endl;

  rawInput* s = new rawInput(datarootname,nscen);
  if(mype==0) cout <<  " raw input created from " << datarootname<< endl;
  PIPSIpmInterface<sFactoryAugSchurLeaf, MehrotraStochSolver> pipsIpm(*s);

  pips_options::setIntParameter("OUTER_SOLVE", outerSolve);
  pips_options::setIntParameter("INNER_SC_SOLVE", innerSolve);

  if(mype==0) cout <<  "PIPSIpmInterface created" << endl;
  delete s;
  if(mype==0) cout <<  "rawInput deleted ... starting to solve" << endl;


#ifdef TIMING_FLOPS
  if(mype==0) cout << "FLOPS are being recorded using HPM library." << endl;
  HPM_Init();
  //HPM_Start("PIPSTotFlops");
#endif

  pipsIpm.go();

#ifdef TIMING_FLOPS
  //HPM_Stop("PIPSTotFlops");
  HPM_Print_Flops();
  HPM_Print_Flops_Agg();
#endif
//   if(mype==0) cout << "solving done" << endl;

//   if(mype==0) cout << "Saving solution" << endl;
//   for(int s=0; s<nscen; s++) {

//     std::vector<double> duals = pipsIpm.getSecondStageDualRowSolution(s);
//     if(duals.size()) {
//       stringstream ss1; ss1<<"out_duals_scen"<<(s+1)<<".txt";
//       cout << "saving duals to " << ss1.str() << endl;
//       ofstream fileduals(ss1.str().c_str());
      
//       for(size_t i=0; i<duals.size(); i++)
// 	fileduals << duals[i] << endl;
//       fileduals.close();
//     }

//     std::vector<double> primals = pipsIpm.getSecondStagePrimalColSolution(s);
//     if(primals.size()) {
//       stringstream ss2; ss2<<"out_primals_scen"<<(s+1)<<".txt";
//       cout << "saving primals to " << ss2.str() << endl;
//       ofstream fileprimals(ss2.str().c_str());
//       std::vector<double> primals = pipsIpm.getSecondStagePrimalColSolution(s);
//       for(size_t i=0; i<primals.size(); i++)
// 	fileprimals << primals[i] << endl;
//       fileprimals.close();
//     }
//   }
//   if(mype==0) {
//     std::vector<double> firstStageSol = pipsIpm.getFirstStagePrimalColSolution();
//     string sFile = "out_primal_1stStage.txt";
//     cout << "saving 1st stage sol to " << sFile << endl;
//     ofstream file1stStg(sFile.c_str());
//     for(size_t i=0; i<firstStageSol.size(); i++)
//       file1stStg << firstStageSol[i] << endl;
//     file1stStg.close();
//   }
//   cout << "Solution saved" << endl;
  MPI_Finalize();
  return 0;
}

