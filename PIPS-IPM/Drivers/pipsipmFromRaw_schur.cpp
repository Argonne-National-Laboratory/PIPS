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

using namespace std;


int main(int argc, char ** argv) {
  MPI_Init(&argc, &argv);
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);

  if(argc<3) {
    if (mype == 0) printf("Usage: %s [rawdump root name] [num scenarios] [solution output root name] [eq mult root name]\n",argv[0]);
    return 1;
  }
  
  string datarootname(argv[1]);
  int nscen = atoi(argv[2]);

  rawInput* s = new rawInput(datarootname,nscen);
  PIPSIpmInterface<sFactoryAugSchurLeaf, MehrotraStochSolver> pipsIpm(*s);
  delete s;

  pipsIpm.go();

  for(int s=0; s<nscen; s++) {

    stringstream ss1; ss1<<"out_duals_scen"<<(s+1)<<".txt";
    cout << "saving duals to " << ss1.str() << endl;
    ofstream fileduals(ss1.str().c_str());
    std::vector<double> duals = pipsIpm.getSecondStageDualRowSolution(s);
    for(size_t i=0; i<duals.size(); i++)
      fileduals << duals[i] << endl;
    fileduals.close();

    stringstream ss2; ss2<<"out_primals_scen"<<(s+1)<<".txt";
    cout << "saving primals to " << ss2.str() << endl;
    ofstream fileprimals(ss2.str().c_str());
    std::vector<double> primals = pipsIpm.getSecondStagePrimalColSolution(s);
    for(size_t i=0; i<primals.size(); i++)
      fileprimals << primals[i] << endl;
    fileprimals.close();
  }

  MPI_Finalize();
  return 0;
}

