/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#include "rawInput.hpp"
#include "sInterfaceCallbacks.h"
#include "PIPSIpmInterface.h"

#include "sFactoryAug.h"
#include "MehrotraStochSolver.h"

using namespace std;

// two ways of solving

// using the latest interface PIPSIpmInterface 
int solve               (const string& datarootname, int nscen);

// using the "C"-callbacks interface sInterfaceCallbacks
int solve_usingCallbacks(const string& datarootname, int nscen);



int main(int argc, char ** argv) {
  MPI_Init(&argc, &argv);
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);

  if(argc<3) {
    if (mype == 0) printf("Usage: %s [rawdump root name] [num scenarios] [solution output root name]\n",argv[0]);
    return 1;
  }
  
  string datarootname(argv[1]);
  int nscen = atoi(argv[2]);

  //solve_usingCallbacks(datarootname, nscen);
  solve(datarootname, nscen);


  MPI_Finalize();
  return 0;
}


int solve(const string& datarootname, int nscen) {
  
  rawInput* s = new rawInput(datarootname,nscen);
  PIPSIpmInterface<sFactoryAug, MehrotraStochSolver> pipsIpm(*s);




  delete s;



  return 0;
}


int solve_usingCallbacks(const string& datarootname, int nscen) {

  rawInput* s = new rawInput(datarootname,nscen);
  sInterfaceCallbacks pipsIpm(*s);

  pipsIpm.go();

  // when using callbacks the input can be deleted only after optimization is finished.
  delete s;

  return 0;
}
