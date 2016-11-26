/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#include "rawInput.hpp"
#include "sInterfaceCallbacks.h"
#include "PIPSIpmInterface.h"

#include "sFactoryAug.h"
#include "MehrotraStochSolver.h"

#include <limits>

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
  if (mype == 0) cout << argv[0] << " starting..." << endl;  
  string datarootname(argv[1]);
  int nscen = atoi(argv[2]);

  //solve_usingCallbacks(datarootname, nscen);
  solve(datarootname, nscen);


  MPI_Finalize();
  return 0;
}


int solve(const string& datarootname, int nscen) {
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  int nprocs; MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if(0==mype) cout << "Using a total of " << nprocs << " MPI processes." << endl;

  rawInput* s = new rawInput(datarootname,nscen,MPI_COMM_WORLD);
  if (mype == 0) cout << "rawInput created .." << endl;
  PIPSIpmInterface<sFactoryAug, MehrotraStochSolver> pipsIpm(*s);
  if (mype == 0) cout << "PIPSIpmInterface created .." << endl;
  delete s;
  if (mype == 0) cout << "rawInput deleted ... solving" << endl;

  pipsIpm.go();

  double obj = pipsIpm.getObjective();
  //cout << "PIPS-IPM: optimal objective: "  << obj << endl;
  if (mype == 0) printf("PIPS-IPM: optimal objective: %.8f \n", obj);

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
