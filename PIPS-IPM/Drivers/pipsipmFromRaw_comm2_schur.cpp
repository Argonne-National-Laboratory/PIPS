/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */
#include <stdio.h>
#include <stdlib.h>

#include "rawInput.hpp"
#include "PIPSIpmInterface.h"

#include "sFactoryAugComm2SchurLeaf.h"
#include "MehrotraStochSolver.h"


#include <string>
#include <sstream>

int main(int argc, char ** argv) {
  MPI_Init(&argc, &argv);
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);

  if(argc<3) {
    if (mype == 0) printf("Usage: %s [rawdump root name] [num scenarios] [outer solve (opt)] [inner solve (opt)]\n",argv[0]);
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
  rawInput* s = new rawInput(datarootname,nscen);
  if(mype==0) cout <<  " raw input created from " << datarootname<< endl;
  PIPSIpmInterface<sFactoryAugComm2SchurLeaf, MehrotraStochSolver> pipsIpm(*s);

  pips_options::setIntParameter("OUTER_SOLVE", outerSolve);
  pips_options::setIntParameter("INNER_SC_SOLVE", innerSolve);

  if(mype==0) cout <<  "PIPSIpmInterface created" << endl;
  delete s;
  if(mype==0) cout <<  "rawInput deleted ... starting to solve" << endl;

  pipsIpm.go();

  MPI_Finalize();
  return 0;
}

