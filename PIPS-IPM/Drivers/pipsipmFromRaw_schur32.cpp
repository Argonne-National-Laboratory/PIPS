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

  if(mype==0) cout << argv[0] << " starting ..." << endl;
  
  string datarootname(argv[1]);
  int nscen = atoi(argv[2]);

  rawInput* s = new rawInput(datarootname,nscen);
  if(mype==0) cout <<  " raw input created ..." << endl;
  PIPSIpmInterface<sFactoryAugSchur32Leaf, MehrotraStochSolver> pipsIpm(*s);
  if(mype==0) cout <<  "PIPSIpmInterface created (32-bit second stage factorization) " << endl;
  delete s;
  if(mype==0) cout <<  "rawInput deleted ... starting to solve" << endl;

  pipsIpm.go();

  MPI_Finalize();
  return 0;
}

