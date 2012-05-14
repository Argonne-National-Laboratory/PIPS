#include "rawInput.hpp"
#include "PIPSIPMInterface.h"

using namespace std;

int main(int argc, char ** argv) {

  MPI_Init(&argc, &argv);
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);

  if(argc<3) {
    if (mype == 0) printf("Usage: %s [rawdump root name] [num scenarios] [solution output root name]\n",argv[0]);
    return 1;
  }
  
  string datarootname(argv[1]);
  int nscen = atoi(argv[2]);

  rawInput* s = new rawInput(datarootname,nscen);
  
  PIPSIPMInterface pipsIpm(*s);
  pipsIpm.loadData();
  pipsIpm.go();

  delete s;

  return 0;
}
