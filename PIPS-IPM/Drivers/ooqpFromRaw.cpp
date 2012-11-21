#include "rawInput.hpp"
#include "OOQPInterface.hpp"
#include "QpGenSparseMa57.h"
#include "QpGenSparseMa27.h"
#include "MehrotraSolver.h"

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
  
  //OOQPInterface<MehrotraSolver,QpGenSparseMa57> ooqp(*s);
  OOQPInterface<MehrotraSolver,QpGenSparseMa27> ooqp(*s);

  delete s;

  ooqp.go();
  

  return 0;
}
