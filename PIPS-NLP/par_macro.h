
#include "mpi.h"

extern int gmyid;
extern int gnprocs;

#define PAR_DEBUG(x) do { \
  if (true) { std::cout<<"["<<gmyid<<"/"<<gnprocs<<"] "<< x << std::endl; } \
} while (0)
