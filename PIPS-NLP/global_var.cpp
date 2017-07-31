
#include "global_var.h"
#include "CoinPackedMatrix.hpp"
#include <sstream>
#include <iostream>
#include <fstream>
#include <mpi.h>

int gmyid;
int gnprocs;
int giterNum;
int gmyid_node;
int gnprocs_node;
MPI_Comm comm_node;
double *gwindow=NULL;
int *gipiv=NULL;
MPI_Win gwin;
MPI_Win gwin_ipiv;
#ifdef NLPTIMING
PerfMetrics gprof = PerfMetrics::getPerfMetrics();
#endif

void convert_to_csr(int m, int n, int* rowidx, int* colptr, double* elts, int nz, double* ret)
{
  if(nz!=0){
    PRINT_ARRAY("rowidx",rowidx,nz);
    PRINT_ARRAY("colptr",colptr,n+1);
    PRINT_ARRAY("elts",elts,nz);
    CoinPackedMatrix mat;
    mat.copyOf(true,m,n,nz,&elts[0],&rowidx[0],&colptr[0],0);
    assert(nz == mat.getNumElements());
    //mat.dumpMatrix();
    mat.reverseOrdering();
    //mat.dumpMatrix();
    const double* csr_elts = mat.getElements();
    PRINT_ARRAY("csr_elts",csr_elts,nz);
    for(int i=0;i<nz;i++){
      ret[i] = csr_elts[i];
    }
  }
}

