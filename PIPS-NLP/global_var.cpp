
#include "par_macro.h"
#include "CoinPackedMatrix.hpp"
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>

int gmyid;
int gnprocs;

void print_array(const std::string& msg, double* data, size_t len)
{
	std::ostringstream oss;
	for(size_t i=0;i<len;i++){
		oss<<data[i]<<", ";
	}
	PAR_DEBUG(""<<msg<<" - "<<"Array [ "<<oss.str()<<"  ]");
}
void print_array(const std::string& msg, int* data, size_t len)
{
	std::ostringstream oss;
	for(size_t i=0;i<len;i++){
		oss<<data[i]<<", ";
	}
	PAR_DEBUG(""<<msg<<" - "<<"Array [ "<<oss.str()<<"  ]");
}

void convert_to_csr(int m, int n, int* rowidx, int* colptr, double* elts, int nz, double* ret)
{
	if(nz!=0){
		print_array("rowidx",rowidx,nz);
		print_array("colptr",colptr,n+1);
		print_array("elts",elts,nz);
		CoinPackedMatrix mat;
		mat.copyOf(true,m,n,nz,&elts[0],&rowidx[0],&colptr[0],0);
		assert(nz == mat.getNumElements());
//		mat.dumpMatrix();
		mat.reverseOrdering();
//		mat.dumpMatrix();
		const double* csr_elts = mat.getElements();
		print_array("csr_elts",csr_elts,nz);
		for(int i=0;i<nz;i++){
			ret[i] = csr_elts[i];
		}
	}
}
