#ifndef PAR_MACRO_H_
#define PAR_MACRO_H_

#include <string>
#include <sstream>
#include <iostream>
#include "./Core/Utilities/PerfMetrics.h"

extern int gmyid;
extern int gnprocs;
extern PerfMetrics gprof;

#ifdef VERBOSE
#define ENABLE_VERBOSE 1
#else
#define ENABLE_VERBOSE 0
#endif

#define IF_VERBOSE_DO(X) do { \
	if (ENABLE_VERBOSE) { X } \
} while (0)

#define MESSAGE(x) do { \
  if (ENABLE_VERBOSE) { std::cout<<"["<<gmyid<<"/"<<gnprocs<<"] "<< x << std::endl; } \
} while (0)


template <class T>
void print_array(const std::string& msg, T* data, size_t len)
{
	std::ostringstream oss;
	for(size_t i=0;i<len;i++){
		oss<<data[i]<<", ";
	}
	MESSAGE(""<<msg<<" - "<<"Array [ "<<oss.str()<<"  ]");
}

#define PRINT_ARRAY(M, DATA, LEN) do { \
  if (ENABLE_VERBOSE) { 	std::ostringstream oss; 	\
				for(size_t i=0;i<LEN;i++){ 	\
					oss<<DATA[i]<<", ";    	\
				} \
				std::cout<<"["<<gmyid<<"/"<<gnprocs<<"] "<< M << "Array [ " <<oss.str() <<" ]"<< std::endl; \
			} \
} while (0)


extern void convert_to_csr(int m, int n, int* rowidx, int* colptr, double* elts, int nz, double* ret);


#endif
