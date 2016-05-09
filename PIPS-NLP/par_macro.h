#ifndef PAR_MACRO_H_
#define PAR_MACRO_H_

#include <string>
#include <sstream>
#include <iostream>

extern int gmyid;
extern int gnprocs;

#define DEBUG(X) do { \
	if (0) { X } \
} while (0)

#define PAR_DEBUG(x) do { \
  if (0) { std::cout<<"["<<gmyid<<"/"<<gnprocs<<"] "<< x << std::endl; } \
} while (0)


template <class T>
void print_array(const std::string& msg, T* data, size_t len)
{
	std::ostringstream oss;
	for(size_t i=0;i<len;i++){
		oss<<data[i]<<", ";
	}
	PAR_DEBUG(""<<msg<<" - "<<"Array [ "<<oss.str()<<"  ]");
}

#define PRINT_ARRAY(M, DATA, LEN) do { \
  if (0) { 	std::ostringstream oss; 	\
				for(size_t i=0;i<LEN;i++){ 	\
					oss<<DATA[i]<<", ";    	\
				} \
				std::cout<<"["<<gmyid<<"/"<<gnprocs<<"] "<< M << "Array [ " <<oss.str() <<" ]"<< std::endl; \
			} \
} while (0)


extern void convert_to_csr(int m, int n, int* rowidx, int* colptr, double* elts, int nz, double* ret);

#endif
