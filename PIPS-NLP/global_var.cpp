
#include "par_macro.h"
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
