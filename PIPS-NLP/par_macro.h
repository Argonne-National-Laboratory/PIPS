#include <string>

extern int gmyid;
extern int gnprocs;

#define PAR_DEBUG(x) do { \
  if (false) { std::cout<<"["<<gmyid<<"/"<<gnprocs<<"] "<< x << std::endl; } \
} while (0)

extern void print_array(const std::string& msg, double* data, size_t len);
extern void print_array(const std::string& msg, int* data, size_t len);
void convert_to_csr(int m, int n, int* rowidx, int* colptr, double* elts, int nz, double* ret);
