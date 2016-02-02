

extern int gmyid;
extern int gnprocs;

#define PAR_DEBUG(x) do { \
  if (false) { std::cout<<"["<<gmyid<<"/"<<gnprocs<<"] "<< x << std::endl; } \
} while (0)
