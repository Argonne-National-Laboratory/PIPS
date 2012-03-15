/* This file is part of OOPS.
 *
 * OOPS is (c) 2003-2009 Jacek Gondzio and Andreas Grothey, 
 *                       University of Edinburgh
 *
 * OOPS is distributed in a restricted form in the hope that it will be a useful
 * example of what can be done with SML, however it is NOT released under a free
 * software license.
 *
 * You may only redistribute this version of OOPS with a version of SML. You
 * may not link OOPS with code which is not part of SML licensed under the
 * LGPL v3.
 *
 * You may NOT modify, disassemble, or otherwise reverse engineer OOPS.
 *
 * OOPS is distributed WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef LOGIVECTOR_H
#define LOGIVECTOR_H

class LogiVector;

#include "oops/Vector.h" 

typedef double (*KrnFct) (double);

class DenseLogiVector {

 public:
  int dim;
  bool mem_alloc;
  bool is_zero;
  bool *elts;
  char *name;
  int newtag;

  DenseLogiVector(const int vdim, const char *vname);
  DenseLogiVector(DenseLogiVector *V, const int begin, const int vdim);
  ~DenseLogiVector ();

};

class LogiVector
{
 public:
  Tree *node;
  LogiVector **subvectors;
  DenseLogiVector *dense;
  bool is_root;

#ifdef WITH_MPI
  CompStatus proc;
  ContStatus sum_stat;
#endif

  LogiVector();
  LogiVector(Tree *T, const char *name);
  ~LogiVector ();

  void setWhereSmaller(Vector *v, const double bound);
  void setWhereLarger(Vector *v, const double bound);
  void setWhereLarger(Vector *v1, Vector *v2);
  int count();
  void copyTo(LogiVector *l2);
  void setToAnd(LogiVector *l2, LogiVector *l3);
  void setToOr(LogiVector *l2, LogiVector *l3);
  void negate();
  void print();
};


#endif
