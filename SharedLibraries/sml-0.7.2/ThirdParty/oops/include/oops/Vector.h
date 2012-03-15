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

#ifndef VECTOR_H
#define VECTOR_H

class Vector;
typedef enum {notcomputed_proc=-3, computed_proc=-1, Contrib} CompStatus;
typedef enum {Exact_st=-1,  NoStat=0, Contr_st=1} ContStatus; 

#include "oops/Tree.h" 
#include "oops/LogiVector.h" 
#include "oops/DenseVector.h"


#define OuterProdVector(outpd, NonzCC,CompCol,c)  OuterProdDenseVector (GetDenseVectorFromVector (outpd), NonzCC, CompCol, c)

#define NbOfSubVect(v)  ((v->node)->nb_sons)
#define SetEjVector(ej, j) SetCompValueVector(ej, (j), 1.);
#define SetMinusEjVector(ej, j)  SetCompValueVector (ej, j, -1.)
#define GetDenseVectorFromVector(v)  v->subvectors[(v->node)->index]->dense
#define FuncToDenseVector(x,y,f) FuncToDenseVectorStatic(x, y->elts- x->node->begin, f)
#define FuncToDenseBcstVector(x,y,f) FuncToDenseBcstVector_(x, y->elts- x->node->begin, f)
#define SubVector(v,i) v->subvectors[v->node->sons[i]->index]
#define MEMORY_ABOVE(Node)          (Node->above != NULL)
#define MEMORY_NOT_ABOVE(Node)      (Node->above == NULL)
#define MEMORY_LOCAL(Node)          (Node->local == true)
#define MEMORY_NOT_LOCAL(Node)      (Node->local == false)
#define FATHER_VECTOR(Vector, Node) (Vector->subvectors[Node->above->index])
#define SUBVECTOR(Vector, Node)     (Vector->subvectors[Node->index])


class Vector
{
 public:
  typedef double (*KrnFct) (double);
  typedef void (*VectorCallBackFunction) (Vector *);
  typedef void (*VectorCallBackFunction2) (Vector *, Vector *);

  Tree *node;
  Vector **subvectors;
  DenseVector *dense;
  bool dense_alloc;
  bool is_root;
  const char *name;
    
#ifdef WITH_MPI
  ContStatus       sum_stat;
#endif

  Vector();
  Vector(Tree *T, const char *name);
  Vector(Tree *T, const char *name, double *dense);
  ~Vector ();

  Vector *
    getPrimalDualParts(Tree *Atree, int *map);

  void setZero();
  void setDoubleValue(const double a);
  void setComponent(const int j, const double value);
  double getComponent(const int j);
  void markZero ();
  void unmarkZero ();
  void checkZeroMark();
  void copy (Vector *targetV);
  void addVector(const double a, Vector *y);
  void setAypbz(const double a, Vector *y, const double b, Vector *z);
  void invertReg(const double reg);
  void absValue ();
  void addScalar(const double a);
  void multScalar(const double a);
  void maxInPlace (Vector *y);
  void minInPlace (Vector *y);
  int countNz();
  void setRandom();
  void multiplyCompwise(Vector *y, KrnFct f);
  void divideCompwise (Vector *y, KrnFct f);
  void addKroneckProd(const double a, Vector *y, Vector *z);
  void addInvKroneckProd(const double a, Vector *y, Vector *z);
  void invert();
  void applyFunc (KrnFct f);
  void setTheta (Vector *y, Vector *z);
  void print ();
  void printFile (FILE *out);
  void printMatlab(FILE *out, const char *name);

  double ddotPar (Vector *y);
  double infNorm(int *j);
  double minInfNorm(int *j, const double bnd);
  double maxComp(int *j);
  double minComp(int *j);
  void copyFromDense (DenseVector *);
  void copyToDense(DenseVector *);
  void compareToDenseVector(DenseVector *y);
  void setExact();

#ifdef WITH_MPI
  void copyToDenseBcastVector (DenseVector *y);
  void reduceVector (MPI_Op op, MPI_Comm comm);
#endif

  static double ddotLinCombPar (Vector *x, Vector *y, Vector *dx, Vector *dy, 
				double a, double b);

  static void GondzioCorr(Vector *x, Vector *dX, Vector *z, Vector *dZ,
		  Vector *r_xz, Vector *rhs_x,
		  const double barr, const double a1p, const double a1d);

  static void
    GondzioCorrTarget(Vector *x, Vector *dX, Vector *z, Vector *dZ,
		  Vector *r_xz, Vector *rhs_x, Vector *target,
		  const double barr, const double a1p, const double a1d);

  static void
    GondzioCorrVector(Vector *x, Vector *dX, Vector *z, Vector *dZ,
		  Vector *r_xz, Vector *rhs_x, Vector *target,
		  const double barr, const double a1p, const double a1d); 

  static void 
    GondzioCorrOnlySmall(Vector *x, Vector *dX, 
				Vector *z, Vector *dZ,
			       Vector *r_xz, Vector *rhs_x, const double barr,
			   const double a1p, const double a1d);

  static void
    BalanceKronProd(Vector *x, Vector *y, const double barr);

  bool comp_node(Vector *y);

  void fillCallBack(VectorCallBackFunction f);

  static void 
    CallBack2(Vector *x, Vector *y, VectorCallBackFunction2 f);

  void copyToWhere(Vector *y, LogiVector *l);
  void addVectorWhere(const double a, Vector *y, LogiVector *l);
  void setAypbzWhere(const double a, Vector *y, const double b, Vector *z,
		     LogiVector *l);
  void multiplyCompwiseWhere(Vector *y, KrnFct f, LogiVector *l);
  void divideCompwiseWhere(Vector *y, KrnFct f, LogiVector *l);
  void addKroneckProdWhere(double a, Vector *y, Vector *z, LogiVector *l);
  void addInvKroneckProdWhere(const double a, Vector *y, Vector *z, LogiVector *l);
  void setValueWhere(const double a, LogiVector *l);
  void maxInPlaceWhere (Vector *y, LogiVector *l);
  void minInPlaceWhere (Vector *y, LogiVector *l);
  double ddotParWhere (Vector *y, LogiVector *l);
  double infNormWhere (int *j, LogiVector *l);
  double minInfNormWhere (int *j, LogiVector *l);
  double minCompWhere(int *j, LogiVector *l);
  double maxCompWhere(int *j, LogiVector *l);
  void invertWhere (LogiVector *l);
  void getThetaWhere (Vector *y, Vector *z, LogiVector *l);
  void applyFuncWhere (KrnFct f, LogiVector *l);

  static double 
    ddotLinCombParWhere(Vector *x, Vector *y, Vector *dx, Vector *dy,
			const double a, const double b, LogiVector *l);

  static void
    GondzioCorrWhere(Vector *x, Vector *dX, Vector *z, Vector *dZ,
		     Vector *r_xz, Vector *rhs_x, const double barr,
		     const double a1p, const double a1d,
		     LogiVector *l);

  static void
    GondzioCorrOnlySmallWhere(Vector *x, Vector *dX, Vector *z, Vector *dZ,
			      Vector *r_xz, Vector *rhs_x, const double barr,
			      const double a1p, const double a1d,
			      LogiVector *l);

  static void
    BalanceKronProdWhere(Vector *x, Vector *y,
			 const double barr, LogiVector *l);

};

#endif
