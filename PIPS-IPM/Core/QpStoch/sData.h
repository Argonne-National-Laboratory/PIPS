#ifndef DATAQPSTOCH
#define DATAQPSTOCH

#include "QpGenData.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "DoubleMatrixHandle.h"

#include <vector>

class sTree;
class LinearAlgebraPackage;
class QpGenVars;

class sData : public QpGenData {
 public:
  /** constructor that makes data objects of the specified dimensions */
  sData( sTree* tree);

  /** constructor that sets up pointers to the data objects that are
      passed as arguments */
  sData( sTree* stochNode,
    OoqpVector * c, SymMatrix * Q,
    OoqpVector * xlow, OoqpVector * ixlow, long long nxlow,
    OoqpVector * xupp, OoqpVector * ixupp, long long nxupp,
    GenMatrix * A, OoqpVector * bA,
    GenMatrix * C,
    OoqpVector * clow, OoqpVector * iclow, long long mclow,
    OoqpVector * cupp, OoqpVector * ciupp, long long mcupp);


  std::vector<sData*> children;
  void AddChild(sData* child);
  sTree* stochNode;
 public:
  long long nxlow, nxupp, mclow, mcupp;

  int getLocalnx();
  int getLocalmy();
  int getLocalmyl();
  int getLocalmz();
  int getLocalmzl();
  int getLocalSizes(int& nx, int& my, int& mz);
  int getLocalSizes(int& nx, int& my, int& mz, int& myl, int& mzl);

  int getLocalNnz(int& nnzQ, int& nnzB, int& nnzD);

  SparseSymMatrix& getLocalQ();
  SparseGenMatrix& getLocalCrossHessian();
  SparseGenMatrix& getLocalA();
  SparseGenMatrix& getLocalB();
  SparseGenMatrix& getLocalF();
  SparseGenMatrix& getLocalC();
  SparseGenMatrix& getLocalD();
  SparseGenMatrix& getLocalG();

  void sync();
 public:
  virtual void writeToStreamDense(ostream& out) const;
  virtual sData* cloneFull() const;
  virtual double objectiveValue( QpGenVars * vars );
  virtual void createScaleFromQ();
  virtual void datainput() {};

  virtual ~sData();

 protected:
  void createChildren();
  void destroyChildren();

};


#endif
