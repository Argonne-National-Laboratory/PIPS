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
	 OoqpVector * xlow, OoqpVector * ixlow, int nxlow,
	 OoqpVector * xupp, OoqpVector * ixupp, int nxupp,
	 GenMatrix * A, OoqpVector * bA,
	 GenMatrix * C,
	 OoqpVector * clow, OoqpVector * iclow, int mclow,
	 OoqpVector * cupp, OoqpVector * ciupp, int mcupp );

  std::vector<sData*> children;
  void AddChild(sData* child);
  sTree* stochNode;
 public:
  int nxlow, nxupp, mclow, mcupp;

  int getLocalnx();
  int getLocalmy();
  int getLocalmz();
  int getLocalSizes(int& nx, int& my, int& mz);

  int getLocalNnz(int& nnzQ, int& nnzB, int& nnzD);

  SparseSymMatrix& getLocalQ();
  SparseGenMatrix& getLocalCrossHessian();
  SparseGenMatrix& getLocalA();
  SparseGenMatrix& getLocalB();
  SparseGenMatrix& getLocalC();
  SparseGenMatrix& getLocalD();

  void sync();
 public:


  virtual double objectiveValue( QpGenVars * vars );
  virtual void createScaleFromQ();
  virtual void datainput() {};

  virtual ~sData();


 protected:
  void createChildren();
  void destroyChildren();

};


#endif
