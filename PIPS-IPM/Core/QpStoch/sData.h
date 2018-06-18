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
	 OoqpVector * cupp, OoqpVector * ciupp, long long mcupp,
	 bool exploit2Links = false);

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
  virtual double objectiveValue( QpGenVars * vars );
  virtual void createScaleFromQ();
  virtual void datainput() {};
  virtual bool with2Links() {return use2Links;};

  virtual ~sData();

 protected:
  std::vector<int> linkRowsPermutationA;
  std::vector<int> linkRowsPermutationC;

  void createChildren();
  void destroyChildren();

 private:
  const static double min2LinksRatio = 0.5;
  static std::vector<unsigned int> getAscending2LinkPermutation(const std::vector<int>& linkStartBlocks, size_t nBlocks);

  bool use2Links;
  std::vector<int> linkStartBlocksA;
  std::vector<int> linkStartBlocksC;

  void init2LinksData(bool exploit2links);
  void permuteLinkingRows(const std::vector<unsigned int>& permvecA, const std::vector<unsigned int>& permvecC);
};


#endif
