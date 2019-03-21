#ifndef DATAQPSTOCH
#define DATAQPSTOCH

#include "QpGenData.h"
#include "StochSymMatrix.h"
#include "SparseSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "DoubleMatrixHandle.h"
#include "pipschecks.h"
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

  // returns permutation vector or empty vector if no permutation has been performed
  std::vector<unsigned int> getLinkVarsPerm() const;
  // returns inverse permutation vector or empty vector if no permutation has been performed
  std::vector<unsigned int> getLinkVarsPermInv() const;
  // returns inverse permutation vector or empty vector if no permutation has been performed
  std::vector<unsigned int> getLinkConsEqPermInv () const;
  // returns inverse permutation vector or empty vector if no permutation has been performed
  std::vector<unsigned int> getLinkConsIneqPermInv() const;


  int getLocalnx();
  int getLocalmy();
  int getLocalmyl();
  int getLocalmz();
  int getLocalmzl();
  int getLocalSizes(int& nx, int& my, int& mz);
  int getLocalSizes(int& nx, int& my, int& mz, int& myl, int& mzl);

  int getLocalNnz(int& nnzQ, int& nnzB, int& nnzD);
  int getN0LinkVars() {return n0LinkVars;}
  // returns upper bound on number of non-zeroes in Schur complement
  int getSchurCompMaxNnz();
  bool exploitingLinkStructure() {return useLinkStructure;};
  SparseSymMatrix* createSchurCompSymbSparseUpper();

  SparseSymMatrix& getLocalQ();
  SparseGenMatrix& getLocalCrossHessian();
  SparseGenMatrix& getLocalA();
  SparseGenMatrix& getLocalB();
  SparseGenMatrix& getLocalF();
  SparseGenMatrix& getLocalC();
  SparseGenMatrix& getLocalD();
  SparseGenMatrix& getLocalG();

  void printLinkVarsStats();
  void printLinkConsStats();

  void activateLinkStructureExploitation();

  void sync();

 public:
  virtual void writeToStreamDense(ostream& out) const;
  void writeMPSformat(ostream& out);
  void writeMPSColumns(ostream& out);
  virtual sData* cloneFull(bool switchToDynamicStorage = false) const;
  virtual double objectiveValue( QpGenVars * vars );
  virtual void createScaleFromQ();
  virtual void datainput() {};

  void cleanUpPresolvedData(const StochVector& rowNnzVecA, const StochVector& rowNnzVecC, const StochVector& colNnzVec);

  virtual ~sData();

 protected:
  void createChildren();
  void destroyChildren();

 private:
  int n0LinkVars;
  const static int nLinkStats = 6;
  const static double minStructuredLinksRatio = 0.5;
  static std::vector<unsigned int> get0VarsRightPermutation(const std::vector<int>& linkVarsNnzCount);
  static std::vector<unsigned int> getAscending2LinkPermutation(std::vector<int>& linkStartBlocks, size_t nBlocks);

  // max nnz in Schur complement diagonal block signified by given vector
  static int getSCdiagBlocksMaxNnz(size_t nRows,
        const std::vector<int>& linkStartBlockLengths);

  // max nnz in Schur complement mixed block signified by given vectors
  static int  getSCmixedBlocksMaxNnz(size_t nRows, size_t nCols,
        const std::vector<int>& linkStartBlockLength_Left,
        const std::vector<int>& linkStartBlockLength_Right);

  // number of sparse 2-link rows
  static int n2linksRows(const std::vector<int>& linkStartBlockLengths);

  static std::vector<int> get2LinkLengthsVec(const std::vector<int>& linkStartBlocks, size_t nBlocks);

  bool useLinkStructure;
  std::vector<int> linkVarsNnz;
  std::vector<int> linkStartBlockIdA;
  std::vector<int> linkStartBlockIdC;
  std::vector<int> linkStartBlockLengthsA;
  std::vector<int> linkStartBlockLengthsC;
  std::vector<unsigned int> linkVarsPermutation;
  std::vector<unsigned int> linkConsPermutationA;
  std::vector<unsigned int> linkConsPermutationC;

  void permuteLinkingVars();
  void permuteLinkingCons();
};


#endif
