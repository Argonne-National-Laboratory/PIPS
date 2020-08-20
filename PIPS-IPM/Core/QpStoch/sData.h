#ifndef DATAQPSTOCH
#define DATAQPSTOCH

#include "QpGenData.h"
#include "sResiduals.h"
#include "sVars.h"
#include "StochSymMatrix.h"
#include "SparseSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "DoubleMatrixHandle.h"
#include "pipschecks.h"
#include "pipsport.h"
#include <vector>

class sTree;
class LinearAlgebraPackage;

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

  
private: 
  // returns permutation vector or empty vector if no permutation has been performed
  std::vector<unsigned int> getLinkVarsPerm() const;
  // returns inverse permutation vector or empty vector if no permutation has been performed
  std::vector<unsigned int> getLinkVarsPermInv() const;
  // returns inverse permutation vector or empty vector if no permutation has been performed
  std::vector<unsigned int> getLinkConsEqPermInv () const;
  // returns inverse permutation vector or empty vector if no permutation has been performed
  std::vector<unsigned int> getLinkConsIneqPermInv() const;

public:
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

  // distributed version
  int getSchurCompMaxNnzDist(int blocksStart, int blocksEnd);
  bool exploitingLinkStructure() {return useLinkStructure;};

  SparseSymMatrix* createSchurCompSymbSparseUpper();

  // distributed version
  SparseSymMatrix* createSchurCompSymbSparseUpperDist(int blocksStart, int blocksEnd);

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
  sResiduals* getResidsUnperm(const sResiduals& resids, const sData& unpermData) const;
  sVars* getVarsUnperm(const sVars& vars, const sData& unpermData) const;

  void sync();
  bool isRootNodeInSync() const;

 public:
  virtual void writeToStreamDense(ostream& out) const;
  void writeMPSformat(ostream& out);
  void writeMPSColumns(ostream& out);
  virtual sData* cloneFull(bool switchToDynamicStorage = false) const;
  double objectiveValue( const QpGenVars * vars ) const override;
  virtual void createScaleFromQ();

  void cleanUpPresolvedData(const StochVectorBase<int>& rowNnzVecA, const StochVectorBase<int>& rowNnzVecC, const StochVectorBase<int>& colNnzVec);

  // marker that indicates whether a Schur complement row is (2-link) local
  const std::vector<bool>& getSCrowMarkerLocal() const;

  // marker that indicates whether a Schur complement row is (2-link) local and owned by this MPI process
  const std::vector<bool>& getSCrowMarkerMyLocal() const;

  // number of sparse 2-link equality rows
  int n2linkRowsEq() const;

  // number of sparse 2-link inequality rows
  int n2linkRowsIneq() const;

  // start and end positions for local 2-links in Schur complement that are non-zero if only
  // blocks greater equal blocksStart and smaller blocksEnd are considered
  void getSCrangeMarkers(int blocksStart, int blocksEnd, int& local2linksStartEq, int& local2linksEndEq,
        int& local2linksStartIneq, int& local2linksEndIneq);

  // start and end positions for local 2-links in Schur complement that are owned by
  // blocks greater equal blocksStart and smaller blocksEnd are considered
  void getSCrangeMarkersMy(int blocksStart, int blocksEnd, int& local2linksStartEq, int& local2linksEndEq,
        int& local2linksStartIneq, int& local2linksEndIneq);

  virtual ~sData();

 protected:
  void createChildren();
  void destroyChildren();

 private:
  int n0LinkVars;
  constexpr static int nLinkStats = 6;
  constexpr static double minStructuredLinksRatio = 0.5;
  static std::vector<unsigned int> get0VarsRightPermutation(const std::vector<int>& linkVarsNnzCount);
  static std::vector<unsigned int> getAscending2LinkPermutation(std::vector<int>& linkStartBlocks, size_t nBlocks);

  // returns number of block rows
  static int getSCdiagBlocksNRows(const std::vector<int>& linkStartBlockLengths);

  // returns number of non-zero block rows within specified range
  static int getSCdiagBlocksNRows(const std::vector<int>& linkStartBlockLengths,
        int blocksStart, int blocksEnd);

  // returns number of owned block rows within specified range
  static int getSCdiagBlocksNRowsMy(const std::vector<int>& linkStartBlockLengths,
        int blocksStart, int blocksEnd);

  // max nnz in Schur complement diagonal block signified by given vector
  static int getSCdiagBlocksMaxNnz(size_t nRows,
        const std::vector<int>& linkStartBlockLengths);

  // distributed version
  static int getSCdiagBlocksMaxNnzDist(size_t nRows,
        const std::vector<int>& linkStartBlockLengths, int blocksStart, int blocksEnd);

  // max nnz in Schur complement mixed block signified by given vectors
  static int  getSCmixedBlocksMaxNnz(size_t nRows, size_t nCols,
        const std::vector<int>& linkStartBlockLength_Left,
        const std::vector<int>& linkStartBlockLength_Right);

  // distributed version
  static int  getSCmixedBlocksMaxNnzDist(size_t nRows, size_t nCols,
        const std::vector<int>& linkStartBlockLength_Left,
        const std::vector<int>& linkStartBlockLength_Right, int blocksStart, int blocksEnd);

  // number of sparse 2-link rows
  static int n2linksRows(const std::vector<int>& linkStartBlockLengths);

  static std::vector<int> get2LinkLengthsVec(const std::vector<int>& linkStartBlocks, size_t nBlocks);

  /* a two link must be in two blocks directly after one another */
  bool useLinkStructure;
  /* number of entries in each linking column */
  std::vector<int> linkVarsNnz;
  /* which blocks do the individual two-links start in */
  std::vector<int> linkStartBlockIdA;
  std::vector<int> linkStartBlockIdC;
  /* how many two-links start in block i */
  std::vector<int> linkStartBlockLengthsA;
  std::vector<int> linkStartBlockLengthsC;
  std::vector<unsigned int> linkVarsPermutation;
  std::vector<unsigned int> linkConsPermutationA;
  std::vector<unsigned int> linkConsPermutationC;
  std::vector<bool> isSCrowLocal;
  std::vector<bool> isSCrowMyLocal;

  void initDistMarker(int blocksStart, int blocksEnd);


  void permuteLinkingVars();
  void permuteLinkingCons();
};


#endif
