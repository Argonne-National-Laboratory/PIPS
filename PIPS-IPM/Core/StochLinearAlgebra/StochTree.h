#ifndef STOCH_TREE
#define STOCH_TREE

#include "StochInputTree.h"
#include "StochResourcesMonitor.h"

#include "mpi.h"

#include <vector>
#include <list>
class Data;
class QpGen;
class QpGenStoch;
class sData;
class StochSymMatrix;
class StochGenMatrix;
class StochVector;
class QpGenStochLinsys;

//#define POOLSCEN 1

class StochTree {
 public:

  StochTree(StochInputTree* root);
  StochTree(const std::vector<StochInputTree::StochInputNode*> &localscens);
  StochTree(StochInputTree::StochInputNode* data_);
  virtual ~StochTree();

  int NumberOfChildren() const { return children.size(); }

  void computeGlobalSizes();
  void GetGlobalSizes(int& NXOut, int& MYOut, int& MZOut);
  void GetLocalSizes(int& nxOut, int& myOut, int& mzOut);


  virtual void assignProcesses  ( );
  virtual void assignProcesses  (MPI_Comm, vector<int>&);
  //void assignProcesses  (MPI_Comm);

  MPI_Comm commWrkrs, myOldMpiComm; //workers only
  vector<int> myProcs, myOldProcs;

  MPI_Comm commP2ZeroW;   // preconditioner (rank P+1) and special (rank 0) worker
  static int rankPrcnd;   // rank of preconditioner
  static int rankZeroW;   // rank of the "special" worker (the root, or 0-rank process)
  static int rankMe;      // rank of the running process 

  void startMonitors(); void startNodeMonitors();
  void stopMonitors();  void stopNodeMonitors();
  void syncMonitoringData(vector<double>& vCPUTotal);
  bool balanceLoad();
  bool balanceLoadPrecond();

  void getSyncInfo(int myRank, int& syncNeeded, int& sendOrRecv, int& toFromCPU );
  void syncPrimalVector(StochVector& vec);
  void syncDualYVector(StochVector& vec);
  void syncDualZVector(StochVector& vec);
  void syncStochVector_old(StochVector& vec, int whatType);
  void syncStochVector(StochVector& vec);

  void syncStochGenMatrix(StochGenMatrix& mat);
  void syncStochSymMatrix(StochSymMatrix& mat);

  StochSymMatrix*   createQ() const;
  StochVector*      createc() const;

  StochVector*      createxlow()  const;
  StochVector*      createixlow() const;
  StochVector*      createxupp()  const;
  StochVector*      createixupp() const;


  StochGenMatrix*   createA() const;
  StochVector*      createb() const;


  StochGenMatrix*   createC() const;
  StochVector*      createclow()  const;
  StochVector*      createiclow() const;
  StochVector*      createcupp()  const;
  StochVector*      createicupp() const;

  StochVector*      newPrimalVector() const;
  StochVector*      newDualYVector()  const;
  StochVector*      newDualZVector()  const;
  StochVector*      newPrimalVectorEmpty() const;
  StochVector*      newDualYVectorEmpty()  const;
  StochVector*      newDualZVectorEmpty()  const;

  StochVector*      newRhs();

  int innerSize(int which);
  int nx() const;
  int my() const; 
  int mz() const; 
  int id() const; 

  //void* user_data() const { return data->user_data; }
  //FNNZ  fnnzQ() const { return data->fnnzQ; }
  //FNNZ  fnnzA() const { return data->fnnzA; }
  //FNNZ  fnnzC() const { return data->fnnzC; }

  //FMAT  fQ() const { return data->fQ; };
  //FMAT  fA() const { return data->fA; };
  //FMAT  fC() const { return data->fC; };

  //returns the global load, i.e. statistic based on the NNZs and
  //dimensions of the node (and subnodes) subproblem before any iteration or the CPU
  //time of this node and its subnodes  after the first iteration.
  double processLoad() const;

 protected:
  StochTree();

  void   toMonitorsList(list<NodeExecEntry>&);
  void fromMonitorsList(list<NodeExecEntry>&);

  void computeNodeTotal();

  void saveCurrentCPUState();

  int isInVector(int elem, const vector<int>& vec);

 protected:
  StochInputTree::StochInputNode* data; //input data
	// in POOLSCEN case, only root node has non-null data
  StochInputTree* tree;
	std::vector<StochInputTree::StochInputNode*> scens;
	StochInputTree::StochInputNode* fakedata; //convenient struct for holding n,my,mz etc
	// holds stoch trees for each of the scenarios that are combined at this node
	// this is just a convenience to reuse the create* and newVector* functions
	std::vector<StochTree*> real_children;	


	//function pointers are invalid
 public:
  int N,MY,MZ; //global sizes
  int NNZA,NNZQ,NNZB,NNZC,NNZD; //global nnz
  int np; //n for the parent

  double IPMIterExecTIME;
  std::vector<StochTree*> children;
  static int numProcs;
  //std::vector<int> processes;

  StochNodeResourcesMonitor    resMon;
  static StochIterateResourcesMonitor iterMon;
#ifdef STOCH_TESTING
  void displayProcessInfo(int onWhichRank=0);
  void displayProcessInfo(char* tab);
  void runTestNGP();
  void displayExecTimes(int onWhichRank=0);
  void displayExecTimes(char* szTabbing);

  void displayVectorVsTreeStructure(StochVector& stVec, int rank, char* szTab);
  void displayVectorVsTreeStructure(StochVector& stVec, int rank);

  void displayMatVsTreeStructure(StochGenMatrix& stVec, int myRank, char* tab);
  void displayMatVsTreeStructure(StochGenMatrix& stVec, int myRank);

  void displayMatVsTreeStructure(StochSymMatrix& stVec, int myRank, char* tab);
  void displayMatVsTreeStructure(StochSymMatrix& stVec, int myRank);
#endif
};

#endif
