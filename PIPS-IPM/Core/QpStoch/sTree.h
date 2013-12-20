/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef STOCH_TREE_BASE
#define STOCH_TREE_BASE

#include "StochResourcesMonitor.h"
#include "StochVector.h"
#include "StochGenMatrix.h"
#include "StochSymMatrix.h"

#include "list"

#include "mpi.h"

class sTree
{
 public:
  virtual ~sTree();

  int NumberOfChildren() const { return children.size(); }

  virtual void computeGlobalSizes() = 0;
  void GetGlobalSizes(long long& NXOut, long long& MYOut, long long& MZOut);
  //void GetLocalSizes(int& nxOut, int& myOut, int& mzOut);

  void assignProcesses  ( MPI_Comm comm = MPI_COMM_WORLD);
  void assignProcesses  ( MPI_Comm, vector<int>&);

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
  void syncStochVector(StochVector& vec);

  void syncStochGenMatrix(StochGenMatrix& mat);
  void syncStochSymMatrix(StochSymMatrix& mat);

  virtual StochSymMatrix*   createQ() const = 0;
  virtual StochVector*      createc() const = 0;

  virtual StochVector*      createxlow()  const = 0;
  virtual StochVector*      createixlow() const = 0;
  virtual StochVector*      createxupp()  const = 0;
  virtual StochVector*      createixupp() const = 0;


  virtual StochGenMatrix*   createA() const = 0;
  virtual StochVector*      createb() const = 0;


  virtual StochGenMatrix*   createC() const = 0;
  virtual StochVector*      createclow()  const = 0;
  virtual StochVector*      createiclow() const = 0;
  virtual StochVector*      createcupp()  const = 0;
  virtual StochVector*      createicupp() const = 0;

  StochVector*      newPrimalVector() const;
  StochVector*      newDualYVector()  const;
  StochVector*      newDualZVector()  const;
  StochVector*      newPrimalVectorEmpty() const;
  StochVector*      newDualYVectorEmpty()  const;
  StochVector*      newDualZVectorEmpty()  const;

  StochVector*      newRhs();

  int innerSize(int which);
  virtual int nx() const = 0;
  virtual int my() const = 0; 
  virtual int mz() const = 0; 
  virtual int id() const = 0; 

  //returns the global load, i.e. statistic based on the NNZs and
  //dimensions of the node (and subnodes) subproblem before any iteration or the CPU
  //time of this node and its subnodes  after the first iteration.
  double processLoad() const;

 protected:
  sTree();

  void   toMonitorsList(list<NodeExecEntry>&);
  void fromMonitorsList(list<NodeExecEntry>&);

  void computeNodeTotal();

  void saveCurrentCPUState();

  int isInVector(int elem, const vector<int>& vec);


 public:
  long long N,MY,MZ; //global sizes
  int np; //n for the parent

  double IPMIterExecTIME;
  std::vector<sTree*> children;

  vector<int> idx_EqIneq_Map;

  static int numProcs;

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
//to be called after assignProcesses
  virtual void loadLocalSizes()=0;
};

#endif 
