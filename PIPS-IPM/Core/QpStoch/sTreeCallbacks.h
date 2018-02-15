/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef STOCH_TREE_CALLBACKS
#define STOCH_TREE_CALLBACKS

#include "sTree.h"
#include "StochInputTree.h"

/** This class creates objects when  the problem is specified by C callbacks.
 *  Obsolete and present only to ensure compatibility with older versions of the code.
 *  The new sTree implementation, C++-like is sTreeImpl.
 */

class sTreeCallbacks : public sTree
{
 public:
  sTreeCallbacks(StochInputTree* root);
  sTreeCallbacks(const std::vector<StochInputTree::StochInputNode*> &localscens);
  sTreeCallbacks(StochInputTree::StochInputNode* data_);
  ~sTreeCallbacks();

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

  int nx() const;
  int my() const;
  int myl() const;
  int mz() const; 
  int mzl() const;
  int id() const; 

  void computeGlobalSizes();
 public:
  int NNZA,NNZQ,NNZB,NNZBl,NNZC,NNZD,NNZDl; //global nnz
  int NNZA_INACTIVE,NNZQ_INACTIVE,NNZB_INACTIVE,NNZBl_INACTIVE,NNZC_INACTIVE,NNZD_INACTIVE,NNZDl_INACTIVE; //global inactive nnz
  long long N_INACTIVE,MY_INACTIVE,MZ_INACTIVE; //global inactive sizes
  int nx_active, my_active, mz_active, myl_active, mzl_active;
  int nx_inactive, my_inactive, mz_inactive, myl_inactive, mzl_inactive;

  void loadLocalSizes();

  virtual void switchToPresolvedData();
  virtual void switchToOriginalData();
  virtual bool isPresolved();
  virtual bool hasPresolved();
  virtual void initPresolvedData(const StochSymMatrix& Q, const StochGenMatrix& A, const StochGenMatrix& C, const StochVector& nxVec, const StochVector& myVec, const StochVector& mzVec);
 protected:
  bool isDataPresolved;
  bool hasPresolvedData;

  sTreeCallbacks();
  StochInputTree::StochInputNode* data; //input data
  // in POOLSCEN case, only root node has non-null data
  StochInputTree* tree;
  std::vector<StochInputTree::StochInputNode*> scens;
  StochInputTree::StochInputNode* fakedata; //convenient struct for holding n,my,mz etc
  // holds stoch trees for each of the scenarios that are combined at this node
  // this is just a convenience to reuse the create* and newVector* functions
  std::vector<sTreeCallbacks*> real_children;
};


#endif
