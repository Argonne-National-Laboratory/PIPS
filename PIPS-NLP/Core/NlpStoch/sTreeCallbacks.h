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

  StochVector*      createCeqBody() const ;
  StochVector*      createCineqBody()  const ;
  StochVector*      createBarrGrad() const ;    

  int nx() const;
  int my() const; 
  int mz() const; 
  int id() const; 

  void computeGlobalSizes();
 public:
  int NNZA,NNZQ,NNZB,NNZC,NNZD; //global nnz
  void loadLocalSizes();
 protected:
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
