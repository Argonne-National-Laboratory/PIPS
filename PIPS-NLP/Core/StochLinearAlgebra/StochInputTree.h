/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef STOCH_INPUT_TREE
#define STOCH_INPUT_TREE

#include <vector>


/**
 * The following types define callback functions passed by user to pass 
 * data for each node to the solver.
 */
extern "C" 
typedef int (*FNNZ)(void* user_data, int id, int* nnz);

extern "C" 
typedef int (*FMAT)(void* user_data, int id, 
		    int* krowM, int* jcolM, double* M);
extern "C"
typedef int (*FVEC)(void* user_data, int id, double* vec, int len);

extern "C"
typedef int (*FLEN)(void* user_data, int id, int* len);


class StochInputTree {
  friend class sTreeCallbacks;
 public:

  /** 
   * Inner class that contains the node related data.
   */
  class StochInputNode {
    friend class StochInputTree; friend class StochTree; friend class sTreeCallbacks;
  public:
  StochInputNode();
	StochInputNode(void* user_data, int id, 
	   int n, int my, int mz,
	   FMAT fQ, FNNZ fnnzQ, FVEC fc,  
	   FMAT fA, FNNZ fnnzA, FMAT fB, FNNZ fnnzB, 
	   FVEC fb, 
	   FMAT fC, FNNZ fnnzC, FMAT fD, FNNZ fnnzD,
	   FVEC fclow, FVEC ficlow, FVEC fcupp, FVEC ficupp,
	   FVEC fxlow, FVEC fixlow, FVEC fxupp, FVEC fixupp,
	   bool deleteUserData=false);
  ~StochInputNode();

  protected:
    
    int id; int n,my,mz; 

    int nnzQ, nnzA, nnzB, nnzC, nnzD;
  protected:
    //call back functions
    FNNZ fnnzQ, fnnzA, fnnzB, fnnzC, fnnzD;
    FMAT fQ, fA, fB, fC, fD;

    FVEC fc, fb;
    FVEC fclow, fcupp, ficlow, ficupp;
    FVEC fxlow, fxupp, fixlow, fixupp;
    
    void *user_data;

    bool deleteUserData;

  }; // end of inner class StochInputNode

  friend class StochTree;

  ////////////////////////////
  // CONSTRUCTOR & DESTRUCTOR
  ////////////////////////////
 public:
  StochInputTree(const StochInputNode& root);
  StochInputTree(StochInputNode* root);
  virtual ~StochInputTree();
 protected:
  StochInputTree();
  /////////////////////////
  // METHODS
  /////////////////////////
 public:
  void AddChild(const StochInputNode  &node);
  void AddChild(StochInputTree *subTree);

  /////////////////////////
  // DATA members
  /////////////////////////
 public:
  std::vector<StochInputTree*> children;
 protected:
  StochInputNode* nodeInput;

};

#endif
