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

using namespace std;

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

	// full callback constructor without constraints
   StochInputNode(void* user_data, int id,
      FNNZ n, FNNZ my, FNNZ mz,
      FMAT fQ, FNNZ fnnzQ, FVEC fc,
      FMAT fA, FNNZ fnnzA, FMAT fB, FNNZ fnnzB,
      FVEC fb,
      FMAT fC, FNNZ fnnzC, FMAT fD, FNNZ fnnzD,
      FVEC fclow, FVEC ficlow, FVEC fcupp, FVEC ficupp,
      FVEC fxlow, FVEC fixlow, FVEC fxupp, FVEC fixupp,
      bool deleteUserData=false);

	// full callback constructor including linking constraints
	StochInputNode(void* user_data, int id,
	   FNNZ n, FNNZ my, FNNZ myl, FNNZ mz, FNNZ mzl,
	   FMAT fQ, FNNZ fnnzQ, FVEC fc,
	   FMAT fA, FNNZ fnnzA, FMAT fB, FNNZ fnnzB, FMAT fBl, FNNZ fnnzBl,
	   FVEC fb, FVEC fbl,
	   FMAT fC, FNNZ fnnzC, FMAT fD, FNNZ fnnzD, FMAT fDl, FNNZ fnnzDl,
	   FVEC fclow, FVEC ficlow, FVEC fcupp, FVEC ficupp,
	   FVEC fdllow, FVEC fidllow, FVEC fdlupp, FVEC fidlupp,
	   FVEC fxlow, FVEC fixlow, FVEC fxupp, FVEC fixupp,
	   bool deleteUserData=false);

  ~StochInputNode();

  protected:

    int id;
    int n, my, myl, mz, mzl;
    int nnzQ, nnzA, nnzB, nnzBl, nnzC, nnzD, nnzDl;

  protected:
    //callback functions

    //callback functions nCall, myCall, mzCall, mylCall, mzlCall can be NULL if data is provided through int n,my,mz,myl,mzl
    FNNZ nCall, myCall, mzCall, mylCall, mzlCall;
    FNNZ fnnzQ, fnnzA, fnnzB, fnnzBl, fnnzC, fnnzD, fnnzDl;
    FMAT fQ, fA, fB, fBl, fC, fD, fDl;

    FVEC fc, fb, fbl;
    FVEC fclow, fcupp, ficlow, ficupp;
    FVEC fdllow, fdlupp, fidllow, fidlupp;
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
  vector<StochInputTree*> children;
 protected:
  StochInputNode* nodeInput;

};

#endif
