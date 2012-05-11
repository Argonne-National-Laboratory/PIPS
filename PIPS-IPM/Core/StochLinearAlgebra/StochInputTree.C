#include "StochInputTree.h"
#include <malloc.h>
//***********************************************************
//************************** TREE ***************************
//***********************************************************
StochInputTree::StochInputTree() 
  : nodeInput(NULL) 
{ }

StochInputTree::StochInputTree(const StochInputNode& root)
{
  nodeInput = new StochInputNode(root);
}

StochInputTree::StochInputTree(StochInputNode* root)
{
  nodeInput = root;
}

StochInputTree::~StochInputTree() 
{
  if (nodeInput) delete nodeInput;
  for (size_t it=0; it<children.size(); it++) delete children[it];
}

void StochInputTree::AddChild(const StochInputNode  &node)
{
  children.push_back(new StochInputTree(node));
}

void StochInputTree::AddChild(StochInputTree *subTree)
{
  children.push_back(subTree);
}

//***********************************************************
//************************** NODE ***************************
//***********************************************************
StochInputTree::StochInputNode::
StochInputNode(void* user_data_, int id_, 
	       int n_, int my_, int mz_,
	       FMAT fQ_, FNNZ fnnzQ_, FVEC fc_, 
	       FMAT fA_, FNNZ fnnzA_, 
	       FMAT fB_, FNNZ fnnzB_, 
	       FVEC fb_, 
	       FMAT fC_, FNNZ fnnzC_,
	       FMAT fD_, FNNZ fnnzD_,
	       FVEC fclow_, FVEC ficlow_, 
	       FVEC fcupp_, FVEC ficupp_,
	       FVEC fxlow_, FVEC fixlow_, 
	       FVEC fxupp_, FVEC fixupp_,
	       bool deleteUserData_/*=false*/)
  : id(id_), n(n_), my(my_), mz(mz_), 
    nnzQ(-1), nnzA(-1), nnzB(-1), nnzC(-1), nnzD(-1),
    fnnzQ(fnnzQ_), fnnzA(fnnzA_), fnnzB(fnnzB_), fnnzC(fnnzC_), fnnzD(fnnzD_),
    fQ(fQ_), fA(fA_), fB(fB_), fC(fC_), fD(fD_),
    fc(fc_), fb(fb_), 
    fclow(fclow_), fcupp(fcupp_), ficlow(ficlow_), ficupp(ficupp_),
    fxlow(fxlow_), fxupp(fxupp_), fixlow(fixlow_), fixupp(fixupp_), 
    user_data(user_data_), 
    deleteUserData(deleteUserData_)
{ };

StochInputTree::StochInputNode::StochInputNode()
  : id(-1), n(-1), my(-1), mz(-1), 
    nnzQ(-1), nnzA(-1), nnzB(-1), nnzC(-1), nnzD(-1),
    user_data(NULL), 
    deleteUserData(false)
{ };


StochInputTree::StochInputNode::~StochInputNode() 
{ 
  if(deleteUserData) free(user_data);
};
