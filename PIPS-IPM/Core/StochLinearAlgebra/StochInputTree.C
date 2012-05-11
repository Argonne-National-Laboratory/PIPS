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
  for (int it=0; it<children.size(); it++) delete children[it];
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
  : user_data(user_data_), 
    id(id_), n(n_), my(my_), mz(mz_), 
    fQ(fQ_), fnnzQ(fnnzQ_), 
    fA(fA_), fnnzA(fnnzA_), fB(fB_), fnnzB(fnnzB_), 
    fC(fC_),  fnnzC(fnnzC_), fD(fD_),  fnnzD(fnnzD_),
    fc(fc_), fb(fb_), 
    fxlow(fxlow_),  fixlow(fixlow_), fxupp(fxupp_),  fixupp(fixupp_), 
    fclow(fclow_),  ficlow(ficlow_), fcupp(fcupp_),  ficupp(ficupp_),
    nnzQ(-1), nnzA(-1), nnzB(-1), nnzC(-1), nnzD(-1),
    deleteUserData(deleteUserData_)
{ };

StochInputTree::StochInputNode::StochInputNode()
  : user_data(NULL), 
    id(-1), n(-1), my(-1), mz(-1), 
    nnzQ(-1), nnzA(-1), nnzB(-1), nnzC(-1), nnzD(-1),
    deleteUserData(false)
{ };


StochInputTree::StochInputNode::~StochInputNode() 
{ 
  if(deleteUserData) free(user_data);
};
