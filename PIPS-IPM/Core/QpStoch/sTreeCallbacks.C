/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#include "sTreeCallbacks.h"
#include "sData.h"

#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "SimpleVector.h"
#include "DoubleMatrixTypes.h"
#include "StochResourcePlanner.h"
#include <cmath>

using namespace std;

#ifndef UCTRANS // see note in smlParDriver.C
#define UCTRANS
#endif

sTreeCallbacks::~sTreeCallbacks()
{
  for(size_t it=0; it<children.size(); it++)
    if (fakedata) delete fakedata;
  for(size_t i = 0; i < real_children.size(); i++)
    delete real_children[i];
}

sTreeCallbacks::sTreeCallbacks() 
  : NNZA(0), NNZQ(0), NNZB(0), NNZC(0), NNZD(0),
    data(NULL), tree(NULL), fakedata(NULL)
{
  if(-1==rankMe) MPI_Comm_rank(MPI_COMM_WORLD, &rankMe);
  if(-1==numProcs) MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
}

sTreeCallbacks::sTreeCallbacks(StochInputTree* inputTree)
  : sTree(),
    NNZA(0), NNZQ(0), NNZB(0), NNZC(0), NNZD(0),
    tree(NULL), fakedata(NULL)
{
  if(-1==rankMe) MPI_Comm_rank(MPI_COMM_WORLD, &rankMe);
  if(-1==numProcs) MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  data = inputTree->nodeInput;
#ifndef POOLSCEN
  for(size_t it=0; it<inputTree->children.size(); it++)
    children.push_back(new sTreeCallbacks(inputTree->children[it]));
#else
	tree = inputTree;
#endif
}

// np==-1 is used to indicate the root node. these can't be root nodes
sTreeCallbacks::sTreeCallbacks(const vector<StochInputTree::StochInputNode*> &localscens)
  : sTree(), 
    NNZA(0), NNZQ(0), NNZB(0), NNZC(0), NNZD(0),
    data(NULL), tree(NULL), scens(localscens)
{
  if(-1==rankMe) MPI_Comm_rank(MPI_COMM_WORLD, &rankMe);
  if(-1==numProcs) MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	fakedata = new StochInputTree::StochInputNode();
	real_children.reserve(scens.size());
	for(size_t i = 0; i < scens.size(); i++) {
		real_children.push_back(new sTreeCallbacks(scens[i]));

	}
}
sTreeCallbacks::sTreeCallbacks(StochInputTree::StochInputNode* data_)
  : sTree(), 
    NNZA(0), NNZQ(0), NNZB(0), NNZC(0), NNZD(0),
    data(data_), tree(NULL), fakedata(NULL)
{
  if(-1==rankMe) MPI_Comm_rank(MPI_COMM_WORLD, &rankMe);
  if(-1==numProcs) MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
}

// this is usually called before assigning processes
void sTreeCallbacks::computeGlobalSizes()
{
  if (data) {
    N  = data->n;
    MY = data->my;
    MZ = data->mz;
    
    NNZQ = data->nnzQ;
    NNZA = data->nnzA;
    NNZB = data->nnzB;
    NNZC = data->nnzC;
    NNZD = data->nnzD;
  } else {
    N = MY = MZ = NNZQ = NNZA = NNZB = NNZC = NNZD = 0;
  }
  if (tree && np == -1) {
    for(size_t it=0; it<tree->children.size();it++) {
      N += tree->children[it]->nodeInput->n;
      MY += tree->children[it]->nodeInput->my;
      MZ += tree->children[it]->nodeInput->mz;
      
      NNZQ += tree->children[it]->nodeInput->nnzQ;
      NNZA += tree->children[it]->nodeInput->nnzA;
      NNZB += tree->children[it]->nodeInput->nnzB;
      NNZC += tree->children[it]->nodeInput->nnzC;
      NNZD += tree->children[it]->nodeInput->nnzD;
    }
  } else if (fakedata) {
    fakedata->n = fakedata->my = fakedata->mz = fakedata->nnzQ = fakedata->nnzA = fakedata->nnzB = fakedata->nnzC = fakedata->nnzD = 0;	
    for(size_t it=0; it<scens.size();it++) {
      fakedata->n += scens[it]->n;
      fakedata->my += scens[it]->my;
      fakedata->mz += scens[it]->mz;
      fakedata->nnzQ += scens[it]->nnzQ;
      fakedata->nnzA += scens[it]->nnzA;
      fakedata->nnzB += scens[it]->nnzB;
      fakedata->nnzC += scens[it]->nnzC;
      fakedata->nnzD += scens[it]->nnzD;
      real_children[it]->np = np;
    }
    N += fakedata->n;
    MY += fakedata->my;
    MZ += fakedata->mz;
    NNZQ += fakedata->nnzQ;
    NNZA += fakedata->nnzA;
    NNZB += fakedata->nnzB;
    NNZC += fakedata->nnzC;
    NNZD += fakedata->nnzD;
  }
  for(size_t it=0; it<children.size(); it++) {
    children[it]->np = this->data->n;
    children[it]->computeGlobalSizes();
    N  += children[it]->N;
    MY += children[it]->MY;
    MZ += children[it]->MZ;
    
    //nnz stuff
    NNZQ += ((sTreeCallbacks*)children[it])->NNZQ;
    NNZA += ((sTreeCallbacks*)children[it])->NNZA;
    NNZB += ((sTreeCallbacks*)children[it])->NNZB;
    NNZC += ((sTreeCallbacks*)children[it])->NNZC;
    NNZD += ((sTreeCallbacks*)children[it])->NNZD;
  }
}

StochSymMatrix* sTreeCallbacks::createQ() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochSymDummyMatrix(id());
  
  if (!fakedata) {
    if(data->nnzQ<0)
      data->fnnzQ(data->user_data, data->id, &data->nnzQ);
    
    StochSymMatrix* Q = 
      new StochSymMatrix(data->id, 
			 N, 
			 data->n, 
			 data->nnzQ,
			 commWrkrs);
    
    data->fQ(data->user_data, data->id, 
	     Q->diag->krowM(), 
	     Q->diag->jcolM(),
	     Q->diag->M());  
    
    for(size_t it=0; it<children.size(); it++) {
      StochSymMatrix* child = children[it]->createQ();
      Q->AddChild(child);
    }
    return Q;
  } else {
    assert(false);
    return NULL;
    /*
    assert(real_children.size() > 0);
    vector<StochSymMatrix*> v(real_children.size());
    for(size_t i = 0; i<real_children.size(); i++) {
      v[i] = real_children[i]->createQ();
    }
    StochSymMatrix *out = new StochSymMatrix(v);
    for(size_t i = 0; i<real_children.size(); i++) delete v[i];
    return out;
    */
  }
}




StochGenMatrix* sTreeCallbacks::createA() const
{

  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL) {
    return new StochGenDummyMatrix(id());
  }

  StochGenMatrix* A = NULL;
  if (!fakedata) {
    if (np==-1) {

      //data->fnnzA(data->user_data, data->id, &nnzA);
      if (data->nnzA<0)
        data->fnnzA(data->user_data, data->id, &data->nnzA);
      data->nnzB=0;

      //this is the root; populate B with A's data
      //B_0 is the A_0 from the theoretical form
      A = new StochGenMatrix(data->id, 
           N, MY, 
           data->my, np, data->nnzB,
           data->my, data->n,  data->nnzA,
           commWrkrs);
      //populate the submatrices A and B
      //data->fA(data->user_data, data->id, A->Amat->krowM(), A->Amat->jcolM(), A->Amat->M());
      data->fA(data->user_data, data->id, A->Bmat->krowM(), A->Bmat->jcolM(), A->Bmat->M());
    } else {

      if (data->nnzA<0)
        data->fnnzA(data->user_data, data->id, &data->nnzA);
      if (data->nnzB<0)
        data->fnnzB(data->user_data, data->id, &data->nnzB);

      A = new StochGenMatrix(data->id, 
           N, MY, 
           data->my, np, data->nnzA, 
           data->my, data->n,  data->nnzB,
           commWrkrs);
      //populate the submatrices A and B
      data->fA(data->user_data, data->id, A->Amat->krowM(), A->Amat->jcolM(), A->Amat->M());
      data->fB(data->user_data, data->id, A->Bmat->krowM(), A->Bmat->jcolM(), A->Bmat->M());

      printf("  -- my=%d nx=%d   1st stg nx=%d nnzA=%d nnzB=%d\n", 
	     data->my, data->n,  np, data->nnzA, data->nnzB);
    }

    for(size_t it=0; it<children.size(); it++) {
      StochGenMatrix* child = children[it]->createA();
      A->AddChild(child);
    }
  } else {
    assert(false);
    return NULL;

//     vector<StochGenMatrix*> v(real_children.size());
// #ifdef UCTRANS
//     v[0] = real_children[0]->createA();
//     for(size_t i = 1; i<real_children.size();i++) v[i] = v[0];
// #else
//     for(size_t i = 0; i<real_children.size(); i++) {
//       v[i] = real_children[i]->createA();
//     }
// #endif
//     A = new StochGenMatrix(v);
// #ifdef UCTRANS
//     delete v[0];
// #else
//     for(size_t i = 0; i<real_children.size(); i++) delete v[i];
// #endif

  }
  return A;
}

StochGenMatrix* sTreeCallbacks::createC() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochGenDummyMatrix(id());

  printf("Create C\n");
  StochGenMatrix* C = NULL;
  if (!fakedata) {
    if (np==-1) {

      if(data->nnzC<0)
        data->fnnzC(data->user_data, data->id, &data->nnzC);
      data->nnzD=0;
      
      //this is the root; populate D with C's data
      //D_0 is the C_0 from the theoretical form
      C = new StochGenMatrix(data->id, 
           N, MZ, 
           data->mz, np, data->nnzD, 
           data->mz, data->n,  data->nnzC, commWrkrs);

      data->fC(data->user_data, data->id, C->Bmat->krowM(), C->Bmat->jcolM(), C->Bmat->M());

      printf("  -- mz=%d nx=%d   1st stg nx=%d nnzD=%d\n", 
	     data->mz, data->n,  np, data->nnzC);
     
    } else {
      if(data->nnzC<0)
        data->fnnzC(data->user_data, data->id, &data->nnzC);
      if(data->nnzD<0)
        data->fnnzD(data->user_data, data->id, &data->nnzD);

      C = new StochGenMatrix(data->id, 
           N, MZ, 
           data->mz, np, data->nnzC, 
           data->mz, data->n,  data->nnzD,
           commWrkrs);

      //populate the submatrices C and D
      data->fC(data->user_data, data->id, C->Amat->krowM(), C->Amat->jcolM(), C->Amat->M());
      data->fD(data->user_data, data->id, C->Bmat->krowM(), C->Bmat->jcolM(), C->Bmat->M());

      printf("  -- mz=%d nx=%d   1st stg nx=%d nnzC=%d nnzD=%d\n", 
	     data->mz, data->n,  np, data->nnzC, data->nnzD);
    }
      
    for(size_t it=0; it<children.size(); it++) {
      StochGenMatrix* child = children[it]->createC();
      C->AddChild(child);
    }
  } else {
    assert(false);
    return NULL;
//     vector<StochGenMatrix*> v(real_children.size());
// #ifdef UCTRANS
//     v[0] = real_children[0]->createC();
//     for(size_t i = 1; i<real_children.size();i++) v[i] = v[0];
// #else
//     for(size_t i = 0; i<real_children.size(); i++) {
//       v[i] = real_children[i]->createC();
//     }
// #endif
//     C = new StochGenMatrix(v);
// #ifdef UCTRANS
//     delete v[0];
// #else
//     for(size_t i = 0; i<real_children.size(); i++) delete v[i];
// #endif
  }
  return C;
}

int sTreeCallbacks::nx() const {
  if (data) return data->n;
  else return fakedata->n;
}

int sTreeCallbacks::my() const {
  if (data) return data->my;
  else return fakedata->my;
}

int sTreeCallbacks::mz() const {
  if (data) return data->mz;
  else return fakedata->mz;
}
// not sure what this is used for
int sTreeCallbacks::id() const {
  if (data) return data->id;
  else return 0;
}

StochVector* sTreeCallbacks::createc() const
{

  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();  

  StochVector* c = new StochVector(nx(), commWrkrs);
  double* vData = ((SimpleVector*)c->vec)->elements();  
  if (!fakedata) {
    // populate the node's data with data from user.

    data->fc(data->user_data, data->id, 
       vData, data->n);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createc();
      c->AddChild(child);
    }
  } else {
#ifdef UCTRANS
    int n = scens[0]->n;
    scens[0]->fc(scens[0]->user_data,scens[0]->id,
        vData, scens[0]->n);
    for(size_t i = 1; i < scens.size(); i++) {
      memcpy(vData+i*n,vData,n*sizeof(double));
    }
#else
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fc(scens[i]->user_data,scens[i]->id,
        vData+pos, scens[i]->n);
      pos += scens[i]->n;
    }
#endif
  }

  return c;
}

StochVector* sTreeCallbacks::createb() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* b = new StochVector(my(), commWrkrs);
  double* vData = ((SimpleVector*)b->vec)->elements();  
  if (!fakedata) {
    data->fb(data->user_data, data->id, 
       vData, data->my);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createb();
      b->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fb(scens[i]->user_data,scens[i]->id,
        vData+pos, scens[i]->my);
      pos += scens[i]->my;
    }
  }

  return b;
}

StochVector* sTreeCallbacks::createxlow() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* xlow = new StochVector(nx(), commWrkrs);
  double* vData = ((SimpleVector*)xlow->vec)->elements();  
  if (!fakedata) {
    data->fxlow(data->user_data, data->id, 
       vData, data->n);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createxlow();
      xlow->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fxlow(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->n);
      pos += scens[i]->n;
    }
  }
  return xlow;
}

StochVector* sTreeCallbacks::createixlow() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* ixlow = new StochVector(nx(), commWrkrs);
  double* vData = ((SimpleVector*)ixlow->vec)->elements();  
  if (!fakedata) {
    data->fixlow(data->user_data, data->id, 
       vData, data->n);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createixlow();
      ixlow->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fixlow(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->n);
      pos += scens[i]->n;
    }
  }
  return ixlow;
}

StochVector* sTreeCallbacks::createxupp() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* xupp = new StochVector(nx(), commWrkrs);
  double* vData = ((SimpleVector*)xupp->vec)->elements();  
  if (!fakedata) {
    data->fxupp(data->user_data, data->id, 
       vData, data->n);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createxupp();
      xupp->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fxupp(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->n);
      pos += scens[i]->n;
    }
  }
  return xupp;
}

StochVector* sTreeCallbacks::createixupp() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* ixupp = new StochVector(nx(), commWrkrs);
  double* vData = ((SimpleVector*)ixupp->vec)->elements();  
  if (!fakedata) {
    data->fixupp(data->user_data, data->id, 
       vData, data->n);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createixupp();
      ixupp->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fixupp(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->n);
      pos += scens[i]->n;
    }

  }

  return ixupp;
}

StochVector* sTreeCallbacks::createclow() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* clow = new StochVector(mz(), commWrkrs);
  double* vData = ((SimpleVector*)clow->vec)->elements();  
  if (!fakedata) {  
    data->fclow(data->user_data, data->id, 
         vData, data->mz);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createclow();
      clow->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fclow(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->mz);
      pos += scens[i]->mz;
    }
  }
  return clow;
}

StochVector* sTreeCallbacks::createiclow() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* iclow = new StochVector(mz(), commWrkrs);
  double* vData = ((SimpleVector*)iclow->vec)->elements();  
  if (!fakedata) {
    data->ficlow(data->user_data, data->id, 
       vData, data->mz);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createiclow();
      iclow->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->ficlow(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->mz);
      pos += scens[i]->mz;
    }
  }
  return iclow;
}

StochVector* sTreeCallbacks::createcupp() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* cupp = new StochVector(mz(), commWrkrs);
  double* vData = ((SimpleVector*)cupp->vec)->elements();  
    if (!fakedata) {
    data->fcupp(data->user_data, data->id, 
       vData, data->mz);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createcupp();
      cupp->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fcupp(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->mz);
      pos += scens[i]->mz;
    }

  }
  return cupp;
}

StochVector* sTreeCallbacks::createicupp() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* icupp = new StochVector(mz(), commWrkrs);
  double* vData = ((SimpleVector*)icupp->vec)->elements();  
  if (!fakedata) {
    data->ficupp(data->user_data, data->id, 
       vData, data->mz);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createicupp();
      icupp->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->ficupp(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->mz);
      pos += scens[i]->mz;
    }
  }
  return icupp;
}

