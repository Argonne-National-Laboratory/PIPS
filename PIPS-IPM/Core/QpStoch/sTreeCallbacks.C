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
#include <cmath>
#include <algorithm>    // std::swap
#include "pipsport.h"

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
  : NNZA(0), NNZQ(0), NNZB(0), NNZBl(0), NNZC(0), NNZD(0), NNZDl(0),
    NNZA_INACTIVE(-1), NNZQ_INACTIVE(-1), NNZB_INACTIVE(-1), NNZBl_INACTIVE(-1), NNZC_INACTIVE(-1), NNZD_INACTIVE(-1), NNZDl_INACTIVE(-1),
    N_INACTIVE(-1), MY_INACTIVE(-1), MZ_INACTIVE(-1),
    nx_active(0),  my_active(0),  mz_active(0),  myl_active(0),  mzl_active(0),
    nx_inactive(-1),  my_inactive(-1),  mz_inactive(-1),  myl_inactive(-1),  mzl_inactive(-1),
    isDataPresolved(false), hasPresolvedData(false), data(nullptr), tree(nullptr), fakedata(nullptr)

{
  if(-1==rankMe) MPI_Comm_rank(MPI_COMM_WORLD, &rankMe);
  if(-1==numProcs) MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
}

sTreeCallbacks::sTreeCallbacks(StochInputTree* inputTree)
  : sTree(),
    NNZA(0), NNZQ(0), NNZB(0), NNZBl(0), NNZC(0), NNZD(0), NNZDl(0),
    NNZA_INACTIVE(-1), NNZQ_INACTIVE(-1), NNZB_INACTIVE(-1), NNZBl_INACTIVE(-1), NNZC_INACTIVE(-1), NNZD_INACTIVE(-1), NNZDl_INACTIVE(-1),
    N_INACTIVE(-1), MY_INACTIVE(-1), MZ_INACTIVE(-1),
    nx_active(0),  my_active(0),  mz_active(0),  myl_active(0),  mzl_active(0),
    nx_inactive(-1),  my_inactive(-1),  mz_inactive(-1),  myl_inactive(-1),  mzl_inactive(-1),
    isDataPresolved(false), hasPresolvedData(false), tree(nullptr), fakedata(nullptr)
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
    NNZA(0), NNZQ(0), NNZB(0), NNZBl(0), NNZC(0), NNZD(0), NNZDl(0),
    NNZA_INACTIVE(-1), NNZQ_INACTIVE(-1), NNZB_INACTIVE(-1), NNZBl_INACTIVE(-1), NNZC_INACTIVE(-1), NNZD_INACTIVE(-1), NNZDl_INACTIVE(-1),
    N_INACTIVE(-1), MY_INACTIVE(-1), MZ_INACTIVE(-1),
    nx_active(0),  my_active(0),  mz_active(0),  myl_active(0),  mzl_active(0),
    nx_inactive(-1),  my_inactive(-1),  mz_inactive(-1),  myl_inactive(-1),  mzl_inactive(-1),
    isDataPresolved(false), hasPresolvedData(false), data(nullptr), tree(nullptr), scens(localscens)
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
    NNZA(0), NNZQ(0), NNZB(0), NNZBl(0), NNZC(0), NNZD(0), NNZDl(0),
    NNZA_INACTIVE(-1), NNZQ_INACTIVE(-1), NNZB_INACTIVE(-1), NNZBl_INACTIVE(-1), NNZC_INACTIVE(-1), NNZD_INACTIVE(-1), NNZDl_INACTIVE(-1),
    N_INACTIVE(-1), MY_INACTIVE(-1), MZ_INACTIVE(-1),
    nx_active(data_->n),  my_active(data_->my),  mz_active(data_->mz),  myl_active(data_->myl),  mzl_active(data_->mzl),
    nx_inactive(-1),  my_inactive(-1),  mz_inactive(-1),  myl_inactive(-1),  mzl_inactive(-1),
    isDataPresolved(false), hasPresolvedData(false), data(data_), tree(nullptr), fakedata(nullptr)
{
  if(-1==rankMe) MPI_Comm_rank(MPI_COMM_WORLD, &rankMe);
  if(-1==numProcs) MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
}

void sTreeCallbacks::loadLocalSizes()
{
  //already in data structure
}

void sTreeCallbacks::switchToPresolvedData()
{
   assert(!isDataPresolved);
   assert(hasPresolvedData);

   std::swap(NNZA_INACTIVE, NNZA);
   std::swap(NNZQ_INACTIVE, NNZQ);
   std::swap(NNZB_INACTIVE, NNZB);
   std::swap(NNZBl_INACTIVE, NNZBl);
   std::swap(NNZC_INACTIVE, NNZC);
   std::swap(NNZD_INACTIVE, NNZD);
   std::swap(NNZDl_INACTIVE, NNZDl);
   std::swap(N_INACTIVE, N);
   std::swap(MY_INACTIVE, MY);
   std::swap(MZ_INACTIVE, MZ);

   std::swap(nx_active, nx_inactive);
   std::swap(my_active, my_inactive);
   std::swap(mz_active, mz_inactive);
   std::swap(myl_active, myl_inactive);
   std::swap(mzl_active, mzl_inactive);

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->np = this->nx_active;
      dynamic_cast<sTreeCallbacks*>(children[it])->switchToPresolvedData();
   }

   isDataPresolved = true;

}


void sTreeCallbacks::switchToOriginalData()
{
   assert(isDataPresolved);

   std::swap(NNZA_INACTIVE, NNZA);
   std::swap(NNZQ_INACTIVE, NNZQ);
   std::swap(NNZB_INACTIVE, NNZB);
   std::swap(NNZBl_INACTIVE, NNZBl);
   std::swap(NNZC_INACTIVE, NNZC);
   std::swap(NNZD_INACTIVE, NNZD);
   std::swap(NNZDl_INACTIVE, NNZDl);
   std::swap(N_INACTIVE, N);
   std::swap(MY_INACTIVE, MY);
   std::swap(MZ_INACTIVE, MZ);

   std::swap(nx_active, nx_inactive);
   std::swap(my_active, my_inactive);
   std::swap(mz_active, mz_inactive);
   std::swap(myl_active, myl_inactive);
   std::swap(mzl_active, mzl_inactive);

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->np = this->nx_active;
      dynamic_cast<sTreeCallbacks*>(children[it])->switchToOriginalData();
   }

   isDataPresolved = false;
}


bool sTreeCallbacks::isPresolved()
{
   return isDataPresolved;
}

bool sTreeCallbacks::hasPresolved()
{
   return hasPresolvedData;
}

void sTreeCallbacks::initPresolvedData(const StochSymMatrix& Q, const StochGenMatrix& A, const StochGenMatrix& C,
      const StochVector& nxVec, const StochVector& myVec, const StochVector& mzVec, int mylParent, int mzlParent)
{
   assert(!hasPresolvedData);

   assert(nxVec.children.size() == children.size());
   assert(myVec.children.size() == children.size());
   assert(mzVec.children.size() == children.size());
   assert(A.children.size() == children.size());
   assert(C.children.size() == children.size());

   assert(fakedata == nullptr);
   assert(tree == nullptr);

   const SimpleVector& nxVecSimple = dynamic_cast<const SimpleVector&>(*nxVec.vec);
   const SimpleVector& myVecSimple = dynamic_cast<const SimpleVector&>(*myVec.vec);
   const SimpleVector& mzVecSimple = dynamic_cast<const SimpleVector&>(*mzVec.vec);
   const SimpleVector* const myVecSimpleLink = dynamic_cast<const SimpleVector*>(myVec.vecl);
   const SimpleVector* const mzVecSimpleLink = dynamic_cast<const SimpleVector*>(mzVec.vecl);

   N_INACTIVE = nxVecSimple.n;
   MY_INACTIVE = myVecSimple.n;
   MZ_INACTIVE = mzVecSimple.n;

   nx_inactive = N_INACTIVE;
   my_inactive = MY_INACTIVE;
   mz_inactive = MZ_INACTIVE;

   if( myVecSimpleLink != nullptr )
   {
      assert(np == -1);
      myl_inactive = myVecSimpleLink->n;
   }
   else
   {
      myl_inactive = mylParent;
   }

   if( mzVecSimpleLink != nullptr )
   {
      assert(np == -1);
      mzl_inactive = mzVecSimpleLink->n;
   }
   else
   {
      mzl_inactive = mzlParent;
   }

   // empty child?
   if( N_INACTIVE == 0 && MY_INACTIVE == 0 && MZ_INACTIVE == 0 && np != -1 )
   {
      myl_inactive = 0;
      mzl_inactive = 0;
   }


   // todo
   NNZQ_INACTIVE = 0;

   // are we at the root?
   if( np == -1 )
   {
      assert(mylParent == -1);
      assert(mzlParent == -1);

      NNZA_INACTIVE = A.Bmat->numberOfNonZeros();
      NNZB_INACTIVE = A.Amat->numberOfNonZeros();
      NNZC_INACTIVE = C.Bmat->numberOfNonZeros();
      NNZD_INACTIVE = C.Amat->numberOfNonZeros();
   }
   else
   {
      NNZA_INACTIVE = A.Amat->numberOfNonZeros();
      NNZB_INACTIVE = A.Bmat->numberOfNonZeros();
      NNZC_INACTIVE = C.Amat->numberOfNonZeros();
      NNZD_INACTIVE = C.Bmat->numberOfNonZeros();
   }

   NNZBl_INACTIVE = A.Blmat->numberOfNonZeros();
   NNZDl_INACTIVE =  C.Blmat->numberOfNonZeros();

   for( size_t it = 0; it < children.size(); it++ )
   {
      assert(children[it]->np == this->nx_active);
      sTreeCallbacks* sTreeCallbacksChild = dynamic_cast<sTreeCallbacks*>(children[it]);

      sTreeCallbacksChild->initPresolvedData(*Q.children[it], *A.children[it], *C.children[it], *nxVec.children[it],
            *myVec.children[it], *mzVec.children[it], myl_inactive, mzl_inactive);

      N_INACTIVE += sTreeCallbacksChild->N_INACTIVE;
      MY_INACTIVE += sTreeCallbacksChild->MY_INACTIVE;
      MZ_INACTIVE += sTreeCallbacksChild->MZ_INACTIVE;

      //nnz stuff
      NNZQ_INACTIVE += sTreeCallbacksChild->NNZQ_INACTIVE;
      NNZA_INACTIVE += sTreeCallbacksChild->NNZA_INACTIVE;
      NNZB_INACTIVE += sTreeCallbacksChild->NNZB_INACTIVE;
      NNZBl_INACTIVE += sTreeCallbacksChild->NNZBl_INACTIVE;
      NNZC_INACTIVE += sTreeCallbacksChild->NNZC_INACTIVE;
      NNZD_INACTIVE += sTreeCallbacksChild->NNZD_INACTIVE;
      NNZDl_INACTIVE += sTreeCallbacksChild->NNZDl_INACTIVE;
   }

   hasPresolvedData = true;
}


void sTreeCallbacks::writeSizes(ostream& sout) const
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

#if 0
    sout << "NNZA : " << NNZA   << "\n";
    sout << "NNZQ : " << NNZQ   << "\n";
    sout << "NNZB : " << NNZB   << "\n";
    sout << "NNZBl : " <<  NNZBl   << "\n";
    sout << "NNZC : " << NNZC   << "\n";
    sout << "NNZD : " << NNZD   << "\n";
    sout << "NNZDl : " <<  NNZDl   << "\n";
#endif

    if( myRank == 0 )
    {
    sout << "N          : "  <<  N           << "\n";
    sout << "MY         : "  <<  MY          << "\n";
    sout << "MZ         : "  <<  MZ          << "\n";
    sout << "nx_active  : "  <<  nx_active    << "\n";
    sout << "my_active  : "  <<  my_active    << "\n";
    sout << "mz_active  : "  <<  mz_active    << "\n";
    sout << "myl_active : "  <<  myl_active   << "\n";
    sout << "mzl_active : "  <<  mzl_active   << "\n";
    }


   for( size_t it = 0; it < children.size(); it++ )
   {
      sout << "child " << it << ": \n\n";
      dynamic_cast<sTreeCallbacks*>(children[it])->writeSizes(sout);
   }
}


// this is usually called before assigning processes
void sTreeCallbacks::computeGlobalSizes()
{
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (data && sTree::isInVector(myrank, myProcs) ) {

    // callback used for sizes?
    if( data->nCall ) {
       assert(data->myCall);
       assert(data->mzCall);

       data->nCall(data->user_data, data->id, &data->n);
       data->myCall(data->user_data, data->id, &data->my);
       data->mzCall(data->user_data, data->id, &data->mz);

       if( data->mylCall )
          data->mylCall(data->user_data, data->id, &data->myl);
       else
          data->myl = -1;

       if( data->mzlCall )
          data->mzlCall(data->user_data, data->id, &data->mzl);
       else
          data->myl = -1;
    }

    N  = data->n;
    MY = data->my;
    MZ = data->mz;

    nx_active = data->n;
    my_active = data->my;
    mz_active = data->mz;
    myl_active = data->myl;
    mzl_active = data->mzl;

    NNZQ = data->nnzQ;
    NNZA = data->nnzA;
    NNZB = data->nnzB;
    NNZBl = data->nnzBl;
    NNZC = data->nnzC;
    NNZD = data->nnzD;
    NNZDl = data->nnzDl;
  } else {
    N = MY = MZ = NNZQ = NNZA = NNZB = NNZBl = NNZC = NNZD = NNZDl = 0;
  }
  if (tree && np == -1) {
    for(size_t it=0; it<tree->children.size();it++) {
      N += tree->children[it]->nodeInput->n;
      MY += tree->children[it]->nodeInput->my;
      MZ += tree->children[it]->nodeInput->mz;
      
      NNZQ += tree->children[it]->nodeInput->nnzQ;
      NNZA += tree->children[it]->nodeInput->nnzA;
      NNZB += tree->children[it]->nodeInput->nnzB;
      NNZBl += tree->children[it]->nodeInput->nnzBl;
      NNZC += tree->children[it]->nodeInput->nnzC;
      NNZD += tree->children[it]->nodeInput->nnzD;
      NNZDl += tree->children[it]->nodeInput->nnzDl;
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
      fakedata->nnzBl += scens[it]->nnzBl;
      fakedata->nnzC += scens[it]->nnzC;
      fakedata->nnzD += scens[it]->nnzD;
      fakedata->nnzDl += scens[it]->nnzDl;
      real_children[it]->np = np;
    }
    N += fakedata->n;
    MY += fakedata->my;
    MZ += fakedata->mz;
    NNZQ += fakedata->nnzQ;
    NNZA += fakedata->nnzA;
    NNZB += fakedata->nnzB;
    NNZBl += fakedata->nnzBl;
    NNZC += fakedata->nnzC;
    NNZD += fakedata->nnzD;
    NNZDl += fakedata->nnzDl;
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
    NNZBl += ((sTreeCallbacks*)children[it])->NNZBl;
    NNZC += ((sTreeCallbacks*)children[it])->NNZC;
    NNZD += ((sTreeCallbacks*)children[it])->NNZD;
    NNZDl += ((sTreeCallbacks*)children[it])->NNZDl;
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
    return nullptr;
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

  StochGenMatrix* A = nullptr;
  if (!fakedata) {

	 if (data->fnnzBl && data->nnzBl < 0)
       data->fnnzBl(data->user_data, data->id, &data->nnzBl);

    if (data->nnzA<0)
      data->fnnzA(data->user_data, data->id, &data->nnzA);

	// are we at the root?
    if (np==-1) {

      data->nnzB=0;

      // are there linking constraints?
      if (data->fnnzBl)
      {
    	// populate B with A's data B_0 is the A_0 from the theoretical form; also fill Bl
    	// (i.e. the first block of linking constraints)
        A = new StochGenMatrix(data->id,
             MY, N,
             data->my, np, data->nnzB,
             data->my, data->n,  data->nnzA,
             data->myl, data->n,  data->nnzBl,
             commWrkrs);
      }
      else
      {
    	// populate B with A's data B_0 is the A_0 from the theoretical form
        A = new StochGenMatrix(data->id,
             MY, N,
             data->my, np, data->nnzB,
             data->my, data->n,  data->nnzA,
             commWrkrs);
      }

      //populate submatrix B
      data->fA(data->user_data, data->id, A->Bmat->krowM(), A->Bmat->jcolM(), A->Bmat->M());


      printf("root  -- my=%d  myl=%d nx=%d   1st stg nx=%d nnzA=%d nnzB=%d, nnzBl=%d\n",
        data->my, data->myl, data->n,  np, data->nnzA, data->nnzB, data->nnzBl);
    } else {

      if (data->nnzB<0)
        data->fnnzB(data->user_data, data->id, &data->nnzB);

      // are there linking constraints?
      if (data->fnnzBl)
      {
        A = new StochGenMatrix(data->id,
             MY, N,
             data->my, np, data->nnzA,
             data->my, data->n,  data->nnzB,
	  	       data->myl, data->n,  data->nnzBl,
             commWrkrs);
      }
      else
      {
        A = new StochGenMatrix(data->id,
             MY, N,
             data->my, np, data->nnzA,
             data->my, data->n,  data->nnzB,
             commWrkrs);
      }
      //populate the submatrices A, B
      data->fA(data->user_data, data->id, A->Amat->krowM(), A->Amat->jcolM(), A->Amat->M());
      data->fB(data->user_data, data->id, A->Bmat->krowM(), A->Bmat->jcolM(), A->Bmat->M());

      printf("  -- my=%d  myl=%d nx=%d   1st stg nx=%d nnzA=%d nnzB=%d, nnzBl=%d\n",
	     data->my, data->myl, data->n,  np, data->nnzA, data->nnzB, data->nnzBl);
    }

    // populate Bl if existent
    if (data->fBl)
      data->fBl(data->user_data, data->id, A->Blmat->krowM(), A->Blmat->jcolM(), A->Blmat->M());

    for(size_t it=0; it<children.size(); it++) {
      StochGenMatrix* child = children[it]->createA();
      A->AddChild(child);
    }
  } else {
    assert(false);
    return nullptr;

  }
  return A;
}

StochGenMatrix* sTreeCallbacks::createC() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochGenDummyMatrix(id());

  StochGenMatrix* C = nullptr;
  if (!fakedata) {

    if (data->fnnzDl && data->nnzDl < 0)
       data->fnnzDl(data->user_data, data->id, &data->nnzDl);

    if(data->nnzC<0)
      data->fnnzC(data->user_data, data->id, &data->nnzC);

    // at the root?
    if (np==-1) {

      data->nnzD=0;
      
      // are there linking constraints?
      if (data->fnnzDl)
      {
      // populate D with C's data D_0 is the C_0 from the theoretical form; also fill Dl
      // (i.e. the first block of linking constraints)
        C = new StochGenMatrix(data->id,
             MZ, N,
             data->mz, np, data->nnzD,
             data->mz, data->n, data->nnzC,
             data->mzl, data->n, data->nnzDl,
             commWrkrs);
      }
      else
      {
        //populate D with C's data
        //D_0 is the C_0 from the theoretical form
        C = new StochGenMatrix(data->id,
             MZ, N,
             data->mz, np, data->nnzD,
             data->mz, data->n,  data->nnzC,
             commWrkrs);
      }

      data->fC(data->user_data, data->id, C->Bmat->krowM(), C->Bmat->jcolM(), C->Bmat->M());

      printf("root  -- mz=%d  mzl=%d nx=%d 1st stg nx=%d nnzD=%d\n",
	     data->mz, data->mzl, data->n,  np, data->nnzC);
     
    } else {

      if(data->nnzD<0)
        data->fnnzD(data->user_data, data->id, &data->nnzD);

      // are there linking constraints?
      if (data->fnnzDl)
      {
        C = new StochGenMatrix(data->id,
             MZ, N,
             data->mz, np, data->nnzC,
             data->mz, data->n, data->nnzD,
             data->mzl, data->n, data->nnzDl,
             commWrkrs);
      }
      else
      {
        C = new StochGenMatrix(data->id,
             MZ, N,
             data->mz, np, data->nnzC,
             data->mz, data->n, data->nnzD,
             commWrkrs);
      }

      //populate the submatrices C and D
      data->fC(data->user_data, data->id, C->Amat->krowM(), C->Amat->jcolM(), C->Amat->M());
      data->fD(data->user_data, data->id, C->Bmat->krowM(), C->Bmat->jcolM(), C->Bmat->M());

      printf("  -- mz=%d mzl=%d nx=%d  1st stg nx=%d nnzC=%d nnzD=%d, nnzDl=%d\n",
	     data->mz, data->mzl, data->n,  np, data->nnzC, data->nnzD, data->nnzDl);
    }
      
    // populate Dl if existent
    if (data->fDl)
      data->fDl(data->user_data, data->id, C->Blmat->krowM(), C->Blmat->jcolM(), C->Blmat->M());

    for(size_t it=0; it<children.size(); it++) {
      StochGenMatrix* child = children[it]->createC();
      C->AddChild(child);
    }
  } else {
    assert(false);
    return nullptr;
  }
  return C;
}

int sTreeCallbacks::nx() const {
   if( data )
      return nx_active;
   else
      return fakedata->n;
}

int sTreeCallbacks::my() const {
   if( data )
      return my_active;
   else
      return fakedata->my;
}

int sTreeCallbacks::myl() const {
   if( data )
      return myl_active;
   else
      return fakedata->myl;
}

int sTreeCallbacks::mz() const {
   if( data )
      return mz_active;
   else
      return fakedata->mz;
}

int sTreeCallbacks::mzl() const {
   if( data )
      return mzl_active;
   else
      return fakedata->mzl;
}

// not sure what this is used for
int sTreeCallbacks::id() const {
  if (data) return data->id;
  else return 0;
}

static double RESCALE=1.0;
StochVector* sTreeCallbacks::createc() const
{

  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();  

  StochVector* c = new StochVector(nx(), commWrkrs);
  double* vData = ((SimpleVector*)c->vec)->elements();  
  if (!fakedata) {
    // populate the node's data with data from user.
    //if(children.size()>0) RESCALE=0.001/1000;//children.size();
#if 0
    if(0==rankMe)
      cout << "RESCALE set to " << RESCALE << endl;
#endif

    data->fc(data->user_data, data->id, 
       vData, data->n);
    
    for(int i=0; i<data->n; i++) vData[i]*=RESCALE;

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

  int yl = (np == -1) ? myl() : -1;

  StochVector* b = new StochVector(my(), yl, commWrkrs, -1);

  double* vData = ((SimpleVector*)b->vec)->elements();
  double* vDataLinkCons = nullptr;

  if (np == -1 && b->vecl )
     vDataLinkCons = ((SimpleVector*)b->vecl)->elements();

  if (!fakedata) {
    data->fb(data->user_data, data->id, 
       vData, data->my);

    // at root and with linking constraints?
    if (np == -1 && data->fbl)
    {
      assert(vDataLinkCons);
      data->fbl(data->user_data, data->id, vDataLinkCons, data->myl);
    }

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

    pos = 0;
    // at root and with linking constraints? todo don't really know whether this is correct
    if (np == -1 && data->fbl) {
        assert(vDataLinkCons);
        scens[0]->fbl(scens[0]->user_data, scens[0]->id, vDataLinkCons, scens[0]->myl);
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

  int zl = (np == -1) ? mzl() : -1;

  StochVector* clow = new StochVector(mz(), zl, commWrkrs, -1);
  double* vData = ((SimpleVector*)clow->vec)->elements();
  double* vDataLinkCons = nullptr;

  if (np == -1 && clow->vecl )
     vDataLinkCons = ((SimpleVector*)clow->vecl)->elements();

  if (!fakedata) {  
    data->fclow(data->user_data, data->id, vData, data->mz);

    // at root and with linking constraints?
    if (np == -1 && data->fdllow)
    {
      assert(vDataLinkCons);
      data->fdllow(data->user_data, data->id, vDataLinkCons, data->mzl);
    }

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

    // at root and with linking constraints? todo don't really know whether this is correct
    if (np == -1 && data->fdllow)
    {
        assert(vDataLinkCons);
        scens[0]->fdllow(scens[0]->user_data, scens[0]->id, vDataLinkCons, scens[0]->mzl);
    }

  }
  return clow;
}

StochVector* sTreeCallbacks::createiclow() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  int zl = (np == -1) ? mzl() : -1;

  StochVector* iclow = new StochVector(mz(), zl, commWrkrs, -1);
  double* vData = ((SimpleVector*)iclow->vec)->elements();  
  double* vDataLinkCons = nullptr;

  if (np == -1 && iclow->vecl )
    vDataLinkCons = ((SimpleVector*)iclow->vecl)->elements();

  if (!fakedata) {
    data->ficlow(data->user_data, data->id,  vData, data->mz);

    // at root and with linking constraints?
    if (np == -1 && data->fidllow)
    {
      assert(vDataLinkCons);
      data->fidllow(data->user_data, data->id, vDataLinkCons, data->mzl);
    }

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

    // at root and with linking constraints? todo don't really know whether this is correct
    if (np == -1 && data->fidllow)
    {
        assert(vDataLinkCons);
        scens[0]->fidllow(scens[0]->user_data, scens[0]->id, vDataLinkCons, scens[0]->mzl);
    }
  }
  return iclow;
}

StochVector* sTreeCallbacks::createcupp() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  int zl = (np == -1) ? mzl() : -1;

  StochVector* cupp = new StochVector(mz(), zl, commWrkrs, -1);
  double* vData = ((SimpleVector*)cupp->vec)->elements();  
  double* vDataLinkCons = nullptr;

  if (np == -1 && cupp->vecl)
    vDataLinkCons = ((SimpleVector*)cupp->vecl)->elements();

  if (!fakedata) {
    data->fcupp(data->user_data, data->id, vData, data->mz);

    // at root and with linking constraints?
    if (np == -1 && data->fdlupp)
    {
      assert(vDataLinkCons);
      data->fdlupp(data->user_data, data->id, vDataLinkCons, data->mzl);
    }

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

    // at root and with linking constraints? todo don't really know whether this is correct
    if (np == -1 && data->fdlupp)
    {
      assert(vDataLinkCons);
      scens[0]->fdlupp(scens[0]->user_data, scens[0]->id, vDataLinkCons, scens[0]->mzl);
    }

  }
  return cupp;
}

StochVector* sTreeCallbacks::createicupp() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  int zl = (np == -1) ? mzl() : -1;

  StochVector* icupp = new StochVector(mz(), zl, commWrkrs, -1);
  double* vData = ((SimpleVector*)icupp->vec)->elements();
  double* vDataLinkCons = nullptr;

  if (np == -1 && icupp->vecl)
     vDataLinkCons = ((SimpleVector*)icupp->vecl)->elements();

  if (!fakedata) {
    data->ficupp(data->user_data, data->id, 
       vData, data->mz);

    // at root and with linking constraints?
    if (np == -1 && data->fidlupp)
    {
      assert(vDataLinkCons);
      data->fidlupp(data->user_data, data->id, vDataLinkCons, data->mzl);
    }

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

    // at root and with linking constraints? todo don't really know whether this is correct
    if (np == -1 && data->fidlupp)
    {
      assert(vDataLinkCons);
      scens[0]->fidlupp(scens[0]->user_data, scens[0]->id, vDataLinkCons, scens[0]->mzl);
    }
  }
  return icupp;
}

