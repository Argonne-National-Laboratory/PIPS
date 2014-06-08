#include "StochInputTree.h"
#include "PIPSIpmInterface.h"
#include "sFactoryAug.h"
#include "MehrotraStochSolver.h"

#include "mpi.h"

extern "C" typedef int (*FNNZ)(void* user_data, int id, int* nnz);

/* Row-major format here */
extern "C" typedef int (*FMAT)(void* user_data, int id, 
			       int* krowM, int* jcolM, double* M);

extern "C" typedef int (*FVEC)(void* user_data, int id, double* vec, int len);

//extern "C" typedef int (*FLEN)(void* user_data, int id, int* len);


//extern "C" 

void PIPSSolve(MPI_Comm comm, int numScens, 
		int nx0, int my0, int mz0, int nx, int my, int mz,
		FMAT fQ0, FNNZ fnnzQ0, FVEC fc0,  
		FMAT fA0, FNNZ fnnzA0, FMAT fB0, FNNZ fnnzB0, 
		FVEC fb0, 
		FMAT fC0, FNNZ fnnzC0, FMAT fD0, FNNZ fnnzD0,
		FVEC fclow0, FVEC ficlow0, FVEC fcupp0, FVEC ficupp0,
		FVEC fxlow0, FVEC fixlow0, FVEC fxupp0, FVEC fixupp0,
		FMAT fQ, FNNZ fnnzQ, FVEC fc,  
		FMAT fA, FNNZ fnnzA, FMAT fB, FNNZ fnnzB, 
		FVEC fb, 
		FMAT fC, FNNZ fnnzC, FMAT fD, FNNZ fnnzD,
		FVEC fclow, FVEC ficlow, FVEC fcupp, FVEC ficupp,
		FVEC fxlow, FVEC fixlow, FVEC fxupp, FVEC fixupp)
{
  //build the tree for specifying the problem
  int globalID=0;
  StochInputTree::StochInputNode data(NULL, globalID, 
				      nx0,my0,mz0,
				      fQ0, fnnzQ0, fc0,
				      fA0, fnnzA0,
				      fB0,  fnnzB0,
				      fb0,
				      fC0, fnnzC0,
				      fD0, fnnzD0,
				      fclow0, ficlow0, fcupp0, ficupp0,
				      fxlow0, fixlow0, fxupp0, fixupp0 );
  globalID++;
  StochInputTree* root = new StochInputTree(data);
  
  for(int i=0; i<numScens; i++) {
    StochInputTree::StochInputNode data(NULL, globalID, 
					nx, my, mz, 
					fQ, fnnzQ, fc,
					fA, fnnzA,
					fB, fnnzB,
					fb,
					fC, fnnzC,
					fD, fnnzD,
					fclow, ficlow, fcupp, ficupp,
					fxlow, fixlow, fxupp, fixupp );
    globalID++;
    
    root->AddChild(new StochInputTree(data));
  }

  int mype; MPI_Comm_rank(comm,&mype);
  int nprocs; MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if(0==mype) cout << "Using a total of " << nprocs << " MPI processes." << endl;

  PIPSIpmInterface<sFactoryAug, MehrotraStochSolver> pipsIpm(root);
  if (mype == 0) cout << "PIPSIpmInterface created .." << endl;

  if (mype == 0) cout << "solving .." << endl;
  pipsIpm.go();
  if (mype == 0) cout << "solving ..solving done." << endl;

  //uninitialization
  delete root;

}
