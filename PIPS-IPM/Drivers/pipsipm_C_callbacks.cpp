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


extern "C"
void PIPSSolve(MPI_Comm comm, void* user_data, int numScens,
		int nx0, int my0, int mz0, int nx, int my, int mz,
		FMAT fQ, FNNZ fnnzQ, FVEC fc,
		FMAT fA, FNNZ fnnzA, FMAT fB, FNNZ fnnzB,
		FVEC fb,
		FMAT fC, FNNZ fnnzC, FMAT fD, FNNZ fnnzD,
		FVEC fclow, FVEC ficlow, FVEC fcupp, FVEC ficupp,
		FVEC fxlow, FVEC fixlow, FVEC fxupp, FVEC fixupp,
        double* obj_val, double* first_primal, double* second_primal,
        double* first_dual, double* second_dual)
{
  //build the tree for specifying the problem
  int globalID=0;
  StochInputTree::StochInputNode data(user_data, globalID,
				      nx0,my0,mz0,
				      fQ, fnnzQ, fc,
				      fA, fnnzA,
				      fB,  fnnzB,
				      fb,
				      fC, fnnzC,
				      fD, fnnzD,
				      fclow, ficlow, fcupp, ficupp,
				      fxlow, fixlow, fxupp, fixupp );
  globalID++;
  StochInputTree* root = new StochInputTree(data);

  for(int i=0; i<numScens; i++) {
    StochInputTree::StochInputNode data(user_data, globalID,
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

  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  int nprocs; MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if(0==mype) cout << "Using a total of " << nprocs << " MPI processes." << endl;

  PIPSIpmInterface<sFactoryAug, MehrotraStochSolver> pipsIpm(root);
  if (mype == 0) cout << "PIPSIpmInterface created .." << endl;

  if (mype == 0) cout << "solving .." << endl;
  pipsIpm.go();
  if (mype == 0) cout << "solving ..solving done." << endl;

  obj_val[0] = pipsIpm.getObjective();
  for (int i=0;i<nx0;i++){ first_primal[i] = pipsIpm.getFirstStagePrimalColSolution()[i]; }
  for (int i=0;i<my+mz0;i++){ first_dual[i] = pipsIpm.getFirstStageDualRowSolution()[i]; }
  int cnt1 = 0;
  int cnt2 = 0;
  for (int i=0;i<numScens;i++){
    for (int j=0;j<nx;j++){
        second_primal[cnt1++] = pipsIpm.getSecondStagePrimalColSolution(i)[j];
    }
  //  for (int k=0;k<my+mz;k++){
  //      second_dual[cnt2++] = pipsIpm.getSecondStageDualRowSolution(i)[k];
  //  }
  }
  // second_primal = &(pipsIpm.getSecondStagePrimalColSolution()[0]);
  // first_dual    = &(pipsIpm.getFirstStageDualRowSolution()[0]);
  // second_dual   = &(pipsIpm.getSecondStageDualRowSolution()[0]);

  //uninitialization
  delete root;

}
