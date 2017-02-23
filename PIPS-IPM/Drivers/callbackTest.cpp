#include "StochInputTree.h"
#include "PIPSIpmInterface.h"
#include "sFactoryAug.h"
#include "MehrotraStochSolver.h"

#include "mpi.h"

#define LINKING_CONS 1

extern "C" typedef int (*FNNZ)(void* user_data, int id, int* nnz);

/* Row-major format */
extern "C" typedef int (*FMAT)(void* user_data, int id, int* krowM, int* jcolM, double* M);

extern "C" typedef int (*FVEC)(void* user_data, int id, double* vec, int len);


/** Problem parameters and data */
struct ProbData
{
public: //data
  int nScenarios;

public: //methods
  ProbData(int nScenarios)
  {
     this->nScenarios = nScenarios;
  };

//  ~ProbData();
};

extern "C" {

int nnzMatEqStage1(void* user_data, int id, int* nnz)
{
	if( id == 0 )
		*nnz = 2;
	else
		*nnz = 2;

	return 0;
}

int nnzMatEqStage2(void* user_data, int id, int* nnz)
{
	if( id == 0 )
		*nnz = 0;
	else
		*nnz = 2;
	return 0;
}

int nnzMatEqLink(void* user_data, int id, int* nnz)
{
	*nnz = 1;

	return 0;
}

int nnzAllZero(void* user_data, int id, int* nnz)
{
	*nnz = 0;
	return 0;
}

int vecAllZero(void* user_data, int id, double* vec, int len)
{
	int i;
	for( i = 0; i < len; i++ )
		vec[i] = 0.0;

	return 0;
}

int vecEqRhs(void* user_data, int id, double* vec, int len)
{
	int i;
	if( id == 0 )
	{
		for( i = 0; i < len; i++ )
			vec[i] = 1.0;
	}
	else
	{
		for( i = 0; i < len; i++ )
			vec[i] =  2.0 * id + (double) i;
	}

	return 0;
}

int vecObj(void* user_data, int id, double* vec, int len)
{
	int i;
	for( i = 0; i < len; i++ )
		vec[i] = 1.0;

	return 0;
}

int vecXlb(void* user_data, int id, double* vec, int len)
{
	int i;
	for( i = 0; i < len; i++ )
		vec[i] = 0.0;

	return 0;
}

int vecXlbActive(void* user_data, int id, double* vec, int len)
{
	int i;
	for( i = 0; i < len; i++ )
		vec[i] = 1.0;

	return 0;
}


int vecLinkRhs(void* user_data, int id, double* vec, int len)
{
	int i;

	for( i = 0; i < len; i++ )
		vec[i] = 6.0;

	return 0;
}

int matAllZero(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
	int i;
    int n = 2;

    for( i = 0; i <= n; i++ )
        krowM[i] = 0;

    return 0;
}

int matEqStage1(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
	int i;
    int n = 2;

    if( id == 0 )
    {

  	   for( i = 0; i < n; i++ )
  		   M[i] = 1.0;

  	   krowM[0] = 0;
  	   krowM[1] = 1;
  	   krowM[2] = 2;

  	   jcolM[0] = 0;
  	   jcolM[1] = 1;
#if 0
  	   for( i = 0; i <= n; i++ )
  		   krowM[i] = 0;
#endif
    }
    else
    {
 	   for( i = 0; i < n; i++ )
 		   M[i] = 1.0;

 	   krowM[0] = 0;
 	   krowM[1] = 1;
 	   krowM[2] = 2;

 	   jcolM[0] = 0;
 	   jcolM[1] = 1;
    }

    return 0;
}

int matEqStage2(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
	int i;
    int n = 2;

    if( id == 0 )
    {
   	   for( i = 0; i <= n; i++ )
   		   krowM[i] = 0;
    }
    else if( id == 1 )
    {
 	   for( i = 0; i < 2; i++ )
 		   M[i] = (double) i + 1.0;

 	   krowM[0] = 0;
 	   krowM[1] = 1;
 	   krowM[2] = 2;

 	   jcolM[0] = 0;
 	   jcolM[1] = 1;
    }
    else
    {
  	   for( i = 0; i < 2; i++ )
  		   M[i] = (double) i + 3.0;

  	   krowM[0] = 0;
  	   krowM[1] = 1;
  	   krowM[2] = 2;

  	   jcolM[0] = 0;
  	   jcolM[1] = 1;
    }

    return 0;
}

int matEqLink(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
#if 1
	if( id == 0 )
	{
       M[0] = 1.0;

       krowM[0] = 0;
       krowM[1] = 1;

       jcolM[0] = 0;
	}
	else
	{
	   M[0] = 1.0 + (double) id;

	   krowM[0] = 0;
	   krowM[1] = 1;

	   jcolM[0] = 1;
	}
#endif
    return 0;
}

};


int main(int argc, char ** argv) {

  MPI_Init(&argc, &argv);

  int nScenarios = 2;

  // set callbacks
  int nx0 = 2;
  int my0 = 2;
  int mz0 = 0;
  int myl0 = 1;
  int mzl0 = 0;

  FNNZ fnnzQ = &nnzAllZero;
  FNNZ fnnzA = &nnzMatEqStage1;
  FNNZ fnnzB = &nnzMatEqStage2;
  FNNZ fnnzBl = &nnzMatEqLink;
  FNNZ fnnzC = &nnzAllZero;
  FNNZ fnnzD = &nnzAllZero;
  FNNZ fnnzDl = &nnzAllZero;

  FVEC fc = &vecObj;
  FVEC fb = &vecEqRhs;
  FVEC fbl = &vecLinkRhs;
  FVEC fclow = &vecAllZero;
  FVEC fcupp = &vecAllZero;
  FVEC fxlow = &vecXlb;
  FVEC fxupp = &vecAllZero;
  FVEC ficlow = &vecAllZero;
  FVEC fixlow = &vecXlbActive;
  FVEC ficupp = &vecAllZero;
  FVEC fixupp = &vecAllZero;

  FVEC fdlupp = &vecAllZero;
  FVEC fdllow = &vecAllZero;
  FVEC fidlupp = &vecAllZero;
  FVEC fidllow = &vecAllZero;

  FMAT fQ = &matAllZero;
  FMAT fA = &matEqStage1;
  FMAT fB = &matEqStage2;
  FMAT fBl = &matEqLink;
  FMAT fC = &matAllZero;
  FMAT fD = &matAllZero;
  FMAT fDl = &matAllZero;

  ProbData probData(nScenarios);

#if LINKING_CONS
  //build the problem tree
  StochInputTree::StochInputNode dataLinkCons(&probData, 0,
				      nx0, my0, myl0, mz0, // mzl0
				      fQ, fnnzQ, fc,
				      fA, fnnzA,
				      fB, fnnzB,
					  fBl, fnnzBl,
				      fb, fbl,
				      fC, fnnzC,
				      fD, fnnzD,
					  //fDL, fnnzDL,
				      fclow, ficlow, fcupp, ficupp,
					  //fdllow, fidllow, fdlupp, fidlupp,
				      fxlow, fixlow, fxupp, fixupp, false );

  StochInputTree* root = new StochInputTree(dataLinkCons);
#else
  //build the problem tree
  StochInputTree::StochInputNode data(&probData, 0,
				      nx0,my0,mz0, //myl0, mzl0
				      fQ, fnnzQ, fc,
				      fA, fnnzA,
				      fB, fnnzB,
					  //fBL, fnnzBL,
				      fb, //fbl
				      fC, fnnzC,
				      fD, fnnzD,
					  //fDL, fnnzDL,
				      fclow, ficlow, fcupp, ficupp,
					  //fdllow, fidllow, fdlupp, fidlupp,
				      fxlow, fixlow, fxupp, fixupp, false );

  StochInputTree* root = new StochInputTree(data);
#endif


  for( int id = 1; id <= nScenarios; id++ ) {
	  int nx = 2;
	  int my = 2;
	  int mz = 0;
	  int myl = 1;
	  int mzl = 0;
#if LINKING_CONS
	  StochInputTree::StochInputNode dataLinkConsChild(&probData, id,
					nx, my, myl, mz, // mzl,
					fQ, fnnzQ, fc,
					fA, fnnzA,
					fB, fnnzB,
					fBl, fnnzBl,
					fb, fbl,
					fC, fnnzC,
					fD, fnnzD,
					//fDL, fnnzDL,
					fclow, ficlow, fcupp, ficupp,
					//fdllow, fidllow, fdlupp, fidlupp,
					fxlow, fixlow, fxupp, fixupp, false);

	  root->AddChild(new StochInputTree(dataLinkConsChild));
#else
	  StochInputTree::StochInputNode data(&probData, id,
					nx, my, mz, //myl, mzl
					fQ, fnnzQ, fc,
					fA, fnnzA,
					fB, fnnzB,
					//fBL, fnnzBL,
					fb, //fbl
					fC, fnnzC,
					fD, fnnzD,
					//fDL, fnnzDL,
					fclow, ficlow, fcupp, ficupp,
					//fdllow, fidllow, fdlupp, fidlupp,
					fxlow, fixlow, fxupp, fixupp, false);

	  root->AddChild(new StochInputTree(data));
#endif

  }

  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if( rank == 0 )
	  cout << "Using a total of " << size << " MPI processes." << endl;

  PIPSIpmInterface<sFactoryAug, MehrotraStochSolver> pipsIpm(root);

  if( rank == 0 )
	  cout << "PIPSIpmInterface created" << endl;

  if( rank == 0 )
	  cout << "solving..." << endl;

  pipsIpm.go();

  if( rank == 0 )
	  cout << "solving finished." << endl;


  // free memory
  delete root;

  MPI_Finalize();
  return 0;
}
