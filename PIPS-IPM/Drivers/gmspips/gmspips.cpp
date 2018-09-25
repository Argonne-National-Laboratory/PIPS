#if defined(GMS_PIPS)
#include "StochInputTree.h"
#include "PIPSIpmInterface.h"
//#include "MehrotraStochSolver.h"
#include "GondzioStochSolver.h"
#include "GondzioStochLpSolver.h"
#include "sFactoryAug.h"
#include "sFactoryAugSchurLeaf.h"
#endif
#if defined(GMS_MPI)
#include "mpi.h"
#endif

#include "gmspipsio.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <iostream>
#include <fstream>
using namespace std;

#define BICGSTAB

#if defined(GMS_PIPS)
extern int gOuterSolve;
extern int gInnerSCsolve;
#endif


extern "C" typedef int (*FNNZ)(void* user_data, int id, int* nnz);

/* Row-major format */
extern "C" typedef int (*FMAT)(void* user_data, int id, int* krowM, int* jcolM, double* M);

extern "C" typedef int (*FVEC)(void* user_data, int id, double* vec, int len);

extern "C" {

static int gmsRank=0;
static int size=1;
static bool allGDX=false;
static char fileName[256];
static char GDXDirectory[256];
static char* pGDXDirectory = NULL;
static int numBlocks=0;
FILE *fLog;

#define checkAndAlloc(blk)                                                        \
if (!blocks[blk])                                                                 \
{                                                                                 \
   int rc;                                                                        \
   fprintf(fLog,"Block %d read on gmsRank %d\n", blk, gmsRank);                         \
   blocks[blk] = (GMSPIPSBlockData_t*) malloc(sizeof(GMSPIPSBlockData_t));        \
   if ( !allGDX )                                                                 \
   {                                                                              \
      char fname[256];                                                            \
      sprintf(fname,"%s%d.gdx", fileName, blk);                                   \
      rc = readBlock(numBlocks,blk,0,1,fname,pGDXDirectory,blocks[blk]);          \
   }                                                                              \
   else                                                                           \
      rc = readBlock(numBlocks,blk,0,1,fileName,pGDXDirectory,blocks[blk]);       \
   assert(rc==0);                                                                 \
}

#define nCB(nType)                                                   \
int fsize##nType(void* user_data, int id, int* nnz)                  \
{                                                                    \
   GMSPIPSBlockData_t** blocks = (GMSPIPSBlockData_t**) user_data;   \
   checkAndAlloc(id);                                                \
   GMSPIPSBlockData_t* blk = blocks[id];                             \
   assert(blk);                                                      \
   *nnz = blk->nType;                                                \
   fprintf(fLog,"nCB blk=%d " #nType " %d\n",id,*nnz);               \
   return 0;                                                         \
}

nCB(ni)
nCB(mA)
nCB(mC)
nCB(mBL)
nCB(mDL)

#define nnzCB(nnzType)                                               \
int fnonzero##nnzType(void* user_data, int id, int* nnz)             \
{                                                                    \
   GMSPIPSBlockData_t** blocks = (GMSPIPSBlockData_t**) user_data;   \
   checkAndAlloc(id);                                                \
   GMSPIPSBlockData_t* blk = blocks[id];                             \
   assert(blk);                                                      \
   *nnz = blk->nnz##nnzType;                                         \
   fprintf(fLog,"nnzCB blk=%d " #nnzType " %d\n",id,*nnz);           \
   return 0;                                                         \
}

nnzCB(A)
nnzCB(B)
nnzCB(C)
nnzCB(D)
nnzCB(BL)
nnzCB(DL)

#define vecCB(vecType, size)                                         \
int fvec##vecType(void* user_data, int id, double* vec, int len)     \
{                                                                    \
   GMSPIPSBlockData_t** blocks = (GMSPIPSBlockData_t**) user_data;   \
   checkAndAlloc(id);                                                \
   GMSPIPSBlockData_t* blk = blocks[id];                             \
   assert(blk);                                                      \
   fprintf(fLog,"vecCB blk=%d " #vecType " len=%d allocLen %d\n",id,blk->size,len); \
   assert(len == blk->size);                                         \
   for( int i = 0; i < len; i++ ) {                                  \
      vec[i] = blk->vecType[i];                                      \
      fprintf(fLog,"  i=%d val=%g\n",i,vec[i]);                      \
   }                                                                 \
   return 0;                                                         \
}

vecCB(c    ,ni )
vecCB(xlow ,ni )
vecCB(ixlow,ni )
vecCB(xupp ,ni )
vecCB(ixupp,ni )
vecCB(b    ,mA )
vecCB(clow ,mC )
vecCB(iclow,mC )
vecCB(cupp ,mC )
vecCB(icupp,mC )
vecCB(bL   ,mBL)
vecCB(dlow ,mDL)
vecCB(idlow,mDL)
vecCB(dupp ,mDL)
vecCB(idupp,mDL)


#define matCB(mat,mmat)                                                   \
int fmat##mat(void* user_data, int id, int* krowM, int* jcolM, double* M) \
{                                                                         \
   GMSPIPSBlockData_t** blocks = (GMSPIPSBlockData_t**) user_data;        \
   checkAndAlloc(id);                                                     \
   GMSPIPSBlockData_t* blk = blocks[id];                                  \
   assert(blk);                                                           \
   fprintf(fLog,"matCB blk=%d " #mat " mLen %d nzLen %ld\n",id,blk->m##mmat,blk->nnz##mat); \
   if ( 0==blk->m##mmat )                                                 \
   {                                                                      \
     fprintf(fLog," empty\n"); \
     krowM[0] = 0;                                                        \
     return 0;                                                            \
   }                                                                      \
                                                                          \
   assert(blk->rm##mat);                                                  \
   for( int i = 0; i <= blk->m##mmat; i++ ) {                             \
      krowM[i] = blk->rm##mat[i];                                         \
      fprintf(fLog,"  i=%d krowM=%d\n",i,krowM[i]);                       \
   }                                                                      \
                                                                          \
   for( int k = 0; k < blk->nnz##mat; k++ ) {                             \
      jcolM[k] = blk->ci##mat[k];                                         \
      M[k] = blk->val##mat[k];                                            \
      fprintf(fLog,"  k=%d jcolM=%d M=%g\n",k,jcolM[k],M[k]);             \
   }                                                                      \
   return 0;                                                              \
}

matCB(A,A)
matCB(B,A)
matCB(C,C)
matCB(D,C)
matCB(BL,BL)
matCB(DL,DL)

int fnonzeroQ(void* user_data, int id, int* nnz)
{
	*nnz = 0;
	return 0;
}

int fmatQ(void* user_data, int id, int* krowM, int* jcolM, double* M)
{
   GMSPIPSBlockData_t* blk = ((GMSPIPSBlockData_t**) user_data)[id]; 
   assert(blk);
    
   for(int i = 0; i <= blk->ni; i++ )
       krowM[i] = 0;
   
   return 0;
}

}

static void setParams(ScalerType& scaler_type, bool& stepDiffLp, bool& presolve, bool& printsol, const char* paramname)
{
   if( strcmp(paramname, "scale") == 0 || strcmp(paramname, "scaleEqui") == 0 )
      scaler_type = SCALER_EQUI_STOCH;
   else if( strcmp(paramname, "scaleGeo") == 0 )
      scaler_type = SCALER_GEO_STOCH;
   else if( strcmp(paramname, "scaleGeoEqui") == 0 )
      scaler_type = SCALER_GEO_EQUI_STOCH;
   else if( strcmp(paramname, "stepLp") == 0 )
      stepDiffLp = true;
   else if( strcmp(paramname, "presolve") == 0 )
      presolve = true;
   else if( strcmp(paramname, "printsol") == 0 )
      printsol = true;
}

int main(int argc, char ** argv) 
{  

#if defined(GMS_MPI)
   MPI_Init(&argc, &argv);
   MPI_Barrier(MPI_COMM_WORLD);

   const double t0 = MPI_Wtime();
#endif

   GMSPIPSBlockData_t** blocks;
   ScalerType scaler_type = SCALER_NONE;
   bool stepDiffLp = false;
   bool presolve = false;
   bool printsol = false;

   if ( (argc<3) || (argc>8) )
   {
      cout << "Usage: " << argv[0] << " numBlocks all.gdx|blockstem [GDXLibDir] [scale] [stepLp] [presolve] [printsol]" << endl;
      exit(1);
   }
   
   allGDX = strstr(argv[2],".gdx")!=NULL;
   numBlocks = atoi(argv[1]);
   strcpy(fileName,argv[2]);
   if ( argc >= 4 )
   {
      strcpy(GDXDirectory,argv[3]);
      pGDXDirectory = &GDXDirectory[0];
   }
   
   for( int i = 5; i <= argc; i++ )
      setParams(scaler_type, stepDiffLp, presolve, printsol, argv[i - 1]);

   blocks = (GMSPIPSBlockData_t**) calloc(numBlocks,sizeof(GMSPIPSBlockData_t*));
#if 0
   int nBlock0 = 0;
   cout << "Start reading data from GDX and preparing the blocks" << endl;
   for (int blk=0; blk<numBlocks; blk++)
   {
      int rc;
      blocks[blk] = (GMSPIPSBlockData_t*) malloc(sizeof(GMSPIPSBlockData_t));
      if ( !allGDX )
      {
         char fname[256];
         sprintf(fname,"%s%d.gdx", argv[2], blk);
         rc = readBlock(numBlocks,blk,0,1,fname,(4==argc)? argv[3]:NULL,blocks[blk]);
      }
      else
         rc = readBlock(numBlocks,blk,0,1,argv[2],(4==argc)? argv[3]:NULL,blocks[blk]);
      if ( 0==blk )
         nBlock0 = blocks[blk]->n0;
      else
         assert(nBlock0 == blocks[blk]->n0);
      assert(0==rc);
   }
   cout << "Done reading data from GDX and preparing the blocks" << endl;
#endif

   FNNZ fsni   = &fsizeni   ;
   FNNZ fsmA   = &fsizemA   ;
   FNNZ fsmC   = &fsizemC   ;
#if defined(LINKCONSTR)   
   FNNZ fsmBL  = &fsizemBL  ;
   FNNZ fsmDL  = &fsizemDL  ;
#endif
   FNNZ fnnzQ  = &fnonzeroQ ;
   FNNZ fnnzA  = &fnonzeroA ;
   FNNZ fnnzB  = &fnonzeroB ;
   FNNZ fnnzC  = &fnonzeroC ;
   FNNZ fnnzD  = &fnonzeroD ;
#if defined(LINKCONSTR)   
   FNNZ fnnzBL = &fnonzeroBL;
   FNNZ fnnzDL = &fnonzeroDL;
#endif
   FVEC fc     = &fvecc    ;
   FVEC fxlow  = &fvecxlow ;
   FVEC fixlow = &fvecixlow;
   FVEC fxupp  = &fvecxupp ;
   FVEC fixupp = &fvecixupp;
   FVEC fb     = &fvecb    ;
   FVEC fclow  = &fvecclow ;
   FVEC ficlow = &fveciclow;
   FVEC fcupp  = &fveccupp ;
   FVEC ficupp = &fvecicupp;
#if defined(LINKCONSTR)   
   FVEC fbL    = &fvecbL   ;
   FVEC fdlow  = &fvecdlow ;
   FVEC fidlow = &fvecidlow;
   FVEC fdupp  = &fvecdupp ;
   FVEC fidupp = &fvecidupp;
#endif
   
   FMAT fA  = &fmatA ;
   FMAT fB  = &fmatB ;
   FMAT fC  = &fmatC ;
   FMAT fD  = &fmatD ;
#if defined(LINKCONSTR)   
   FMAT fBL = &fmatBL;
   FMAT fDL = &fmatDL;
#endif
   FMAT fQ  = &fmatQ ;

#if defined(GMS_PIPS)
   //build the problem tree
   StochInputTree::StochInputNode data(blocks, 0,
#if defined(LINKCONSTR)   
      fsni, fsmA, fsmBL, fsmC, fsmDL,
#else
      fsni, fsmA, fsmC,
#endif      
      fQ,  fnnzQ, fc,
      fA,  fnnzA,
      fB,  fnnzB,
#if defined(LINKCONSTR)   
      fBL, fnnzBL,
#endif      
      fb,  
#if defined(LINKCONSTR)   
      fbL,
#endif      
      fC,  fnnzC,
      fD,  fnnzD,
#if defined(LINKCONSTR)   
      fDL, fnnzDL,
#endif      
      fclow, ficlow, fcupp, ficupp,
#if defined(LINKCONSTR)   
      fdlow, fidlow, fdupp, fidupp,
#endif      
      fxlow, fixlow, fxupp, fixupp, false );
   StochInputTree* root = new StochInputTree(data);
#endif
   for( int blk = 1; blk < numBlocks; blk++ ) 
   {

#if defined(GMS_PIPS)
      StochInputTree::StochInputNode data(blocks, blk,
#if defined(LINKCONSTR)
         fsni, fsmA, fsmBL, fsmC, fsmDL,
#else
         fsni, fsmA, fsmC,
#endif
         fQ,  fnnzQ, fc,
         fA,  fnnzA,
         fB,  fnnzB,
#if defined(LINKCONSTR)   
         fBL, fnnzBL,
#endif      
         fb,  
#if defined(LINKCONSTR)   
         fbL,
#endif      
         fC,  fnnzC,
         fD,  fnnzD,
#if defined(LINKCONSTR)   
         fDL, fnnzDL,
#endif      
         fclow, ficlow, fcupp, ficupp,
#if defined(LINKCONSTR)   
         fdlow, fidlow, fdupp, fidupp,
#endif      
         fxlow, fixlow, fxupp, fixupp, false );

       root->AddChild(new StochInputTree(data));
#endif       
   }

#if defined(GMS_MPI)   
   MPI_Comm_rank(MPI_COMM_WORLD, &gmsRank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   char fbuf[256];
#if defined (GMS_LOG)
   sprintf(fbuf,"log%d.txt", gmsRank);
#else
   sprintf(fbuf,"/dev/null");
#endif
   fLog = fopen(fbuf, "w+");
   fprintf(fLog, "PIPS Log for gmsRank %d\n", gmsRank);
#endif  
   if( gmsRank == 0 )
#if defined(LINKCONSTR)   
      cout << "Using version with linking constraint." << endl;
#else
      cout << "Using version without linking constraint." << endl;
#endif      
   if( gmsRank == 0 )
      cout << "Using a total of " << size << " MPI processes." << endl;

#if defined(GMS_PIPS)
#ifdef BICGSTAB
   if( gmsRank == 0 )
      cout << "using outer BICGSTAB" << endl;
   gOuterSolve=2;
#else
   gOuterSolve=0;
#endif

#ifdef INNER_BICGSTAB
   gInnerSCsolve=2;

   if( gmsRank == 0 )
      cout << "using inner BICGSTAB" << endl;
#else
   gInnerSCsolve=0;
#endif

   std::vector<double> primalSolVec;
   double objective = 0.0;

	if (stepDiffLp)
	{
	   if( gmsRank == 0 )
	      cout << "Different steplengths in primal and dual direction are used." << endl;

#ifdef WITH_PARDISO
      PIPSIpmInterface<sFactoryAugSchurLeaf, GondzioStochLpSolver> pipsIpm(root, MPI_COMM_WORLD,
            scaler_type, presolve ? PRESOLVER_STOCH : PRESOLVER_NONE );
#else
      PIPSIpmInterface<sFactoryAug, GondzioStochLpSolver> pipsIpm(root, MPI_COMM_WORLD,
            scaler_type, presolve ? PRESOLVER_STOCH : PRESOLVER_NONE );
#endif

		if( gmsRank == 0 )
		   cout << "PIPSIpmInterface created" << endl;

		if( gmsRank == 0 )
		   cout << "solving..." << endl;

		pipsIpm.go();
      if( printsol )
      {
         primalSolVec = pipsIpm.gatherPrimalSolution();
         objective = pipsIpm.getObjective();
      }
	}

	else {
#ifdef WITH_PARDISO
      PIPSIpmInterface<sFactoryAugSchurLeaf, GondzioStochSolver> pipsIpm(root, MPI_COMM_WORLD,
            scaler_type, presolve ? PRESOLVER_STOCH : PRESOLVER_NONE );
#else
      PIPSIpmInterface<sFactoryAug, GondzioStochSolver> pipsIpm(root, MPI_COMM_WORLD,
            scaler_type, presolve ? PRESOLVER_STOCH : PRESOLVER_NONE );
#endif

		if( gmsRank == 0 )
		   cout << "PIPSIpmInterface created" << endl;

		if( gmsRank == 0 )
		   cout << "solving..." << endl;

		pipsIpm.go();
      if( printsol )
      {
         primalSolVec = pipsIpm.gatherPrimalSolution();
         objective = pipsIpm.getObjective();
      }
	}
   if( gmsRank == 0 )
      cout << "solving finished." << endl;

   if( printsol && gmsRank == 0 )
   {
      double* varl;
      int rc;
      
      varl = (double*) malloc(primalSolVec.size()*sizeof(double));
      for( size_t i = 0; i < primalSolVec.size(); i++ )
         varl[i] = primalSolVec[i];

      rc = writeSolution(fileName,primalSolVec.size(),0,objective,varl,NULL,NULL,NULL,pGDXDirectory);
      
      free(varl);
      if (0==rc)
         std::cout << "Solution written to " << fileName << "_sol.gdx" << std::endl;
      else if (-1==rc)
         std::cout << "Could not access " << fileName << ".map" << std::endl;
      else
         std::cout << "Other error writing solution: rc=" << rc << std::endl;
   }
   
   // free memory
  delete root;

#endif


  for (int blk=0; blk<numBlocks; blk++)
  {
     freeBlock(blocks[blk]);
     free(blocks[blk]);
  }
  free(blocks);
  
#if defined(GMS_MPI)   
  MPI_Barrier(MPI_COMM_WORLD);
  const double t1 = MPI_Wtime();

  if( gmsRank == 0 )
     std::cout << "---total time (in sec.): " << t1 - t0 << std::endl;

  MPI_Finalize();
#endif 
  fclose(fLog);
 
  return 0;
}
