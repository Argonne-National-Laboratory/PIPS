#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "gmspipsio.h"
#include "gclgms.h"
#include "gdxcc.h"
#include "gmomcc.h"
#include "gevmcc.h"

#if defined(__cplusplus)
extern "C" {
#endif

//#define debug(msg)   printf("%s\n", msg);
#define debug(msg)   

#define GDXSAVECALLX(gxf,f) {                                                  \
                                int rc;                                        \
                                rc=f;                                          \
                                if (!rc)                                       \
                                {                                              \
                                  char s[GMS_SSSIZE];                          \
                                  gdxErrorStr(gxf, gdxGetLastError(gxf), s); \
                                  printf("*** Fatal GDX error (" __FILE__ " line: %d): %s\n",__LINE__,s);      \
                                  exit(1);                                     \
                                }                                              \
                            }
 
#define OFFSET 1
 
#define PRINTV1(n,v) \
if ( printLevel > 1 ) { int i; printf("%s:\n",#v);for (i=0; i<blk->n; i++) printf("  %3d: %15.3f\n", i, blk->v[i]); }
#define PRINTV4(n,i1,v1,i2,v2)                                                 \
if ( printLevel > 1 ) { int i; char s1[64], s2[64];                            \
  printf("%s,%s:\n",#v1,#v2);                                                  \
  for (i=0; i<blk->n; i++)                                                     \
  {                                                                            \
     sprintf(s1,"%15.3f",blk->v1[i]); sprintf(s2,"%15.3f",blk->v2[i]);         \
     printf("  %3d: %15s %15s\n", i, blk->i1[i]?s1:"-", blk->i2[i]?s2:"-");    \
  }                                                                            \
}                                    
#define PRINTMAT(mat,n)                                                        \
if ( printLevel > 1) { int i,j;                                                \
  printf("%s:\n",#mat);                                                        \
  for (i=0; i<blk->n; i++)                                                     \
     for (j=blk->rm##mat[i]; j<blk->rm##mat[i+1]; j++)                         \
        printf("  %3d %3d: %15.3f\n", i, blk->ci##mat[j], blk->val##mat[j]);   \
}                                    
 
int writeBlock(const char* scrFilename,  /** < scratch file name. If NULL write ASCII to stdout */
               GMSPIPSBlockData_t* blk,  /** < block structure to write */
               int printLevel)
{
   if (NULL==scrFilename)
   {
      printf("numBlocks: %5d\n", blk->numBlocks);
      printf("blockID:   %5d\n", blk->blockID);
      printf("n0:        %5d\n", blk->n0);
      printf("ni:        %5d\n", blk->ni);
      PRINTV1(ni,c);
      PRINTV4(ni,ixlow,xlow,ixupp,xupp);
      printf("mA:        %5d\n", blk->mA);
      if ( blk->mA )
      {
         PRINTV1(mA,b);
         PRINTMAT(A, mA);
         if (blk->blockID)
         {
            PRINTMAT(B, mA);
         }
      }
      printf("mC:        %5d\n", blk->mC);
      if ( blk->mC )
      {
         PRINTV4(mC,iclow,clow,icupp,cupp);
         PRINTMAT(C, mC);
         if (blk->blockID)
         {
            PRINTMAT(D, mC);
         }
      }
      printf("mBL:       %5d\n", blk->mBL);
      if ( blk->mBL )
      {
         PRINTV1(mBL,bL);
         PRINTMAT(BL, mBL);
      }
      printf("mDL:       %5d\n", blk->mDL);
      if ( blk->mDL )
      {
         PRINTV4(mBL,idlow,dlow,idupp,dupp);
         PRINTMAT(DL, mDL);
      }
   }

   return 0;
}               
 
int writeSolution(const char* gdxFilename,  /** < GDX output file name */
                  const char* dmFilename,   /** < DictMap file name */
                  GMSPIPSBlockData_t* blk)  /** < block structure to write */
{
   return 0;
}

void freeBlock(GMSPIPSBlockData_t* blk)
{
   if ( NULL==blk)
      return;
   
   if ( blk->ni )
   {
      free(blk->c);
      free(blk->xlow);
      free(blk->xupp);
      free(blk->ixlow);
      free(blk->ixupp);
   }
   if (blk->mA)
   {
      free(blk->b);
      free(blk->rmA);
      if (blk->nnzA)
      {
         free(blk->ciA);
         free(blk->valA);
      }
      if ( blk->blockID != 0 )
      {
         free(blk->rmB);      
         if (blk->nnzB)
         {
            free(blk->ciB);
            free(blk->valB);
         }
      }
   }
   if (blk->mC)
   {
      free(blk->clow);
      free(blk->cupp);
      free(blk->iclow);
      free(blk->icupp);
      free(blk->rmC);      
      if (blk->nnzC)
      {
         free(blk->ciC);
         free(blk->valC);
      }
      if ( blk->blockID != 0 )
      {
         free(blk->rmD);
         if (blk->nnzD)
         {
            free(blk->ciD);
            free(blk->valD);
         }
      }
   }
   if (blk->mBL)
   {
      free(blk->bL);
      free(blk->rmBL);
      if (blk->nnzBL)
      {
         free(blk->ciBL);
         free(blk->valBL);
      }
   }
   if (blk->mDL)
   {
      free(blk->dlow);
      free(blk->dupp);
      free(blk->idlow);
      free(blk->idupp);
      free(blk->rmDL);
      if (blk->nnzDL)
      {
         free(blk->ciDL);
         free(blk->valDL);
      }
   }   
}

#define PRINTANDEXIT(...) { printf (__VA_ARGS__); return 1; }
#define MATALLOC(mat)                                                  \
if ( blk->nnz##mat )                                                   \
{                                                                      \
   blk->ci##mat = (int32_t *) malloc(blk->nnz##mat * sizeof(int32_t)); \
   blk->val##mat = (double *) malloc(blk->nnz##mat * sizeof(double));  \
}

#if 1 
int doColumnPermutation(const gmoHandle_t gmo, const int strict, const int n, const int stageI, const int stage0, 
                     int32_t* n0, int32_t* ni, int perm[] )
{
   /* Permute the column variables */
   int j,rc;
   char colname[GMS_SSSIZE];

   double* stages = (double *) malloc(n*sizeof(double));
   rc = gmoGetVarScale(gmo,stages);
   if ( rc )
      PRINTANDEXIT("Problems accessing variable stages");

   for ( j=0; j<n; j++ )
   {
      if ( stageI == stages[j] )
      {
         (*ni)++;
      }
      else if ( stage0 == stages[j] )
      {
         (*n0)++;
      }
      else if ( strict )
         PRINTANDEXIT("Unmatched stage %d of column %s while reading stage %d", (int) stages[j], gmoGetVarNameOne(gmo,j,colname), stageI);
   }
   *ni = *n0;
   *n0 = 0;
   for ( j=0; j<n; j++ )
   {
      if ( stageI == stages[j] )
      {
         perm[(*ni)++] = j;
      }
      else if ( stage0 == stages[j] )
      {
         perm[(*n0)++] = j;
      }
   }
   
   return 0;   
}

int doRowPermutation(const gmoHandle_t gmo, const int strict, const int m, const int stageI, const int stageN, 
                     int32_t* meq, int32_t* mleq, int32_t* mLeq, int32_t* mLleq, int perm[] )
{
   /* Permute the row =,<=,linking =, linking <= */
   int i,rc;
   char rowname[GMS_SSSIZE];

   double* stages = (double *) malloc(m*sizeof(double));
   int* eTypes = (int *) malloc(m*sizeof(int));
   rc = gmoGetEquScale(gmo,stages);
   if ( rc )
      PRINTANDEXIT("Problems accessing equation stages");
   rc = gmoGetEquType(gmo, eTypes);
   if ( rc )
      PRINTANDEXIT("Problems accessing equation types");
   
   for ( i=0; i<m; i++ )
   {
      if ( gmoequ_N == eTypes[i] )
         continue;
      if ( eTypes[i] > gmoequ_N )
         PRINTANDEXIT("Unknown row type %d in row %s", eTypes[i], gmoGetEquNameOne(gmo,i,rowname));

      if ( stageI == stages[i] )
      {
         if ( gmoequ_E == eTypes[i] )
            (*meq)++;
         else
            (*mleq)++;
      }
      else if ( stageN == stages[i] )
      {
         if ( gmoequ_E == eTypes[i] )
            (*mLeq)++;
         else
            (*mLleq)++;
      }
      else if ( strict )
         PRINTANDEXIT("Unmatched stage %d of row %s while reading stage %d of %d stages", (int) stages[i], gmoGetEquNameOne(gmo,i,rowname), stageI, stageN);
   }
   *mLleq = *meq + *mleq + *mLeq;
   *mLeq = *meq + *mleq;
   *mleq = *meq;
   *meq = 0;

   for ( i=0; i<m; i++ )
   {
      if ( gmoequ_N == eTypes[i] )
         continue;
      if ( stageI == stages[i] )
      {
         if ( gmoequ_E == eTypes[i] )
            perm[(*meq)++] = i;
         else
            perm[(*mleq)++] = i;
      }
      else if ( stageN == stages[i] )
      {
         if ( gmoequ_E == eTypes[i] )
            perm[(*mLeq)++] = i;
         else
            perm[(*mLleq)++] = i;
      }
   }
   
   return 0;
}

int fillMatrix(const gmoHandle_t gmo, const int mStart, const int mEnd, int32_t rm[], int col[], double val[], int64_t* nnz)
{
   int m = mEnd - mStart;
   int i, nz, nlnz, n=gmoN(gmo);
   int* colidx = (int *) malloc(n*sizeof(int));
   double* jacval = (double *) malloc(n*sizeof(double));

   assert(m);
   rm[0] = 0;
   for ( i=0; i<m; i++)
   {
      if ( gmoGetRowSparse (gmo,i+mStart,colidx,jacval,NULL,&nz,&nlnz) )
         PRINTANDEXIT("Cannot get row %d in sparse format", i+mStart);
      assert(0==nlnz);
      rm[i+1] = rm[i] + nz;
      memcpy(col+rm[i],colidx,nz*sizeof(int));
      memcpy(val+rm[i],jacval,nz*sizeof(double));
   }
   *nnz = rm[m];   
   free(colidx);
   free(jacval);
   return 0;
} 

int readBlockSqueezed(int numBlocks,         /** < total number of blocks n in problem 0..n */
              int actBlock,                  /** < number of block to read 0..n */
              int strict,                    /** < indicator for clean blocks */
              const char* cntrFilename,      /** < Control filename */
              const char* blk0DictFilename,  /** < Dictionary file name for block 0 */
              const char* GAMSSysDir,        /** < GAMS system directory to locate shared libraries (can be NULL) */
              GMSPIPSBlockData_t* blk)       /** < block structure to be filled */
{
   gevHandle_t  fGEV=NULL;
   gmoHandle_t  fGMO=NULL;
   char msg[GMS_SSSIZE];
   int rc=0;

   assert(blk);
   assert(numBlocks>0);   
   assert(actBlock>=0 && actBlock<numBlocks);
   
   memset(blk,0,sizeof(GMSPIPSBlockData_t)); /* Initialize everything to 0/NULL */
   blk->numBlocks = numBlocks;
   blk->blockID = actBlock;   
   
   if ( GAMSSysDir )
   {
      rc = gevCreateD (&fGEV, GAMSSysDir, msg, sizeof(msg));
      if (!rc) 
         PRINTANDEXIT("Could not create gev object: %s\n", msg);
      rc = gmoCreateD (&fGMO, GAMSSysDir, msg, sizeof(msg));
      if (!rc) 
         PRINTANDEXIT("Could not create gmo object: %s\n", msg);
   }
   else
   {
      rc = gevCreate (&fGEV, msg, sizeof(msg));
      if (!rc)
         PRINTANDEXIT("Could not create gev object: %s\n", msg);
      rc = gmoCreate (&fGMO, msg, sizeof(msg));
      if (!rc) 
         PRINTANDEXIT("Could not create gmo object: %s\n", msg);
   }
   
   assert(cntrFilename);
   rc = gevInitEnvironmentLegacy(fGEV,cntrFilename);
   if (rc) 
      PRINTANDEXIT("Failed gevInitEnvironmentLegacy\n");

   rc = gmoRegisterEnvironment(fGMO, fGEV, msg);
   if (rc)
      PRINTANDEXIT("Could not register gev environment: %s\n", msg);
   
   rc = gmoLoadDataLegacy(fGMO, msg);
   if (rc)
      PRINTANDEXIT("Could not load data: %s\n", msg);

   /* Make sure objective variable and objective row are at the end */
   gmoIndexBaseSet(fGMO, 0);
   if ( gmoN(fGMO)!=gmoObjVar(fGMO)-1 )
      PRINTANDEXIT("Objective variable not at the of columns\n");
   
   if ( gmoM(fGMO)!=gmoObjRow(fGMO)-1 )
      PRINTANDEXIT("Objective equation not at the of rows\n");

   /* First we set some GMO objective function flavors */
   gmoObjStyleSet(fGMO, gmoObjType_Fun);
   gmoObjReformSet(fGMO, 1);
   
   if ( 0 == actBlock )
   {
      int mAStart, mAEnd;
      int mCStart, mCEnd;
      int mBLStart, mBLEnd;
      int mDLStart, mDLEnd;
      int gmom = gmoM(fGMO), gmon = gmoN(fGMO), n0, pipsm, pipsn;
      int* rowPerm;
      int* colPerm;
      
      rowPerm = (int *) malloc(gmom * sizeof(int)); 
      colPerm = (int *) malloc(gmon * sizeof(int)); 

      rc = doRowPermutation(fGMO, gmom, strict, actBlock+OFFSET, numBlocks+OFFSET, &(blk->mA), &(blk->mC), &(blk->mBL), &(blk->mDL), rowPerm);
      assert(0==rc);
      
      rc = doColumnPermutation(fGMO, gmon, strict, actBlock+OFFSET, OFFSET, &(blk->ni), &n0, colPerm);
      assert(0==rc);

      pipsn = blk->ni;
      pipsm = blk->mA + blk->mC + blk->mBL + blk->mDL;     

      rc = gmoSetRvEquPermutation(fGMO, rowPerm, pipsm);
      if (rc)
         PRINTANDEXIT("Could not install row permutation\n");
      free(rowPerm);

      rc = gmoSetRvVarPermutation(fGMO, colPerm, pipsn);
      if (rc)
         PRINTANDEXIT("Could not install column permutation\n");
      free(colPerm);
      
      /* Now fill the block */
      
      mAStart  = 0;      mAEnd  = mAStart  + blk->mA;
      mCStart  = mAEnd;  mCEnd  = mCStart  + blk->mC;
      mBLStart = mCEnd;  mBLEnd = mBLStart + blk->mBL;
      mDLStart = mBLEnd; mDLEnd = mDLStart + blk->mDL;

#if defined(GMSGENPIPSINPUTMAIN)
      /* Write row permutation for linking constraints to disk */
      {
         FILE* fpLinkingRowsPerm;
         int i;
         if ( (fpLinkingRowsPerm = fopen("block0LCPerm.txt", "w")) == NULL ) 
            PRINTANDEXIT("Error opening block0LCPerm.txt for writing\n");
         fprintf(fpLinkingRowsPerm,"%d %d\n",gmom,mDLEnd-mBLStart);
         for ( i=mBLStart; i<mDLEnd; i++ )
            fprintf(fpLinkingRowsPerm,"%d %d\n",i,rowPerm[i]);
         fclose(fpLinkingRowsPerm);
      }
      /* Write column permutation for linking variables to disk */
      {
         FILE* fpLinkingColumnsPerm;
         int j;
         if ( (fpLinkingColumnsPerm = fopen("block0LVPerm.txt", "w")) == NULL ) 
            PRINTANDEXIT("Error opening block0LVPerm.txt for writing\n");
         fprintf(fpLinkingColumnsPerm,"%d %d\n",gmon,pipsn);
         for ( j=0; j<pipsn; j++ )
            fprintf(fpLinkingColumnsPerm,"%d %d\n",j,colPerm[j]);
         fclose(fpLinkingColumnsPerm);
      }
#endif      

      /* Variable allocation */
      blk->c     = (double *)  calloc(blk->ni, sizeof(double)); 
      blk->xlow  = (double *)  calloc(blk->ni, sizeof(double)); 
      blk->xupp  = (double *)  calloc(blk->ni, sizeof(double)); 
      blk->ixlow = (int16_t *) calloc(blk->ni, sizeof(int16_t)); 
      blk->ixupp = (int16_t *) calloc(blk->ni, sizeof(int16_t));
   
      if (!gmoGetObjVector(fGMO,blk->c,NULL))
         PRINTANDEXIT("Could not get objective vector\n");
      if (!gmoGetVarLower(fGMO, blk->xlow))
         PRINTANDEXIT("Could not get variable lower bound\n");
      if (!gmoGetVarUpper(fGMO, blk->xupp))
         PRINTANDEXIT("Could not get variable lower bound\n");
      /* Set indicator vectors */
      {
         int j;
         for ( j=0; j<gmon; j++)
         {
            if (blk->xlow[j]!=GMS_SV_MINF)
               blk->ixlow[j] = 1;
            if (blk->xupp[j]!=GMS_SV_PINF)
               blk->ixupp[j] = 1;
         }
      }
      if ( blk->mA )
         blk->b     = (double *)  calloc(blk->mA,   sizeof(double)); 
      if ( blk->mC )
      {
         blk->clow  = (double *)  calloc(blk->mC,   sizeof(double)); 
         blk->cupp  = (double *)  calloc(blk->mC,   sizeof(double)); 
         blk->iclow = (int16_t *) calloc(blk->mC,   sizeof(int16_t)); 
         blk->icupp = (int16_t *) calloc(blk->mC,   sizeof(int16_t));
      }
      if ( blk->mBL )
         blk->bL    = (double *)  calloc(blk->mBL,  sizeof(double)); 
      if ( blk->mDL )
      {
         blk->dlow  = (double *)  calloc(blk->mDL,  sizeof(double)); 
         blk->dupp  = (double *)  calloc(blk->mDL,  sizeof(double)); 
         blk->idlow = (int16_t *) calloc(blk->mDL,  sizeof(int16_t)); 
         blk->idupp = (int16_t *) calloc(blk->mDL,  sizeof(int16_t));
      }

      {
         int i;
         double* rhs = (double *)  malloc(gmom*sizeof(double)); 
         int* eType = (int *)  malloc(gmom*sizeof(double)); 

         if (!gmoGetRhs(fGMO,rhs))
            PRINTANDEXIT("Could not get rhs vector\n");
         if (!gmoGetEquType(fGMO,eType))
            PRINTANDEXIT("Could not get equation type vector\n");
         
         if ( blk->mA )
            memcpy(blk->b,   rhs+mAStart ,blk->mA *sizeof(double));
         if ( blk->mC )
         {
            memcpy(blk->clow,rhs+mCStart ,blk->mC *sizeof(double));
            memcpy(blk->cupp,rhs+mCStart ,blk->mC *sizeof(double));
            memcpy(blk->bL,  rhs+mBLStart,blk->mBL*sizeof(double));
         }
         if ( blk->mDL )
         {
            memcpy(blk->dlow,rhs+mDLStart,blk->mDL*sizeof(double));
            memcpy(blk->dupp,rhs+mDLStart,blk->mDL*sizeof(double));
         }
         
         for ( i=mCStart; i<mCEnd; i++ )
         {
            assert(eType[i]==gmoequ_L || eType[i]==gmoequ_G);
            if ( gmoequ_L == eType[i] )
               blk->icupp[i-mCStart] = 1;
            else
               blk->iclow[i-mCStart] = 1;
         }
         for ( i=mDLStart; i<mDLEnd; i++ )
         {
            assert(eType[i]==gmoequ_L || eType[i]==gmoequ_G);
            if ( gmoequ_L == eType[i])
               blk->idupp[i-mDLStart] = 1;
            else
               blk->idlow[i-mDLStart] = 1;
         }
         free(rhs);
         free(eType);
      }

#define FILLMATRIX(mat)                                                                              \
         if ( blk->m##mat )                                                                          \
         {                                                                                           \
            blk->rm##mat  = (int32_t *) calloc(blk->m##mat+1, sizeof(int32_t));                      \
            if ( fillMatrix(fGMO,m##mat##Start,m##mat##End,blk->rm##mat,jcol,jval,&(blk->nnz##mat)) )\
               return 1;                                                                             \
            MATALLOC(mat);                                                                           \
            memcpy(blk->ci##mat,jcol,sizeof(int)*blk->nnz##mat);                                     \
            memcpy(blk->val##mat,jcol,sizeof(double)*blk->nnz##mat);                                 \
         }

      /* Fill matrices */
      {
         int gmonz = gmoNZ(fGMO);
         int* jcol = (int *) malloc(gmonz * sizeof(int));
         double* jval = (double *) malloc(gmonz * sizeof(double));

         FILLMATRIX(A);
         FILLMATRIX(C);
         FILLMATRIX(BL);
         FILLMATRIX(DL);

         free(jcol);
         free(jval);
      }
   }
      
   gmoFree(&fGMO);
   gevFree(&fGEV);
   return 0;
}              
#endif 

void copyGDXSymbol(int         numBlocks,
                   gdxHandle_t bGDX[],
                   gdxHandle_t fGDX,
                   const char* symName,
                   const int   nUelOffSet,
                   const int   offSet,
                   const int   stage[],
                   const int   vstage[],
                   const int   estage[],
                   const int   linkingBlock,
                   const int   readType)
{
   int k, rc, objVarUel;
   int symNr=0, symType=0, symDim=0, recNr=0, userInfo=0;
   int dimFirst=0;
   gdxValues_t   vals;
   gdxUelIndex_t keyInt;   
   char symText[GMS_SSSIZE];
 
   printf("Copying %s\n", symName); fflush(stdout);
   if ( 2 == readType )
   {
      GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "jobj", &symNr));
      GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &recNr));
      gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst);
      GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
      objVarUel = keyInt[0];
      assert(objVarUel);
   }
   rc = gdxFindSymbol(fGDX, symName, &symNr);
   if (!rc && 0==strcmp(symName,"ANl"))
   {
      for (k=0; k<numBlocks; k++)
      {
         GDXSAVECALLX(bGDX[k],gdxDataWriteRawStart(bGDX[k], symName, "Non-linear Jacobian indicator", 2, dt_par, 0));
         GDXSAVECALLX(bGDX[k],gdxDataWriteDone(bGDX[k]));
      }
      return;      
   }
   GDXSAVECALLX(fGDX,gdxSymbolInfo(fGDX, symNr, symText, &symDim, &symType));
   GDXSAVECALLX(fGDX,gdxSymbolInfoX(fGDX, symNr, &recNr, &userInfo, symText));

   for (k=0; k<numBlocks; k++)
   {
      GDXSAVECALLX(bGDX[k],gdxDataWriteRawStart(bGDX[k], symName, symText, symDim, symType, userInfo));
   }
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &recNr));
   while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
   {
      if ( 0 == readType || 3 == readType )
      {
         int blk;
         if ( 0 == readType )
            blk = stage[keyInt[0]-1];
         else
            blk = stage[keyInt[0]-1-nUelOffSet];
         if (blk==linkingBlock)
         {
            for (k=0; k<numBlocks; k++)
            {
               GDXSAVECALLX(bGDX[k],gdxDataWriteRaw (bGDX[k], keyInt, vals));      
            }
         }
         else
         {
            assert(blk<numBlocks);
            GDXSAVECALLX(bGDX[blk],gdxDataWriteRaw (bGDX[blk], keyInt, vals));
         }
      }
      else if ( 1 == readType )
      {
         for (k=0; k<numBlocks; k++)
         {
            GDXSAVECALLX(bGDX[k],gdxDataWriteRaw (bGDX[k], keyInt, vals));
         }
      }
      else if ( 2 == readType )
      {
         int row = keyInt[0]-1, col = keyInt[1]-nUelOffSet-1; 
         int jblk = vstage[col], iblk = estage[row];
         if (iblk<numBlocks)
         {
            if ( jblk!=0 && jblk!=iblk) /* Bad matrix element */
               printf("*** Unexpected matrix coefficient %f of equation e%d (stage=%d) and variable x%d (stage=%d)\n", vals[GMS_VAL_LEVEL], row+1,estage[row]+offSet,col+1,vstage[col]+offSet);
            else if (iblk > 0)
            {
               GDXSAVECALLX(bGDX[iblk],gdxDataWriteRaw (bGDX[iblk], keyInt, vals));      
            }
            else /* iblk == 0 */
            {
               assert(jblk<numBlocks);
               GDXSAVECALLX(bGDX[jblk],gdxDataWriteRaw (bGDX[jblk], keyInt, vals));      
            }
         }
         else /* linking constraint */
         {
            if ( keyInt[1] == objVarUel )
            {
               for (k=0; k<numBlocks; k++)
               {
                  GDXSAVECALLX(bGDX[k],gdxDataWriteRaw (bGDX[k], keyInt, vals));      
               }
            }
            else
            {
               assert(jblk<numBlocks);
               GDXSAVECALLX(bGDX[jblk],gdxDataWriteRaw (bGDX[jblk], keyInt, vals));      
            }
         }
      }
   }
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   for (k=0; k<numBlocks; k++)
   {
      GDXSAVECALLX(bGDX[k],gdxDataWriteDone(bGDX[k]));
   }
}

int gdxSplitting(const int numBlocks,        /** < total number of blocks n in problem 0..n */
              const int offset,              /** < indicator for clean blocks */
              const int skipStrings,         /** < indicator for not registering uels and strings */
              const char* gdxFilename,       /** < GDX file name with CONVERTD jacobian structure */
              const char* GAMSSysDir,        /** < GAMS system directory to locate shared libraries (can be NULL) */
	           const char* DictMap)           /** < GDX file name for dict structure (can be NULL) */
{
   gdxHandle_t  fGDX=NULL;
   gdxHandle_t* bGDX=NULL;
   gdxHandle_t  fDMGDX=NULL;
   gdxHandle_t* bDMGDX=NULL;
   
   char msg[GMS_SSSIZE];
   int rc=0;
   int symNr=0;
   int gdxN=0, gdxM=0;
   int dmM=0, dmN=0, dmTTLBLK=0;
   int dimFirst=0, numUels;
   gdxValues_t  vals;
   gdxUelIndex_t keyInt;
   int i=0,j=0,k=0;
   char bFileStem[GMS_SSSIZE], bDMFileStem[GMS_SSSIZE], fileName[GMS_SSSIZE];
   int* varstage = NULL;
   int* rowstage = NULL;


   assert(numBlocks>0);   
   assert(gdxFilename);
   
   if ( GAMSSysDir )
      rc = gdxCreateD (&fGDX, GAMSSysDir, msg, sizeof(msg));
   else
      rc = gdxCreate (&fGDX, msg, sizeof(msg));

   if ( !rc )
   {
      printf("Could not create gdx object: %s\n", msg);
      return 1;
   }

   gdxOpenRead(fGDX, gdxFilename, &rc);
   if (rc)
   {
      printf("Could not open GDX file %s (errNr=%d)\n", gdxFilename, rc);
      return 1;
   }

   if (DictMap)
   {
      if ( GAMSSysDir )
         rc = gdxCreateD (&fDMGDX, GAMSSysDir, msg, sizeof(msg));
      else
         rc = gdxCreate (&fDMGDX, msg, sizeof(msg));
      if ( !rc ) 
      {
         printf("Could not create gdx object: %s\n", msg);
         return 1;
      }
      gdxOpenRead(fDMGDX, DictMap, &rc);
      if (rc) 
      {
         printf("Could not open GDX file %s (errNr=%d)\n", DictMap, rc);
         return 1;
      }
   }
   
   printf("Reading equations stages\n");fflush(stdout);
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "i", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &gdxM));
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   rowstage = (int *)calloc(gdxM,sizeof(int));
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "e", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &rc));
   while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
      rowstage[keyInt[0]-1] = (int) vals[GMS_VAL_SCALE] - offset;
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));

   printf("Reading variable stages\n");fflush(stdout);
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "j", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &gdxN));
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   varstage = (int *) calloc(gdxN,sizeof(int));
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "x", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &rc));
   while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
      varstage[keyInt[0]-1-gdxM] = (int) vals[GMS_VAL_SCALE] - offset;
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));

   if (DictMap)
   {
     int nrRecs, k;
     GDXSAVECALLX(fDMGDX,gdxDataReadRawStart(fDMGDX, 1, &nrRecs));
     if (3!=nrRecs)
     {
        printf("expect 3 record in dictionary symbol 1, have %d\n", nrRecs);
        return 1;
     }
     gdxDataReadRaw(fDMGDX, keyInt, vals, &dimFirst); // ttblk
     dmTTLBLK = (int) vals[0];
     gdxDataReadRaw(fDMGDX, keyInt, vals, &dimFirst); // mincolcnt
     dmN = (int) vals[0];
     gdxDataReadRaw(fDMGDX, keyInt, vals, &dimFirst); // minrowcnt
     dmM = (int) vals[0];
     GDXSAVECALLX(fDMGDX,gdxDataReadDone(fDMGDX));
     if ((gdxN!=dmN) || (gdxM!=dmM))
     {
        printf("jacobian N=%d <> dict N=%d or jacobian M=%d <> dict M=%d\n", gdxN, dmN,gdxM,dmM);
        return 1;
     }
   }

   {
      char* lastdot;     
         strcpy(bFileStem,gdxFilename);
      lastdot = strrchr (bFileStem, '.');
      if (lastdot != NULL)
           *lastdot = '\0';

      if (DictMap)
      {
         strcpy(bDMFileStem,DictMap);
         lastdot = strrchr (bDMFileStem, '.');
         if (lastdot != NULL)
              *lastdot = '\0';
      }
   }

   GDXSAVECALLX(fGDX,gdxSystemInfo (fGDX, &symNr, &numUels));
   bGDX = (gdxHandle_t*) calloc(numBlocks, sizeof(gdxHandle_t));
   for (k=0; k<numBlocks; k++)
   {
      int nUel;
      if ( GAMSSysDir )
         rc = gdxCreateD (&(bGDX[k]), GAMSSysDir, msg, sizeof(msg));
      else
         rc = gdxCreate (&(bGDX[k]), msg, sizeof(msg));
      
      if ( !rc ) 
      {
         printf("Could not create %dth gdx object: %s\n", k, msg);
         return 1;
      }
      sprintf(fileName,"%s%d.gdx",bFileStem,k);
      gdxOpenWrite (bGDX[k], fileName, "gdxSplitter", &rc);
      assert(!rc);

      if ( !skipStrings )
      {
         printf("UEL Registration block %d\n",k);fflush(stdout);
         GDXSAVECALLX(bGDX[k],gdxUELRegisterRawStart(bGDX[k]));
         for (nUel=1; nUel<=numUels; nUel++)
         {
            char uel[GMS_SSSIZE];
            int map;
            GDXSAVECALLX(fGDX,gdxUMUelGet (fGDX,nUel,uel,&map));
            GDXSAVECALLX(bGDX[k],gdxUELRegisterRaw(bGDX[k],uel));
         }
         
         for (nUel=1; nUel<=numUels; nUel++)
         {
            char elemText[GMS_SSSIZE];
            int node, rc;
            rc = gdxGetElemText (fGDX, nUel, elemText, &node);
            if (0==rc) break;
            GDXSAVECALLX(bGDX[k],gdxAddSetText(bGDX[k], elemText, &node));
         }
         GDXSAVECALLX(bGDX[k],gdxUELRegisterDone(bGDX[k]));
      }
      else
      {
          GDXSAVECALLX(bGDX[k],gdxDataWriteRawStart(bGDX[k], "numUel", "Number of UELS", 0, GMS_DT_PAR, 0));
          vals[GMS_VAL_LEVEL] = numUels;
          GDXSAVECALLX(bGDX[k],gdxDataWriteRaw (bGDX[k], keyInt, vals));      
          GDXSAVECALLX(bGDX[k],gdxDataWriteDone(bGDX[k]));
      }
   }

   if (DictMap)
   {  
      int *nBlock, *mBlock;
      GDXSAVECALLX(fDMGDX,gdxSystemInfo (fDMGDX, &symNr, &numUels));
      bDMGDX = (gdxHandle_t*) calloc(numBlocks, sizeof(gdxHandle_t));
      nBlock = (int *) calloc(numBlocks, sizeof(int));
      mBlock = (int *) calloc(numBlocks, sizeof(int));
      for (j=0; j<gdxN; j++)
	      if (0==varstage[j])
	         for (k=0; k<numBlocks; k++)
	            nBlock[k]++;
	      else
	         nBlock[varstage[j]]++;
      for (i=0; i<gdxM; i++)
	      if (numBlocks==rowstage[i])
	         for (k=0; k<numBlocks; k++)
	            mBlock[k]++;
	      else
	         mBlock[rowstage[i]]++;
	
      for (k=0; k<numBlocks; k++)
      {
         int nUel;
         if ( GAMSSysDir )
            rc = gdxCreateD (&(bDMGDX[k]), GAMSSysDir, msg, sizeof(msg));
         else
            rc = gdxCreate (&(bDMGDX[k]), msg, sizeof(msg));
      
         if ( !rc ) 
         {
            printf("Could not create %dth gdx object: %s\n", k, msg);
            return 1;
         }
         sprintf(fileName,"%s%d.gdx",bDMFileStem,k);
         gdxOpenWrite (bDMGDX[k], fileName, "gdxSplitter", &rc);
         assert(!rc);

         printf("UEL Registration block %d\n",k);fflush(stdout);
         GDXSAVECALLX(bDMGDX[k],gdxUELRegisterRawStart(bDMGDX[k]));
         for (nUel=1; nUel<=numUels; nUel++)
         {
            char uel[GMS_SSSIZE];
            int map;
            GDXSAVECALLX(fDMGDX,gdxUMUelGet (fDMGDX,nUel,uel,&map));
            GDXSAVECALLX(bDMGDX[k],gdxUELRegisterRaw(bDMGDX[k],uel));
         }
         
         for (nUel=1; nUel<=numUels; nUel++)
         {
            char elemText[GMS_SSSIZE];
            int node, rc;
            rc = gdxGetElemText (fDMGDX, nUel, elemText, &node);
            if (0==rc) break;
            GDXSAVECALLX(bDMGDX[k],gdxAddSetText(bDMGDX[k], elemText, &node));
         }
         GDXSAVECALLX(bDMGDX[k],gdxUELRegisterDone(bDMGDX[k]));

         GDXSAVECALLX(bDMGDX[k],gdxDataWriteRawStart(bDMGDX[k], "model", "basic counts", 1, dt_par, 0));
         vals[GMS_VAL_LEVEL] = (double) dmTTLBLK; keyInt[0] = numUels-2;
         GDXSAVECALLX(bDMGDX[k],gdxDataWriteRaw (bDMGDX[k], keyInt, vals));      
         vals[GMS_VAL_LEVEL] = (double) nBlock[i]; keyInt[0] = numUels-1;
         GDXSAVECALLX(bDMGDX[k],gdxDataWriteRaw (bDMGDX[k], keyInt, vals));      
         vals[GMS_VAL_LEVEL] = (double) mBlock[i]; keyInt[0] = numUels;
         GDXSAVECALLX(bDMGDX[k],gdxDataWriteRaw (bDMGDX[k], keyInt, vals));      
         GDXSAVECALLX(bDMGDX[k],gdxDataWriteDone(bDMGDX[k]));
	 
         GDXSAVECALLX(bDMGDX[k],gdxDataWriteRawStart(bDMGDX[k], "display", "uels and their quotes", 1, dt_par, 0));
         GDXSAVECALLX(bDMGDX[k],gdxDataWriteDone(bDMGDX[k]));
	 
      }
   }
   
   /* Copy symbols */
   copyGDXSymbol(numBlocks,bGDX,fGDX,"i",      gdxM,offset,rowstage,NULL,     NULL,    numBlocks,0);
   copyGDXSymbol(numBlocks,bGDX,fGDX,"j",      gdxM,offset,varstage,NULL,     NULL,    0        ,3);
   copyGDXSymbol(numBlocks,bGDX,fGDX,"jobj",   gdxM,offset,NULL,    NULL,     NULL,    0        ,1);
   copyGDXSymbol(numBlocks,bGDX,fGDX,"objcoef",gdxM,offset,NULL,    NULL,     NULL,    0        ,1);
   copyGDXSymbol(numBlocks,bGDX,fGDX,"e",      gdxM,offset,rowstage,NULL,     NULL,    numBlocks,0);
   copyGDXSymbol(numBlocks,bGDX,fGDX,"x",      gdxM,offset,varstage,NULL,     NULL,    0        ,3);
   copyGDXSymbol(numBlocks,bGDX,fGDX,"A",      gdxM,offset,NULL,    varstage, rowstage,0        ,2);
   copyGDXSymbol(numBlocks,bGDX,fGDX,"ANl",    gdxM,offset,NULL,    NULL,     NULL,    0        ,1);
   
   for (k=0; k<numBlocks; k++)
   {
      int errNr = gdxGetLastError(bGDX[k]);
      if (errNr)
      {
         char s[GMS_SSSIZE];
         gdxErrorStr(bGDX[k], errNr, s);
         printf("GDX Error for GDX file %d: %s\n",k,s);
      }
      assert(0==errNr);
      gdxClose(bGDX[k]);
      gdxFree(&(bGDX[k]));
   }
   gdxClose(fGDX);
   gdxFree(&fGDX);

   if (DictMap)
   {
      int sym=2, cntSymRec=-1, *wroteSym;
      int symDim=0, symType=0, userInfo=0;
      char symText[GMS_SSSIZE], symName[GMS_SSSIZE];
      wroteSym = (int *) calloc(numBlocks, sizeof(int));

      for (i=0; i<gdxM; i++)
      {
	      if (cntSymRec <= 0)
         {
	         if (0==cntSymRec)
	         {
               GDXSAVECALLX(fDMGDX,gdxDataReadDone(fDMGDX));
	           for (k=0; k<numBlocks; k++)
	           {
	              if (wroteSym[k])
                    GDXSAVECALLX(bDMGDX[k],gdxDataWriteDone(bDMGDX[k]));
	              wroteSym[k] = 0;
	           }	    
	         }
	         while (0==cntSymRec)
	         {
	            sym++;
               GDXSAVECALLX(fDMGDX,gdxDataReadRawStart(fDMGDX, sym, &cntSymRec));
	            if (0==cntSymRec)
                  GDXSAVECALLX(fDMGDX,gdxDataReadDone(fDMGDX));
	         }
            GDXSAVECALLX(fDMGDX,gdxSymbolInfo(fDMGDX, sym, symName, &symDim, &symType));
            GDXSAVECALLX(fDMGDX,gdxSymbolInfoX(fDMGDX, sym, &cntSymRec, &userInfo, symText));
	      }
         gdxDataReadRaw(fDMGDX, keyInt, vals, &dimFirst); cntSymRec--;
	      if (numBlocks==rowstage[i])
	      {
	         for (k=0; k<numBlocks; k++)
	         {
	            if (!wroteSym[k])
                  GDXSAVECALLX(bDMGDX[k],gdxDataWriteRawStart(bDMGDX[k], symName, symText, symDim, dt_par, 0));
               GDXSAVECALLX(bDMGDX[k],gdxDataWriteRaw (bDMGDX[k], keyInt, vals));      
	         }
	      }
	      else
	      {
	         k = rowstage[i];
            if (!wroteSym[k])
               GDXSAVECALLX(bDMGDX[k],gdxDataWriteRawStart(bDMGDX[k], symName, symText, symDim, dt_par, 0));
            GDXSAVECALLX(bDMGDX[k],gdxDataWriteRaw (bDMGDX[k], keyInt, vals));      
	      }
      }
      assert(0==cntSymRec);
      for (j=0; j<gdxN; j++)
      {
	      if (cntSymRec <= 0)
         {
	         if (0==cntSymRec)
	         {
               GDXSAVECALLX(fDMGDX,gdxDataReadDone(fDMGDX));
	            for (k=0; k<numBlocks; k++)
	            {
	               if (wroteSym[k])
                     GDXSAVECALLX(bDMGDX[k],gdxDataWriteDone(bDMGDX[k]));
	               wroteSym[k] = 0;
	            }	    
	         }
	         while (0==cntSymRec)
	         {
	            sym++;
               GDXSAVECALLX(fDMGDX,gdxDataReadRawStart(fDMGDX, sym, &cntSymRec));
	            if (0==cntSymRec)
                  GDXSAVECALLX(fDMGDX,gdxDataReadDone(fDMGDX));
	         }
            GDXSAVECALLX(fDMGDX,gdxSymbolInfo(fDMGDX, sym, symName, &symDim, &symType));
            GDXSAVECALLX(fDMGDX,gdxSymbolInfoX(fDMGDX, sym, &cntSymRec, &userInfo, symText));
	      }
         gdxDataReadRaw(fDMGDX, keyInt, vals, &dimFirst); cntSymRec--;
	      if (0==varstage[j])
	      {
	         for (k=0; k<numBlocks; k++)
	         {
	            if (!wroteSym[k])
                  GDXSAVECALLX(bDMGDX[k],gdxDataWriteRawStart(bDMGDX[k], symName, symText, symDim, dt_par, 0));
               GDXSAVECALLX(bDMGDX[k],gdxDataWriteRaw (bDMGDX[k], keyInt, vals));      
	         }
	      }
	      else
	      {
	         k = varstage[j];
            if (!wroteSym[k])
               GDXSAVECALLX(bDMGDX[k],gdxDataWriteRawStart(bDMGDX[k], symName, symText, symDim, dt_par, 0));
            GDXSAVECALLX(bDMGDX[k],gdxDataWriteRaw (bDMGDX[k], keyInt, vals));      
	      }
      }
      GDXSAVECALLX(fDMGDX,gdxDataReadDone(fDMGDX));
      for (k=0; k<numBlocks; k++)
      {
         if (wroteSym[k])
            GDXSAVECALLX(bDMGDX[k],gdxDataWriteDone(bDMGDX[k]));
      }	    
      
      for (k=0; k<numBlocks; k++)
      {
         int errNr = gdxGetLastError(bDMGDX[k]);
         if (errNr)
         {
            char s[GMS_SSSIZE];
            gdxErrorStr(bDMGDX[k], errNr, s);
            printf("GDX Error for GDX file %d: %s\n",k,s);
         }
         assert(0==errNr);
         gdxClose(bDMGDX[k]);
         gdxFree(&(bDMGDX[k]));
      }
      gdxClose(fDMGDX);
      gdxFree(&fDMGDX);
   }
   return 0;   
}              


int readBlock(const int numBlocks,       /** < total number of blocks n in problem 0..n */
              const int actBlock,        /** < number of block to read 0..n */
              const int strict,          /** < indicator for clean blocks */
              const int offset,          /** < indicator for clean blocks */
              const char* gdxFilename,   /** < GDX file name with CONVERTD jacobian structure */
              const char* GAMSSysDir,    /** < GAMS system directory to locate shared libraries (can be NULL) */
              GMSPIPSBlockData_t* blk)   /** < block structure to be filled */
{
   gdxHandle_t  fGDX=NULL;
   char msg[GMS_SSSIZE];
   int rc=0;
   int symNr=0, idummy;
   int gdxN=0, gdxM=0, gdxNNZ=0;
   int dimFirst=0;
   gdxValues_t  vals;
   gdxUelIndex_t keyInt;
   int j, i;
   int* varPerm = NULL;
   int* equTypeNr = NULL;
   int objRowUel=0;
   int objVarUel=0;
   int* cIdxUel=NULL;
   double* cVal=NULL;
   int cCnt=0, badCnt=0, numUels=0;
   int objDirection=1;
   int objCoef=0.0;
   char** varname = NULL;
   char** rowname = NULL;
   int* varstage = NULL;
   int* rowstage = NULL;
   int* vemap = NULL;

   assert(blk);
   assert(numBlocks>0);   
   assert(actBlock>=0 && actBlock<numBlocks);
   assert(gdxFilename);
   
   if ( GAMSSysDir )
      rc = gdxCreateD (&fGDX, GAMSSysDir, msg, sizeof(msg));
   else
      rc = gdxCreate (&fGDX, msg, sizeof(msg));

   if ( !rc ) 
   {
      printf("Could not create gdx object: %s\n", msg);
      return 1;
   }
   
   gdxOpenRead(fGDX, gdxFilename, &rc);
   if (rc) 
   {
      printf("Could not open GDX file %s (errNr=%d)\n", gdxFilename, rc);
      return 1;
   }
   
   GDXSAVECALLX(fGDX,gdxSystemInfo (fGDX, &symNr, &numUels));
   if ( 0 == numUels )
   {
      GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "numUel", &symNr));
      GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &idummy));
      gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst);
      numUels = (int) vals[GMS_VAL_LEVEL];
      assert(numUels);
      GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   }
   vemap = (int*) calloc(numUels,sizeof(int));
   
   /* Objective variable UEL */
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "objcoef", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &idummy));
   gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst);
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   objDirection = (int) vals[GMS_VAL_LEVEL]; /* 1 for min, -1 for max */
   assert(objDirection==-1 || objDirection==1);

   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "i", &symNr));
   GDXSAVECALLX(fGDX,gdxSymbolInfoX (fGDX, symNr, &gdxM, &idummy, msg));
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "j", &symNr));
   GDXSAVECALLX(fGDX,gdxSymbolInfoX (fGDX, symNr, &gdxN, &idummy, msg));

   if ( strict )
   {
      int i=0,j=0;
      GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "j", &symNr));
      GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &idummy));
      varname = (char **) calloc(gdxN, sizeof(char*));
      while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
      {
         int node, textNr = (int) vals[GMS_VAL_LEVEL];
         char buf[GMS_SSSIZE];
         gdxGetElemText(fGDX, textNr, buf, &node);
         varname[j] = (char *)malloc(sizeof(char)*(strlen(buf)+1));
         strcpy(varname[j],buf);
         j++;
      }      
      GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));

      j = 0;
      GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "x", &symNr));
      GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &idummy));
      varstage = (int *) malloc(gdxN*sizeof(int));
      while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
         varstage[j++] = (int) vals[GMS_VAL_SCALE];
      GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));

      GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "i", &symNr));
      GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &idummy));
      rowname = (char **) calloc(gdxM, sizeof(char*));
      while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
      {
         int node, textNr = (int) vals[GMS_VAL_LEVEL];
         char buf[GMS_SSSIZE];
         gdxGetElemText(fGDX, textNr, buf, &node);
         rowname[i] = (char *)malloc(sizeof(char)*(strlen(buf)+1));
         strcpy(rowname[i],buf);
         i++;
      }      
      GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
      
      i = 0;
      GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "e", &symNr));
      GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &idummy));
      rowstage = (int *)malloc(gdxM*sizeof(int));
      while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
         rowstage[i++] = (int) vals[GMS_VAL_SCALE];
      GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   }
   
   /* Objective variable UEL */
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "jobj", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &idummy));
   gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst);
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   objVarUel = keyInt[0];
   assert(objVarUel);

   memset(blk,0,sizeof(GMSPIPSBlockData_t)); /* Initialize everything to 0/NULL */
   
   blk->numBlocks = numBlocks;
   blk->blockID = actBlock;
   
   debug("First pass over the variables to get variable counts right");
   /* First pass over the variables to get variable counts right */
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "j", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &idummy));
   {  int n = 1;
      while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
         vemap[keyInt[0]-1] = n++;
   }
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));

   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "x", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &idummy));
   varPerm = (int *) calloc(gdxN, sizeof(gdxN)); 
   while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
   {
      int blockNr = (int) vals[GMS_VAL_SCALE] - offset;
      int n = vemap[keyInt[0]-1]-1;
      if ( objVarUel==keyInt[0] ) /* skip objective variable */
         continue;
      if ( 0 == blockNr ) 
      {
         varPerm[n] = 1;
         blk->n0++;
      }
      else if ( actBlock == blockNr )
      {
         varPerm[n] = blockNr+1;
         blk->ni++;
      }
      else if ( strict > 1 )
         printf("*** Variable %s with block index %d while scanning for block index %d\n", varname[n],blockNr+offset, actBlock+offset);
   }
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));

   {
      int n0=0, ni=blk->n0;
      for ( j=0; j<gdxN; j++)
      {
         if (0==varPerm[j])
            continue;
         if (1==varPerm[j])
            varPerm[j] = ++n0;
         else
            varPerm[j] = ++ni;
      }
      assert(blk->n0==n0);
      assert(n0+blk->ni==ni);
   }
   if ( 0 == actBlock )
      blk->ni = blk->n0;

   if ( 0==blk->n0 || 0==blk->ni )
   {
      printf("Zero variable count: 0:%d %d:%d\n", blk->n0, actBlock, blk->ni);
      return 1;
   }

   /* Variable allocation */
   blk->c     = (double *)  calloc(blk->ni, sizeof(double)); 
   blk->xlow  = (double *)  calloc(blk->ni, sizeof(double)); 
   blk->xupp  = (double *)  calloc(blk->ni, sizeof(double)); 
   blk->ixlow = (int16_t *) calloc(blk->ni, sizeof(int16_t)); 
   blk->ixupp = (int16_t *) calloc(blk->ni, sizeof(int16_t));
   
   debug("Second pass over the variables to get the bounds");
   /* Second pass over the variables to get the bounds */
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &idummy));
   {  int n=0;
      while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
      {
         int blockNr = (int) vals[GMS_VAL_SCALE] - offset;
         if ( objVarUel==keyInt[0] ) /* skip objective variable */
            continue;
         if ( actBlock == blockNr )
         {
            if ( GMS_SV_MINF!=vals[GMS_VAL_LOWER] )
            {
               blk->xlow[n] = vals[GMS_VAL_LOWER];
               blk->ixlow[n] = 1;
            }
            if ( GMS_SV_PINF!=vals[GMS_VAL_UPPER] )
            {
               blk->xupp[n] = vals[GMS_VAL_UPPER];
               blk->ixupp[n] = 1;
            }
            n++;
            assert(n<=blk->ni);
         }
      }
   }
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   
   debug("First pass over the matrix to identify objective function");
   /* First pass over the matrix to identify objective function */
   cVal = (double *)malloc(gdxN*sizeof(double));
   cIdxUel = (int *)malloc(gdxN*sizeof(int));
   
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "A", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &gdxNNZ));
   while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
   {
      if ( 1==dimFirst && 0==objRowUel ) /* New row */
         cCnt = 0;
      if (objVarUel == keyInt[1])
      {
         if ( objRowUel )
         {
            printf("Objective variable used in more than one row: e%d e%d\n", objRowUel, keyInt[0]);
            return 1;           
         }
         objRowUel = keyInt[0];
         objCoef = vals[GMS_VAL_LEVEL];
      }
      else if ( 0==objRowUel || keyInt[0] == objRowUel )
      {
         cVal[cCnt] = vals[GMS_VAL_LEVEL];
         cIdxUel[cCnt] = keyInt[1];
         cCnt++;
      }
   }
   /* No objective coefficients */
   if (0==objRowUel) 
      cCnt = 0;
      
   assert(cCnt==0 || objCoef>1e-9 || objCoef<-1e-9);
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));

   debug("First pass over the equations to get equation counts right");
   /* First pass over the equations to get equation counts right */
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "i", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &idummy));
   {  int m = 1;
      while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
         vemap[keyInt[0]-1] = m++;
   }
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "e", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &idummy));
   equTypeNr = (int *) calloc(gdxM, sizeof(int)); 
   while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
   {
      int blockNr = (int) vals[GMS_VAL_SCALE] - offset;
      int m = vemap[keyInt[0]-1]-1;
      if ( GMS_SV_MINF==vals[GMS_VAL_LOWER] && GMS_SV_PINF==vals[GMS_VAL_UPPER] ) /* =n= */
         continue;
      if ( objRowUel==keyInt[0] ) /* skip objective defining row */
         continue;
      if ( actBlock == blockNr ) 
      {
         if ( GMS_SV_MINF==vals[GMS_VAL_LOWER] || GMS_SV_PINF==vals[GMS_VAL_UPPER] ) /* =l= or =g= */
            equTypeNr[m] = 2;
         else /* =e= */
            equTypeNr[m] = 1;
      }
      else if ( numBlocks == blockNr ) /* linking constraint */
      {
         if ( GMS_SV_MINF==vals[GMS_VAL_LOWER] || GMS_SV_PINF==vals[GMS_VAL_UPPER] ) /* =l= or =g= */
            equTypeNr[m] = 4;
         else /* =e= */
            equTypeNr[m] = 3;
      }
      else if ( strict > 1 )
         printf("*** Equation %s with block index %d while scanning for block index %d\n", rowname[m],blockNr+offset, actBlock+offset);
         
   }
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   
   for ( j=0; j<cCnt; j++ )
   {
      int col = vemap[cIdxUel[j]-1]-1;
      if ( 0==varPerm[col] )
      {
         if ( strict > 1)
            printf("*** Objective (%s) coefficient of variable %s with block index %d while scanning for block index %d\n", rowname[vemap[objRowUel-1]-1], varname[col], varstage[col], actBlock+offset);
         continue;
      }
      if ( varPerm[col] <= blk->n0 && actBlock > 0 )
         continue;

      col = varPerm[col] - ((varPerm[col] <= blk->n0)? 1:(blk->n0+1));
      blk->c[col] = objDirection*(-cVal[j]/objCoef);
   }
   
   for ( i=0; i<gdxM; i++)
   {
      switch (equTypeNr[i])
      {
         case 1: blk->mA++; break;
         case 2: blk->mC++; break;
         case 3: blk->mBL++; break;
         case 4: blk->mDL++; break;
      }
   }

   if ( blk->mA )
   {
      blk->b   = (double *)  calloc(blk->mA,   sizeof(double)); 
      blk->rmA = (int32_t *) calloc(blk->mA+1, sizeof(int32_t));
      if ( 0!=actBlock )
         blk->rmB = (int32_t *) calloc(blk->mA+1, sizeof(int32_t));
   }
   if ( blk->mC )
   {
      blk->clow  = (double *)  calloc(blk->mC,   sizeof(double)); 
      blk->cupp  = (double *)  calloc(blk->mC,   sizeof(double)); 
      blk->iclow = (int16_t *) calloc(blk->mC,   sizeof(int16_t)); 
      blk->icupp = (int16_t *) calloc(blk->mC,   sizeof(int16_t));
      blk->rmC   = (int32_t *) calloc(blk->mC+1, sizeof(int32_t));
      if ( 0!=actBlock )
         blk->rmD = (int32_t *) calloc(blk->mC+1, sizeof(int32_t));
   }      
   if ( blk->mBL )
   {
      blk->bL   = (double *)  calloc(blk->mBL,   sizeof(double)); 
      blk->rmBL = (int32_t *) calloc(blk->mBL+1, sizeof(int32_t));
   }

   if ( blk->mDL )
   {
      blk->dlow  = (double *)  calloc(blk->mDL,   sizeof(double)); 
      blk->dupp  = (double *)  calloc(blk->mDL,   sizeof(double)); 
      blk->idlow = (int16_t *) calloc(blk->mDL,   sizeof(int16_t)); 
      blk->idupp = (int16_t *) calloc(blk->mDL,   sizeof(int16_t));
      blk->rmDL  = (int32_t *) calloc(blk->mDL+1, sizeof(int32_t));
   }      
   
   debug("Second pass over the equations to get lhs/rhs");
   /* Second pass over the equations to get lhs/rhs */
   {
      int mA=0, mC=0, mBL=0, mDL=0;

      GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &idummy));
      while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
      {
         int blockNr = (int) vals[GMS_VAL_SCALE] - offset;
         if ( GMS_SV_MINF==vals[GMS_VAL_LOWER] && GMS_SV_PINF==vals[GMS_VAL_UPPER] ) /* =n= */
            continue;
         if ( objRowUel==keyInt[0] ) /* skip objective defining row */
            continue;
         if ( actBlock == blockNr ) 
         {
            if ( GMS_SV_MINF==vals[GMS_VAL_LOWER] ) /* =l= */
            {
               blk->cupp[mC] = vals[GMS_VAL_UPPER];
               blk->icupp[mC++] = 1;
            }
            else if ( GMS_SV_PINF==vals[GMS_VAL_UPPER] ) /* =g= */
            {
               blk->clow[mC] = vals[GMS_VAL_LOWER];
               blk->iclow[mC++] = 1;
            }
            else /* =e= */
               blk->b[mA++] = vals[GMS_VAL_LOWER];
         }
         else if ( numBlocks == blockNr ) /* linking constraint */
         {
            if ( GMS_SV_MINF==vals[GMS_VAL_LOWER] ) /* =l= */
            {
               blk->dupp[mDL] = vals[GMS_VAL_UPPER];
               blk->idupp[mDL++] = 1;
            }
            else if ( GMS_SV_PINF==vals[GMS_VAL_UPPER] ) /* =g= */
            {
               blk->dlow[mDL] = vals[GMS_VAL_LOWER];
               blk->idlow[mDL++] = 1;
            }
            else /* =e= */
               blk->bL[mBL++] = vals[GMS_VAL_LOWER];
         }
      }
      GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
      assert(mA==blk->mA);
      assert(mC==blk->mC);
      assert(mBL==blk->mBL);
      assert(mDL==blk->mDL);
   }
   
   /* For now */
   //assert(0==blk->mBL);
   //assert(0==blk->mDL);
   
   debug("Second pass over the matrix to get nnz counts right");
   /* Second pass over the matrix to get nnz counts right */
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "A", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &gdxNNZ));
   while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
   {
      int row = vemap[keyInt[0]-1]-1;
      int col = vemap[keyInt[1]-1]-1;

      //printf("actblk %d gdxN %d gdxM %d col %d coluel %d row %d rowuel %d\n", actBlock, gdxN, gdxM, keyInt[1], col, row, keyInt[0]);
      if ( objRowUel==keyInt[0] ) /* skip objective defining row */
         continue;
      if ( 0 == varPerm[col] && 0 == equTypeNr[row] ) /* skip block not relevant to actBlock */
         continue;
      if ( varPerm[col] <= blk->n0 && 0==equTypeNr[row] ) /* skip block not relevant to actBlock */
         continue;
      if ( equTypeNr[row] > 2 &&  0 == varPerm[col] ) /* skip nz in BL/DL not relevant to actBlock */
         continue;
      /* Skip the nz in BL/DL for block 0 variables if actBlock!=0 */
      if ( varPerm[col] <= blk->n0 && equTypeNr[row] > 2 && 0!=actBlock) 
         continue;
      if ( !(varPerm[col] && equTypeNr[row]) )
      {
         badCnt++;
         //printf("i=%d j=%d row=%d col=%d  varPerm[col]=%d blk->n0=%d equTypeNr[row]=%d\n",keyInt[0],keyInt[1],row,col, varPerm[col], blk->n0, equTypeNr[row]);
         if ( strict )
            printf("*** Unexpected matrix coefficient %f of equation %s (stage=%d) and variable %s (stage=%d) while scanning for block index %d [col %d coluel %d vp %d row %d rowuel %d et %d]\n", vals[GMS_VAL_LEVEL], rowname[row],rowstage[row], varname[col], varstage[col], actBlock+offset,col,keyInt[1],varPerm[col],row,keyInt[0],equTypeNr[row]);
      }
      
      if ( varPerm[col] <= blk->n0 )
      {
         switch (equTypeNr[row])
         {
            case 1: blk->nnzA++; break;
            case 2: blk->nnzC++; break;
            case 3: blk->nnzBL++; break;
            case 4: blk->nnzDL++; break;
         }
      }
      else
      {
         switch (equTypeNr[row])
         {
            case 1: blk->nnzB++; break;
            case 2: blk->nnzD++; break;
            case 3: blk->nnzBL++; break;
            case 4: blk->nnzDL++; break;
         }
      }
   }
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   assert( 0 == badCnt );
  
   MATALLOC(A);
   MATALLOC(B);
   MATALLOC(C);
   MATALLOC(D);
   MATALLOC(BL);
   MATALLOC(DL);
   
   debug("Third pass over the matrix to setup the matrix structures");
   //printf("mA %d mC %d mBL %d mDL %d\n", blk->mA, blk->mC, blk->mBL, blk->mDL);
   //printf("n0 %d mi %d\n", blk->n0, blk->ni);
   /* Third pass over the matrix to setup the matrix structures */
   {
      
      int mA=0, mC=0, mBL=0, mDL=0;
      int lastA=0, lastB=0, lastC=0, lastD=0, lastBL=0, lastDL=0;
      int irow=0, icol=0, lastrow=-1;
      
      GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &gdxNNZ));
      while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
      {
         int row = vemap[keyInt[0]-1]-1;
         int col = vemap[keyInt[1]-1]-1;
         int doX; /* 1:A 2:C 3:BL 4:DL 5:B 6:D*/
         
         //printf("uel %d %d val %f col %d varPerm %d row %d etype %d\n", keyInt[0], keyInt[1], vals[GMS_VAL_LEVEL], col, varPerm[col], row, equTypeNr[row] );
         if ( objRowUel==keyInt[0] ) /* skip objective defining row */
            continue;
         if ( 0 == varPerm[col] && 0 == equTypeNr[row] ) /* skip block not relevant to actBlock */
            continue;
         if ( varPerm[col] <= blk->n0 && 0==equTypeNr[row] ) /* skip block not relevant to actBlock */
            continue;
         
         if (1==dimFirst) /* new row */
         {
            int i;
            for (i=lastrow+1; i<row; i++)
            {
               if ( 3 == equTypeNr[i] )
                  ++mBL;
               else if ( 4 == equTypeNr[i] )
                  ++mDL;
            }
            lastrow = row;
            switch (equTypeNr[row])
            {
               case 1: irow = ++mA;  break;
               case 2: irow = ++mC;  break;
               case 3: irow = ++mBL; break;
               case 4: irow = ++mDL; break;
            }
         }
         if ( equTypeNr[row] > 2 &&  0 == varPerm[col] ) /* skip nz in BL/DL not relevant to actBlock */
            continue;
         /* Skip the nz in BL/DL for block 0 variables if actBlock!=0 */
         if ( varPerm[col] <= blk->n0 && equTypeNr[row] > 2 && 0!=actBlock) 
            continue;
         if ( varPerm[col] <= blk->n0 )
         {
            doX = equTypeNr[row];
            icol = varPerm[col]-1;
         }
         else
         {
            doX = (equTypeNr[row]<3)? equTypeNr[row]+4:equTypeNr[row];
            icol = varPerm[col] - blk->n0-1;
         }         
         
#define MATSTORE(mat)                                                                                \
{                                                                                                    \
   while (last##mat<irow) { blk->rm##mat[last##mat+1] = blk->rm##mat[last##mat]; last##mat++; }      \
   blk->ci##mat[blk->rm##mat[irow]] = icol; blk->val##mat[blk->rm##mat[irow]] = vals[GMS_VAL_LEVEL]; \
   blk->rm##mat[irow]++;                                                                             \
}
         //printf("dox %d irow %d icol %d\n", doX, irow, icol);fflush(stdout);
         switch (doX)
         {
            case 1: MATSTORE(A); break;
            case 2: MATSTORE(C); break;
            case 3: MATSTORE(BL); break;
            case 4: MATSTORE(DL); break;
            case 5: MATSTORE(B); break;
            case 6: MATSTORE(D); break;
         }
      }
      GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));

      
#define FILLMAT(mat,mmat)                                                                              \
if (blk->rm##mat)                                                                                      \
{                                                                                                      \
   while (last##mat<blk->m##mmat) {blk->rm##mat[last##mat+1] = blk->rm##mat[last##mat]; last##mat++; } \
   assert(blk->nnz##mat==blk->rm##mat[blk->m##mmat]);                                                  \
}

      debug("Finalize matrix structures");
      FILLMAT(A,A);
      FILLMAT(C,C);
      FILLMAT(BL,BL);
      FILLMAT(DL,DL);
      FILLMAT(B,A);
      FILLMAT(D,C);
   }

   gdxClose(fGDX);
   gdxFree(&fGDX);
   free(cVal);
   free(cIdxUel);
   free(varPerm);
   free(equTypeNr);
   free(vemap);

   if (strict)
   {
      int i,j;
      for (j=0; j<gdxN; j++)
         free(varname[j]);
      free(varname);
      free(varstage);
      for (i=0; i<gdxM; i++)
         free(rowname[i]);
      free(rowname);
      free(rowstage);
   }
   debug("Returning from readBlock");
   
   return 0;   

}

#if defined(__cplusplus)
}
#endif

