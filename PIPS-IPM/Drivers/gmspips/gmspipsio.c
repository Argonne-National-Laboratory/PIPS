#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "gmspipsio.h"
#include "gclgms.h"
#if defined(GDXSOURCE)
#include "gdxstatic.h"
#else
#include "gdxcc.h"
#include "gmomcc.h"
#include "gevmcc.h"
#endif   

#if defined(__cplusplus)
extern "C" {
#endif

//#define DEBUG(msg)   printf("%s\n", msg);
#define DEBUG(msg)   

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
         PRINTV4(mDL,idlow,dlow,idupp,dupp);
         PRINTMAT(DL, mDL);
      }
   }

   return 0;
}               
 
int writeSolution(const char* gdxFileStem,  /** < GDX file stem */
                  const int numcol,         /** < length of varl/varmlo/up array */
                  const int numErow,        /** < length of equEm array */
                  const int numIrow,        /** < length of equIl/equIm array */
                  const double objval,      /** < objective value */
                  double* varl,             /** < variable level (can be NULL) */
                  double* varmlo,           /** < variable marginals (can be NULL) */
                  double* varmup,           /** < variable marginals (can be NULL) */
                  double* equEl,            /** < equation =e= level (can be NULL) */
                  double* equIl,            /** < equation =lg= level (can be NULL) */
                  double* equEm,            /** < equation =e= marginals */
                  double* equIm,            /** < equation =lg= marginals */
                  const char* GAMSSysDir)   /** < GAMS system directory to locate shared libraries (can be NULL) */                  
{
   FILE *fmap;
   int* p2gvmap=NULL;
   double* p2gvlo=NULL;
   double* p2gvup=NULL;
   int* p2gemap=NULL;
   int* p2isE=NULL;
   char fileName[GMS_SSSIZE];
   int i, j, rc, p2gN, p2gM, objvar, objrow, gdxN, gdxM;
   gdxHandle_t fDCT=NULL;
   gdxHandle_t fSOL=NULL;
   char msg[GMS_SSSIZE];
   double* gvarl=NULL;
   double* gvarm=NULL;
   double* gequl=NULL;
   double* gequm=NULL;
   double objcoef = 1.0; /* TODO carry through the objective coefficient */
   gdxValues_t   vals, wvals;
   gdxUelIndex_t keyInt;   
   char symName[GMS_SSSIZE], symText[GMS_SSSIZE];
   int numrow, dimFirst, nrRecs, numUels, symDim, symType, userInfo, symCnt, symStart=3;
   int colsSeen=0, rowsSeen=0; 

#if !defined(GDXSOURCE)   
   if ( GAMSSysDir )
      rc = gdxCreateD (&fDCT, GAMSSysDir, msg, sizeof(msg));
   else
#endif      
      rc = gdxCreate (&fDCT, msg, sizeof(msg));
   if ( !rc ) 
   {
      printf("Could not create gdx object (dct): %s\n", msg);
      return -2;
   }
   strcpy(fileName,gdxFileStem);
   gdxOpenRead(fDCT, strcat(fileName,"_dict.gdx"), &rc);
   if (rc) 
   {
      printf("Could not open GDX file %s (errNr=%d)\n", fileName, rc);
      return -2;
   }

   GDXSAVECALLX(fDCT,gdxSystemInfo (fDCT, &symCnt, &numUels));
   numrow = numErow + numIrow;
   
   GDXSAVECALLX(fDCT,gdxDataReadRawStart(fDCT, 1, &nrRecs));
   if (3!=nrRecs)
   {
      printf("expect 3 record in dictionary symbol 1, have %d\n", nrRecs);
      return -2;
   }
   GDXSAVECALLX(fDCT,gdxDataReadRaw(fDCT, keyInt, vals, &dimFirst)); // ttblk
   GDXSAVECALLX(fDCT,gdxDataReadRaw(fDCT, keyInt, vals, &dimFirst)); // mincolcnt
   if (varl || varmlo || varmup)
   {
      if ( numcol+1!=(int) vals[0] )
      {
         printf("dictionary has %d variables PIPS only provides %d. Wrong dictionary file?\n", (int) vals[0], numcol+1);
         return -2;
      }
   }
   gdxN = (int) vals[0];
   
   GDXSAVECALLX(fDCT,gdxDataReadRaw(fDCT, keyInt, vals, &dimFirst)); // minrowcnt
   if (equEl || equEm || equIl || equIm)
   {
      if ( numErow+numIrow+1!=(int) vals[0] )
      {
         printf("dictionary has %d equations PIPS only provides %d. Wrong dictionary file?\n", (int) vals[0], numrow+1);
         return -2;
      }
   }
   gdxM = (int) vals[0];
   GDXSAVECALLX(fDCT,gdxDataReadDone(fDCT));
   
#if !defined(GDXSOURCE)   
   if ( GAMSSysDir )
      rc = gdxCreateD (&fSOL, GAMSSysDir, msg, sizeof(msg));
   else
#endif      
      rc = gdxCreate (&fSOL, msg, sizeof(msg));
   if ( !rc ) 
   {
      printf("Could not create gdx object (sol): %s\n", msg);
      return -2;
   }
   strcpy(fileName,gdxFileStem);
   gdxOpenWrite (fSOL, strcat(fileName,"_sol.gdx"), "GMSPIPS", &rc);
   if (rc) 
   {
      printf("Could not open GDX file %s for writing (errNr=%d)\n", fileName, rc);
      return -2;
   }

   strcpy(fileName,gdxFileStem);
   fmap = fopen(strcat(fileName,".map"), "r");
   if (!fmap)
      return -1;
   
   if(fscanf(fmap,"%d%d%d%d",&p2gN,&p2gM,&objvar,&objrow) != 4)
      return -1;

   if (varl || varmlo && varmup)
   {
      assert(p2gN==numcol);
      p2gvmap = (int *) malloc(numcol*sizeof(int));
      assert(p2gvmap);
      p2gvlo = (double *) malloc(numcol*sizeof(double));
      assert(p2gvlo);
      p2gvup = (double *) malloc(numcol*sizeof(double));
      assert(p2gvup);
   }
   if (equEl || equEm || equIl || equIm)
   {
      assert(p2gM==numrow);
      p2gemap = (int *) malloc(numrow*sizeof(int));
      assert(p2gemap);
      p2isE = (int *) malloc(numrow*sizeof(int));
      assert(p2isE);
   }
   
   for (j=0; j<numcol; j++)
   {
      int col;
	   double lo, up;
      if( fscanf(fmap,"%d %lf %lf",&col,&lo,&up) != 3);
         return -1;

      if (p2gvmap)
      {
         p2gvmap[j] = col;
         p2gvlo[j] = lo;
         p2gvup[j] = up;
	   }
   }
      
   for (i=0; i<numrow; i++)
   {
      int row, isE;
      if(fscanf(fmap,"%d %d",&row,&isE) != 2);
         return -1;
      
      if (p2gemap)
	   {
         p2gemap[i] = row;
         p2isE[i] = isE;
	   }
   }
   fclose(fmap);
   
   if (varl)
   {
      gvarl = (double *) malloc(gdxN*sizeof(double));
      assert(gvarl);
      gvarl[objvar] = objval; 
      for (j=0; j<numcol; j++)
         gvarl[p2gvmap[j]] = varl[j];
   }
   if (varl && varmlo && varmup)
   {
      gvarm = (double *) malloc(gdxN*sizeof(double));
      assert(gvarm);
      gvarm[objvar] = 0.0; 
      for (j=0; j<numcol; j++)
	  {
		  if (fabs(varl[j]-p2gvlo[j]) < 1e-6)
			  gvarm[p2gvmap[j]] = varmlo[j];
		  else if (fabs(varl[j]-p2gvup[j]) < 1e-6)
			  gvarm[p2gvmap[j]] = -varmup[j];
		  else
			  gvarm[p2gvmap[j]] = 0.0;
	  }
   }
   if (equEl && equIl)
   {
      int mE=0, mI=0;	   
      gequl = (double *) malloc(gdxM*sizeof(double));
      assert(gequl);
      gequl[objrow] = 0.0; 
      for (i=0; i<numrow; i++)
		  if (p2isE[i])
			  gequl[p2gemap[i]] = equEl[mE++];
		  else
			  gequl[p2gemap[i]] = equIl[mI++];
   }
   if (equEm && equIm)
   {
      int mE=0, mI=0;	   
      gequm = (double *) malloc(gdxM*sizeof(double));
      assert(gequm);
      gequm[objrow] = 1.0/objcoef;
      for (i=0; i<numrow; i++) 
	  {
		  if (p2isE[i])
			  gequm[p2gemap[i]] = equEm[mE++];
		  else
			  gequm[p2gemap[i]] = equIm[mI++];
	  }
   }
   if (p2gvmap) free(p2gvmap);
   if (p2gvlo)  free(p2gvlo);
   if (p2gvup)  free(p2gvup);
   if (p2gemap) free(p2gemap);
   if (p2isE)   free(p2isE);
   
   /* Initialize values to write to some useful defaults */  
   wvals[GMS_VAL_LEVEL] = 0.0;
   wvals[GMS_VAL_MARGINAL] = 0.0;
   wvals[GMS_VAL_LOWER] = GMS_SV_MINF;
   wvals[GMS_VAL_UPPER] = GMS_SV_PINF;
   wvals[GMS_VAL_SCALE] = 1.0;

   /* Now walk dictionary and variable and equation vector side-by-side */
   if (equEl && equIl || equEm && equIm)
   {
      for(; symStart<=symCnt; symStart++ )
      {
         if (rowsSeen==gdxM)
            break;
         GDXSAVECALLX(fDCT,gdxSymbolInfo(fDCT, symStart, symName, &symDim, &symType));
         GDXSAVECALLX(fDCT,gdxSymbolInfoX(fDCT, symStart, &nrRecs, &userInfo, symText));
         GDXSAVECALLX(fDCT,gdxDataReadRawStart(fDCT, symStart, &nrRecs));
         GDXSAVECALLX(fSOL,gdxDataWriteRawStart(fSOL, symName, symText, symDim, dt_equ, GMS_EQUEOFFSET+GMS_EQUTYPE_E));
         while ( gdxDataReadRaw(fDCT, keyInt, vals, &dimFirst) )
         {
            if (gequl) wvals[GMS_VAL_LEVEL] = gequl[rowsSeen];
            if (gequm) wvals[GMS_VAL_MARGINAL] = gequm[rowsSeen];
            rowsSeen++;
            GDXSAVECALLX(fSOL,gdxDataWriteRaw (fSOL, keyInt, wvals));      
         }
         GDXSAVECALLX(fDCT,gdxDataReadDone(fDCT));
         GDXSAVECALLX(fSOL,gdxDataWriteDone(fSOL));        
      }      
   }
   else /* fast forward to variable symbols */
   {  
      for(; symStart<=symCnt; symStart++ )
      {
         if (rowsSeen==gdxM)
            break;
         GDXSAVECALLX(fDCT,gdxSymbolInfo(fDCT, symStart, symName, &symDim, &symType));
         GDXSAVECALLX(fDCT,gdxSymbolInfoX(fDCT, symStart, &nrRecs, &userInfo, symText));
         rowsSeen += nrRecs;
      }
   }
   assert(rowsSeen==gdxM);   
   if (gequl) free(gequl);
   if (gequm) free(gequm);

   if (varl || varmlo && varmup)
   {
      for(; symStart<=symCnt; symStart++ )
      {
         GDXSAVECALLX(fDCT,gdxSymbolInfo(fDCT, symStart, symName, &symDim, &symType));
         GDXSAVECALLX(fDCT,gdxSymbolInfoX(fDCT, symStart, &nrRecs, &userInfo, symText));
         GDXSAVECALLX(fDCT,gdxDataReadRawStart(fDCT, symStart, &nrRecs));
         GDXSAVECALLX(fSOL,gdxDataWriteRawStart(fSOL, symName, symText, symDim, dt_var, 0));
         while ( gdxDataReadRaw(fDCT, keyInt, vals, &dimFirst) )
         {
            if (gvarl) wvals[GMS_VAL_LEVEL] = gvarl[colsSeen];
            if (gvarm) wvals[GMS_VAL_MARGINAL] = gvarm[colsSeen];
            colsSeen++;
            GDXSAVECALLX(fSOL,gdxDataWriteRaw (fSOL, keyInt, wvals));      
         }
         GDXSAVECALLX(fDCT,gdxDataReadDone(fDCT));
         GDXSAVECALLX(fSOL,gdxDataWriteDone(fSOL));        
      }
   }
   assert(colsSeen==gdxN);   
   if (gvarl) free(gvarl);
   if (gvarm) free(gvarm);
   
   /* Now regsiter uels */
   GDXSAVECALLX(fSOL,gdxUELRegisterRawStart(fSOL));
   for (i=1; i<=numUels; i++)
   {
      char uel[GMS_SSSIZE];
      int map;
      GDXSAVECALLX(fDCT,gdxUMUelGet (fDCT,i,uel,&map));
      GDXSAVECALLX(fSOL,gdxUELRegisterRaw(fSOL,uel));
   }
   GDXSAVECALLX(fSOL,gdxUELRegisterDone(fSOL));
   
   gdxClose(fDCT);
   gdxFree(&fDCT);
   gdxClose(fSOL);
   gdxFree(&fSOL);

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

#define MATALLOC(mat)                                                  \
if ( blk->nnz##mat )                                                   \
{                                                                      \
   blk->ci##mat = (int32_t *) malloc(blk->nnz##mat * sizeof(int32_t)); \
   blk->val##mat = (double *) malloc(blk->nnz##mat * sizeof(double));  \
}

void copyGDXSymbol(int         numBlocks,
                   int         actBlock,
                   gdxHandle_t bGDX[],
                   gdxHandle_t fGDX,
                   const char* symName,
                   const int   nUelOffSet,
                   const int   offSet,
                   const int   stage[],
                   const int   vstage[],
                   const int   estage[],
                   const int   linkingBlock,
                   const int   readType,
                   const int   objVarUel,
                   const int   objRowUel)
{
   int k, rc;
   int symNr=0, symType=0, symDim=0, recNr=0, userInfo=0;
   int dimFirst=0;
   gdxValues_t   vals;
   gdxUelIndex_t keyInt;   
   char symText[GMS_SSSIZE];
 
   rc = gdxFindSymbol(fGDX, symName, &symNr);
   if (!rc && 0==strcmp(symName,"ANl"))
   {
      printf("Copying %s\n", symName); fflush(stdout);
      for (k=(actBlock < 0)? 0:actBlock; k<numBlocks; k++)
      {
         GDXSAVECALLX(bGDX[k],gdxDataWriteRawStart(bGDX[k], symName, "Non-linear Jacobian indicator", 2, dt_par, 0));
         GDXSAVECALLX(bGDX[k],gdxDataWriteDone(bGDX[k]));
		 if (actBlock >= 0) break;
      }
      return;      
   }
   if (!rc && 0==strcmp(symName,"iobj"))
   {
      printf("Copying %s\n", symName); fflush(stdout);
      vals[GMS_VAL_LEVEL] = 0;
      keyInt[0] = objRowUel;
      for (k=(actBlock < 0)? 0:actBlock; k<numBlocks; k++)
      {
         GDXSAVECALLX(bGDX[k],gdxDataWriteRawStart(bGDX[k], symName, "Objective row (if reformulation possible, empty otherwise)", 1, dt_set, 0));
         GDXSAVECALLX(bGDX[k],gdxDataWriteRaw (bGDX[k], keyInt, vals));
         GDXSAVECALLX(bGDX[k],gdxDataWriteDone(bGDX[k]));
 		 if (actBlock >= 0) break;
      }
      return;      
   }
   GDXSAVECALLX(fGDX,gdxSymbolInfo(fGDX, symNr, symText, &symDim, &symType));
   GDXSAVECALLX(fGDX,gdxSymbolInfoX(fGDX, symNr, &recNr, &userInfo, symText));
   printf("Copying %s (#recs=%d)\n", symName, recNr); fflush(stdout);

   for (k=(actBlock < 0)? 0:actBlock; k<numBlocks; k++)
   {
      GDXSAVECALLX(bGDX[k],gdxDataWriteRawStart(bGDX[k], symName, symText, symDim, symType, userInfo));
      if (actBlock >= 0) break;
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
            for (k=(actBlock < 0)? 0:actBlock; k<numBlocks; k++)
            {
               GDXSAVECALLX(bGDX[k],gdxDataWriteRaw (bGDX[k], keyInt, vals));      
               if (actBlock >= 0) break;
            }
         }
         else
         {
            assert(blk<numBlocks);
			if (actBlock<0 || actBlock==blk)
			{
				GDXSAVECALLX(bGDX[blk],gdxDataWriteRaw (bGDX[blk], keyInt, vals));
			}
         }
      }
      else if ( 1 == readType )
      {
         for (k=(actBlock < 0)? 0:actBlock; k<numBlocks; k++)
         {
            GDXSAVECALLX(bGDX[k],gdxDataWriteRaw (bGDX[k], keyInt, vals));
            if (actBlock >= 0) break;
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
   			   if (actBlock<0 || actBlock==iblk)
			   {
				   GDXSAVECALLX(bGDX[iblk],gdxDataWriteRaw (bGDX[iblk], keyInt, vals));      
			   }
            }
            else /* iblk == 0 */
            {
               assert(jblk<numBlocks);
   			   if (actBlock<0 || actBlock==jblk)
			   {
				   GDXSAVECALLX(bGDX[jblk],gdxDataWriteRaw (bGDX[jblk], keyInt, vals));
			   }
            }
         }
         else /* linking constraint */
         {
            if ( keyInt[1] == objVarUel )
            {
               for (k=(actBlock < 0)? 0:actBlock; k<numBlocks; k++)
               {
                  GDXSAVECALLX(bGDX[k],gdxDataWriteRaw (bGDX[k], keyInt, vals));      
                  if (actBlock >= 0) break;
               }
            }
            else
            {
               assert(jblk<numBlocks);
   			   if (actBlock<0 || actBlock==jblk)
			   {
				   GDXSAVECALLX(bGDX[jblk],gdxDataWriteRaw (bGDX[jblk], keyInt, vals));      
			   }
            }
         }
      }
   }
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   for (k=(actBlock < 0)? 0:actBlock; k<numBlocks; k++)
   {
      GDXSAVECALLX(bGDX[k],gdxDataWriteDone(bGDX[k]));
      if (actBlock >= 0) break;
   }
}

int gdxSplitting(const int numBlocks,        /** < total number of blocks n in problem 0..n */
              const int actBlock,            /** < block to split from big GDX file, -1 split all */
              const int offset,              /** < indicator for clean blocks */
              const int skipStrings,         /** < indicator for not registering uels and strings */
              const char* gdxFilename,       /** < GDX file name with CONVERTD jacobian structure */
              const char* GAMSSysDir)        /** < GAMS system directory to locate shared libraries (can be NULL) */
{
   gdxHandle_t  fGDX=NULL;
   gdxHandle_t* bGDX=NULL;
   
   char msg[GMS_SSSIZE];
   int rc=0;
   int symNr=0, recNr=0;
   int gdxN=0, gdxM=0;
   int objVarUel=0, objRowUel=0;
   int dimFirst=0, numUels;
   gdxValues_t  vals;
   gdxUelIndex_t keyInt;
   int i=0,j=0,k=0;
   char bFileStem[GMS_SSSIZE], fileName[GMS_SSSIZE];
   double* varlo = NULL;
   double* varup = NULL;
   int* varstage = NULL;
   int* rowstage = NULL;
   int* isE = NULL;
   int intMaxBlock = 0;


   assert(numBlocks>0);   
   assert(gdxFilename);
   
#if !defined(GDXSOURCE)   
   if ( GAMSSysDir )
      rc = gdxCreateD (&fGDX, GAMSSysDir, msg, sizeof(msg));
   else
#endif
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
   
   printf("Reading equations stages\n");fflush(stdout);
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "i", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &gdxM));
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   rowstage = (int *)calloc(gdxM,sizeof(int));
   if (actBlock <= 0)
	   isE = (int *)calloc(gdxM,sizeof(int));
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "e", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &rc));
   while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
   {
      rowstage[keyInt[0]-1] = (int) vals[GMS_VAL_SCALE] - offset;
	  if ((actBlock <= 0) && (vals[GMS_VAL_LOWER] == vals[GMS_VAL_UPPER]))
		  isE[keyInt[0]-1] = 1;
   }
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));

   printf("Reading variable stages\n");fflush(stdout);
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "j", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &gdxN));
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   varstage = (int *) calloc(gdxN,sizeof(int));
   if (actBlock <= 0)
   {
	   varlo = (double *) malloc(gdxN*sizeof(double));
       varup = (double *) malloc(gdxN*sizeof(double));
       for (j=0; j<gdxN; j++)
       {
    	   varlo[j] = GMS_SV_MINF;
    	   varup[j] = GMS_SV_PINF;
       }
   }
   
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "x", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &rc));
   while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
   {
      int blk;
      blk = (int) vals[GMS_VAL_SCALE] - offset;
      if (actBlock <= 0)
	  {
		  varlo[keyInt[0]-1-gdxM] = vals[GMS_VAL_LOWER];
          varup[keyInt[0]-1-gdxM] = vals[GMS_VAL_UPPER];
	  }
      varstage[keyInt[0]-1-gdxM] = blk;
      if (blk>intMaxBlock)
         intMaxBlock = blk;
   }
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));

   if (intMaxBlock+offset != numBlocks)
   {
      printf("numblocks from command line (%d) does not match blocks in GDX file (%d)\n", numBlocks, intMaxBlock+offset);
      return 1;
   }
   {
      char* lastdot;     
         strcpy(bFileStem,gdxFilename);
      lastdot = strrchr (bFileStem, '.');
      if (lastdot != NULL)
           *lastdot = '\0';

   }

   GDXSAVECALLX(fGDX,gdxSystemInfo (fGDX, &symNr, &numUels));
   bGDX = (gdxHandle_t*) calloc(numBlocks, sizeof(gdxHandle_t));
   for (k=(actBlock < 0)? 0:actBlock;k<numBlocks; k++)
   {
      int nUel;
#if !defined(GDXSOURCE)   
      if ( GAMSSysDir )
         rc = gdxCreateD (&(bGDX[k]), GAMSSysDir, msg, sizeof(msg));
      else
#endif
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
         if (0==k)
            printf("#UELs: %d\n",numUels);fflush(stdout);

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
      if (actBlock >= 0)
		  break;
	  
   }

   /* Get objective uels */
   GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "jobj", &symNr));
   GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &recNr));
   gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst);
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   objVarUel = keyInt[0];
   assert(objVarUel);

   rc = gdxFindSymbol(fGDX, "iobj", &symNr);
   if (!rc) /* old jacobian without iobj */
   {
      GDXSAVECALLX(fGDX,gdxFindSymbol(fGDX, "A", &symNr));
      GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &recNr));
      while ( gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst) )
         if (objVarUel==keyInt[1])
         {  
            if (objRowUel)
            {
               objRowUel = 0;
               break;
            }
            else
            {
               objRowUel = keyInt[0];
            }
         }

      GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   }
   else
   {
      GDXSAVECALLX(fGDX,gdxDataReadRawStart(fGDX, symNr, &recNr));
      gdxDataReadRaw(fGDX, keyInt, vals, &dimFirst);
      GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
      objRowUel = keyInt[0];      
   }
   assert(objRowUel);
   
   /* Create PIPS2GAMS mapping file */
   if ( actBlock <= 0 )
   {
      FILE *fmap;
      int* p2gmap;
      int* p2gblkmap;
      int i,j,k,cnt,start;
      char fileName[GMS_SSSIZE];
      strcpy(fileName,bFileStem);
      fmap = fopen(strcat(fileName,".map"), "w");
      assert(fmap);
      fprintf(fmap,"%d %d %d %d\n", gdxN-1, gdxM-1, objVarUel-gdxM-1, objRowUel-1);

      p2gblkmap = (int*) calloc(numBlocks+1,sizeof(int));
      /* Counts by block */ 
      for (j=0; j<gdxN; j++)
         if (j!=objVarUel-gdxM-1)
            p2gblkmap[varstage[j]]++;
      /* Calculate start into map array by block */
      start = 0;
      for (k=0; k<numBlocks; k++)
      {
         cnt = p2gblkmap[k];
         p2gblkmap[k] = start;
         start += cnt;         
      }
      assert(start==gdxN-1);
      /* Fill map array */ 
      p2gmap = (int*) malloc(gdxN*sizeof(int));
      for (j=0; j<gdxN; j++)
         if (j!=objVarUel-gdxM-1)
            p2gmap[p2gblkmap[varstage[j]]++] = j;
      for (j=0; j<gdxN-1; j++)
         fprintf(fmap,"%d %g %g\n", p2gmap[j], varlo[p2gmap[j]], varup[p2gmap[j]]);
      free(p2gmap);
      
      /* Now the same for rows */
      memset(p2gblkmap,0,(numBlocks+1)*sizeof(int));
      /* Counts by block */ 
      for (i=0; i<gdxM; i++)
         if (i!=objRowUel-1)
            p2gblkmap[rowstage[i]]++;
      /* Calculate start into map array by block */
      start = 0;
      for (k=0; k<=numBlocks; k++)
      {
         cnt = p2gblkmap[k];
         p2gblkmap[k] = start;
         start += cnt;         
      }
      assert(start==gdxM-1);
      /* Fill map array */ 
      p2gmap = (int*) malloc(gdxM*sizeof(int));
      for (i=0; i<gdxM; i++)
         if (i!=objRowUel-1)
            p2gmap[p2gblkmap[rowstage[i]]++] = i;
      for (i=0; i<gdxM-1; i++)
         fprintf(fmap,"%d %d\n", p2gmap[i], isE[p2gmap[i]]);
      free(p2gmap);
      free(p2gblkmap);      
      fclose(fmap);
      free(varlo);
      free(varup);
      free(isE);
  }
   
   /* Copy symbols */
   copyGDXSymbol(numBlocks,actBlock,bGDX,fGDX,"i",      gdxM,offset,rowstage,NULL,     NULL,    numBlocks,0,objVarUel,objRowUel);
   copyGDXSymbol(numBlocks,actBlock,bGDX,fGDX,"j",      gdxM,offset,varstage,NULL,     NULL,    0        ,3,objVarUel,objRowUel);
   copyGDXSymbol(numBlocks,actBlock,bGDX,fGDX,"jobj",   gdxM,offset,NULL,    NULL,     NULL,    0        ,1,objVarUel,objRowUel);
   copyGDXSymbol(numBlocks,actBlock,bGDX,fGDX,"iobj",   gdxM,offset,NULL,    NULL,     NULL,    0        ,1,objVarUel,objRowUel);
   copyGDXSymbol(numBlocks,actBlock,bGDX,fGDX,"objcoef",gdxM,offset,NULL,    NULL,     NULL,    0        ,1,objVarUel,objRowUel);
   copyGDXSymbol(numBlocks,actBlock,bGDX,fGDX,"e",      gdxM,offset,rowstage,NULL,     NULL,    numBlocks,0,objVarUel,objRowUel);
   copyGDXSymbol(numBlocks,actBlock,bGDX,fGDX,"x",      gdxM,offset,varstage,NULL,     NULL,    0        ,3,objVarUel,objRowUel);
   copyGDXSymbol(numBlocks,actBlock,bGDX,fGDX,"A",      gdxM,offset,NULL,    varstage, rowstage,0        ,2,objVarUel,objRowUel);
   copyGDXSymbol(numBlocks,actBlock,bGDX,fGDX,"ANl",    gdxM,offset,NULL,    NULL,     NULL,    0        ,1,objVarUel,objRowUel);
  
   if (actBlock <= 0)
	   printf("gmspipscall: gmspips %d %s %s [scale] ...\n",numBlocks,bFileStem,GAMSSysDir);
   for (k=(actBlock < 0)? 0:actBlock;k<numBlocks; k++)
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
      if (actBlock >= 0)
		  break;
 
   }
   gdxClose(fGDX);
   gdxFree(&fGDX);

   if (varstage) free(varstage);
   if (rowstage) free(rowstage);

   return 0;   
}              


int readBlock(const int numBlocks,       /** < total number of blocks n in problem 0..n */
              const int actBlock,        /** < number of block to read 0..n */
              const int debugMode,       /** < indicator for clean blocks */
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
   double objCoef=0.0;
   char** varname = NULL;
   char** rowname = NULL;
   int* varstage = NULL;
   int* rowstage = NULL;
   int* vemap = NULL;
   int zjv=0;

   assert(blk);
   assert(numBlocks>0);   
   assert(actBlock>=0 && actBlock<numBlocks);
   assert(gdxFilename);
   
#if !defined(GDXSOURCE)   
   if ( GAMSSysDir )
      rc = gdxCreateD (&fGDX, GAMSSysDir, msg, sizeof(msg));
   else
#endif
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

   if ( debugMode )
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
   
   DEBUG("First pass over the variables to get variable counts right");
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
      else if ( debugMode > 1 )
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

   if ( 0==blk->ni )
	   if ( 0==actBlock )
       {
          printf("Zero joint variable count!\n");
		  zjv = 1;
       }
	   else
       {
          printf("Zero variable count for block %d!\n", actBlock);
          return 1;
       }

   if (!zjv)
   {
      /* Variable allocation */
      blk->c     = (double *)  calloc(blk->ni, sizeof(double)); 
      blk->xlow  = (double *)  calloc(blk->ni, sizeof(double)); 
      blk->xupp  = (double *)  calloc(blk->ni, sizeof(double)); 
      blk->ixlow = (int16_t *) calloc(blk->ni, sizeof(int16_t)); 
      blk->ixupp = (int16_t *) calloc(blk->ni, sizeof(int16_t));
      
      DEBUG("Second pass over the variables to get the bounds");
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
      
      DEBUG("First pass over the matrix to identify objective function");
      /* First pass over the matrix to identify objective function */
      cVal = (double *)malloc(gdxN*sizeof(double));
      cIdxUel = (int *)malloc(gdxN*sizeof(int));
   }
   
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
      else if ((0==objRowUel || keyInt[0] == objRowUel) && !zjv )
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

   DEBUG("First pass over the equations to get equation counts right");
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
      else if ( debugMode > 1 )
         printf("*** Equation %s with block index %d while scanning for block index %d\n", rowname[m],blockNr+offset, actBlock+offset);
         
   }
   GDXSAVECALLX(fGDX,gdxDataReadDone(fGDX));
   
   for ( j=0; j<cCnt; j++ )
   {
      int col = vemap[cIdxUel[j]-1]-1;
      if ( 0==varPerm[col] )
      {
         if ( debugMode > 1)
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
   
   DEBUG("Second pass over the equations to get lhs/rhs");
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
   
   if (zjv) goto skipMatrix;
   /* For now */
   //assert(0==blk->mBL);
   //assert(0==blk->mDL);
   
   DEBUG("Second pass over the matrix to get nnz counts right");
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
         if ( debugMode )
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
   if (badCnt)
      printf("*** %d unexpected matrix coefficient var.stage <> equ.stage and not linking. %s\n", badCnt, debugMode?"":"Add -d to gmschk call to see details.");
   assert( 0 == badCnt );
  
   MATALLOC(A);
   MATALLOC(B);
   MATALLOC(C);
   MATALLOC(D);
   MATALLOC(BL);
   MATALLOC(DL);
   
   DEBUG("Third pass over the matrix to setup the matrix structures");
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

      DEBUG("Finalize matrix structures");
      FILLMAT(A,A);
      FILLMAT(C,C);
      FILLMAT(BL,BL);
      FILLMAT(DL,DL);
      FILLMAT(B,A);
      FILLMAT(D,C);
   }

   free(cVal);
   free(cIdxUel);

   skipMatrix: 
   gdxClose(fGDX);
   gdxFree(&fGDX);
   free(varPerm);
   free(equTypeNr);
   free(vemap);

   if (debugMode)
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
   DEBUG("Returning from readBlock");
   
   return 0;   

}

int initGMSPIPSIO()
{
#if defined(GDXSOURCE)
   _P3_DllInit();   
#endif   
   return 0;
}

#if defined(__cplusplus)
}
#endif

