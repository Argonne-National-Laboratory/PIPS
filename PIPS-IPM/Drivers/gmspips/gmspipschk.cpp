#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "gmspipsio.h"
#include "gclgms.h"

void printUsage(void)
{
   printf("Usage: [-hsSptTxXdD] [-g GAMSSysDir] [-b actBlock] [-o n] numBlocks dirStem\n");
   printf("  -h        print usage\n");
   printf("  -s        strict mode\n");
   printf("  -S        very strict mode\n");
   printf("  -p        GNUPlot output to stdout, cmd= 'plot \"file\" using 2:1 with dots'\n");
   printf("  -t        split GDX file into multiple GDX files\n");
   printf("  -T        split GDX file into multiple GDX files without uels and strings\n");
   printf("  -w        output of block structure counts to stdout\n");
   printf("  -W        output of block structure to stdout\n");
   printf("  -x        dirStem is GDX file stem\n");
   printf("  -X        dirStem is GDX file\n");
   printf("  -d        dirStem is scratch directory stem (default)\n");
   printf("  -D        dirStem is scratch directory\n");
   printf("  -g        specify GAMS system directory\n");
   printf("  -b        specify single block\n");
   printf("  -o        specify stage offset (default 1)\n");
   printf("  numblocks total number of blocks\n");
   printf("  dirStem   scratch directory name stem\n");
}

#define SCRDIRSTEM   1
#define SCRDIR       2
#define GDXFILESTEM 10
#define GDXFILE     20

int  readOneBlock(const int numBlocks, 
                  const int actBlock, 
                  const int strict, 
                  const int offset, 
                  const int fType, 
                  const char* pDirStem,
                  const char* pGAMSSysDir, 
                  GMSPIPSBlockData_t* block)
{
   int rc;
   if ( SCRDIRSTEM == fType || SCRDIR == fType )
   {
      char fname[GMS_SSSIZE], dname[GMS_SSSIZE];
      
      if ( SCRDIRSTEM == fType )
      {
         sprintf(fname,"%s%d\\gamscntr.dat", pDirStem, actBlock);
         sprintf(dname,"%s%d\\gamsmat.dat", pDirStem, actBlock);
      }
      else
      {
         sprintf(fname,"%s\\gamscntr.dat", pDirStem);
         sprintf(dname,"%s\\gamsmat.dat", pDirStem);
      }
      rc = readBlockSqueezed(numBlocks,actBlock,strict,fname,dname,pGAMSSysDir,block);
   }
   else
   {
      char fname[GMS_SSSIZE];
      if ( GDXFILESTEM == fType )
         sprintf(fname,"%s%d.gdx", pDirStem, actBlock);
      else
         strcpy(fname,pDirStem);
      rc = readBlock(numBlocks,actBlock,strict,offset,fname,pGAMSSysDir,block);
   }
   return rc;
}   

int main(int argc, char* argv[])
{
   int numBlocks = 0;
   int actBlock = -1;
   int strict = 0;
   int offset = 1;
   int gnuplot = 0;
   int printMat = 0;
   int gdxSplit = 0;
   int fType = SCRDIRSTEM;
   char* pGAMSSysDir = NULL;
   char* pDirStem = NULL;
   
   while ( --argc > 0 )
   {
      char* p = *++argv;
      if ('-' == *p)
      {
         switch (*(p+1))
         {
            case 'g': pGAMSSysDir = *++argv; argc--; break;
            case 'b': actBlock = atoi(*++argv); argc--; break;
            case 'o': offset = atoi(*++argv); argc--; break;
            case 's': strict = 1; break;
            case 'S': strict = 2; break;
            case 't': gdxSplit = 1; break;
            case 'T': gdxSplit = 2; break;
            case 'p': gnuplot = 1; break;
            case 'w': printMat = 1; break;
            case 'W': printMat = 2; break;
            case 'd': fType = SCRDIRSTEM; break;
            case 'D': fType = SCRDIR; break;
            case 'x': fType = GDXFILESTEM; break;
            case 'X': fType = GDXFILE; break;
            case 'h': printUsage(); exit(0); break;
         }
      }
      else
         break;
   }
   if ( argc != 2 )
   {
      printUsage();
      exit(1);
   }
   numBlocks = atoi(*argv); 
   pDirStem = *++argv; 

   assert(pDirStem);
   assert(numBlocks>0);
   assert(actBlock<numBlocks);
   
   if ( gdxSplit && (fType == GDXFILE) )
   {
      int rc = gdxSplitting(numBlocks, offset, gdxSplit==2, pDirStem, pGAMSSysDir);
      assert(0==rc);
      return 0;
   }
   
   if ( actBlock>=0 )
   {
      GMSPIPSBlockData_t block;
      int rc;
      
      rc = readOneBlock(numBlocks, actBlock, strict, offset, fType, pDirStem, pGAMSSysDir, &block);
      assert(0==rc);
      if ( printMat )
      {
         rc = writeBlock(NULL,&block,printMat);
         assert(0==rc);
      }
      else
         printf("All done\n");

   }      
   else
   {
      int rc, nBlock0=-1;
      GMSPIPSBlockData_t** blocks = (GMSPIPSBlockData_t**) calloc(numBlocks,sizeof(GMSPIPSBlockData_t*));
      for (int blk=0; blk<numBlocks; blk++)
      {
         blocks[blk] = (GMSPIPSBlockData_t*) malloc(sizeof(GMSPIPSBlockData_t));
         rc = readOneBlock(numBlocks, blk, strict, offset, fType, pDirStem, pGAMSSysDir, blocks[blk]);
         assert(0==rc);
         if ( 0==blk )
            nBlock0 = blocks[blk]->n0;
         else
            assert(nBlock0 == blocks[blk]->n0);
      }
      
      if ( gnuplot )
      {
         /* Block 0 is the reference for all linking variables and constraints */
          int nSum = 0;
          int mSum = blocks[0]->mBL + blocks[0]->mDL;

          
#define plotmat(mat,matX,msum,nsum)                                                          \
             if (blocks[blk]->nnz ## mat)                                                    \
             {                                                                               \
                for (int i=0; i<blocks[blk]->m ## matX; i++)                                 \
                   for (int k=blocks[blk]->rm ## mat[i]; k<blocks[blk]->rm ## mat[i+1]; k++) \
                      printf("-%d %d\n", i+1+msum, blocks[blk]->ci ## mat[k]+1+nsum);        \
             }
             
          for (int blk=0; blk<numBlocks; blk++)                                              
          {
             //fprintf(stderr,"blk %d nSum=%d mSum=%d\n",blk,nSum,mSum);
             //printf("BL\n");
             plotmat(BL,BL,0, nSum);
             //printf("DL\n");
             plotmat(DL,DL,0, nSum);
             //printf("A %d\n",blocks[blk]->mA);
             plotmat(A,A,mSum,0);
             //printf("B\n");
             plotmat(B,A,mSum,nSum);
             //printf("C %d\n",blocks[blk]->mC);
             plotmat(C,C,mSum+blocks[blk]->mA,0);
             //printf("D\n");
             plotmat(D,C,mSum+blocks[blk]->mA,nSum);
             mSum += blocks[blk]->mA + blocks[blk]->mC;
             nSum += blocks[blk]->ni;
          }
          fprintf(stderr,"# All done\n");
      }
      else if ( printMat )
      {
         for (int blk=0; blk<numBlocks; blk++)
         {
            rc = writeBlock(NULL,blocks[blk],printMat);
            assert(0==rc);
         }
         printf("All done\n");
      }
      else 
         printf("All done\n");

      exit(0);
      for (int blk=0; blk<numBlocks; blk++)
      {
         freeBlock(blocks[blk]);
         free(blocks[blk]);
      }
      free(blocks);
   }      
   return 0;
}
