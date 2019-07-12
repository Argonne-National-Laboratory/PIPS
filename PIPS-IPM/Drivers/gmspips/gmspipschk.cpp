#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "gmspipsio.h"
#include "gclgms.h"

void printUsage(void)
{
   printf("Usage: [-dhtTwWx] [-b actBlock] [-g GAMSSysDir] [-o n] numBlocks file[Stem]\n");
   printf("  -h          print usage\n\n");
   printf("  Splitting operation:\n");
   printf("    -t        split GDX file into multiple GDX files\n");
   printf("    -T        split GDX file into multiple GDX files without uels and strings\n");
   printf("    -g        specify GAMS system directory\n");
   printf("    -b        specify single block\n");
   printf("    -o        specify stage offset (default 1)\n");
   printf("    numblocks total number of blocks\n");
   printf("    file      GDX file\n\n");
   printf("  Analysis operation:\n");
   printf("    -d        debugging mode (unmatched vars, equs, and matrix elements with good names)\n");
   printf("    -w        output of block structure counts to stdout\n");
   printf("    -W        output of block structure to stdout\n");
   printf("    -x        fileStem is GDX file stem\n");
   printf("    -g        specify GAMS system directory\n");
   printf("    -b        specify single block\n");
   printf("    -o        specify stage offset (default 1)\n");
   printf("    numblocks total number of blocks\n");
   printf("    fileStem  GDX file or file stem\n");
}

#define GDXFILESTEM 10
#define GDXFILE     20

int  readOneBlock(const int numBlocks, 
                  const int actBlock, 
                  const int debug, 
                  const int offset, 
                  const int fType, 
                  const char* pDirStem,
                  const char* pGAMSSysDir, 
                  GMSPIPSBlockData_t* block)
{
   char fname[GMS_SSSIZE];
   if ( GDXFILESTEM == fType )
      sprintf(fname,"%s%d.gdx", pDirStem, actBlock);
   else
      strcpy(fname,pDirStem);
   return readBlock(numBlocks,actBlock,debug,offset,fname,pGAMSSysDir,block);
}   

int main(int argc, char* argv[])
{
   int numBlocks = 0;
   int actBlock = -1;
   int debug = 0;
   int offset = 1;
   // int gnuplot = 0;
   int printMat = 0;
   int gdxSplit = 0;
   int rc = 0;
   int fType = GDXFILE;
   char* pGAMSSysDir = NULL;
   char* pFileStem = NULL;
   
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
            case 'd': debug = 1; break;
            case 't': gdxSplit = 1; break;
            case 'T': gdxSplit = 2; break;
            case 'w': printMat = 1; break;
            case 'W': printMat = 2; break;
            case 'x': fType = GDXFILESTEM; break;
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
   if (debug && fType == GDXFILESTEM)
	   debug = 2;
   numBlocks = atoi(*argv); 
   pFileStem = *++argv; 

   assert(pFileStem);
   assert(numBlocks>0);
   assert(actBlock<numBlocks);
   
   initGMSPIPSIO();   
   
   if ( gdxSplit && (fType == GDXFILE) )
   {
      rc = gdxSplitting(numBlocks, actBlock, offset, gdxSplit==2, pFileStem, pGAMSSysDir);
      if (rc)
         printf("gdxSplitting failed (rc=%d)\n", rc);
      return rc;
   }
   
   if ( actBlock>=0 )
   {
      GMSPIPSBlockData_t block;
      
      rc = readOneBlock(numBlocks, actBlock, debug, offset, fType, pFileStem, pGAMSSysDir, &block);
      if (rc)
      {
         printf("readOneBlock with actBlock=%d failed (rc=%d)\n", actBlock, rc);
         return rc;
      }
      if ( printMat )
      {
         rc = writeBlock(NULL,&block,printMat);
         if (rc)
         {
            printf("writeBlock with actBlock=%d failed (rc=%d)\n", actBlock, rc);
            return rc;
         }
      }
      else
         printf("All done\n");

   }      
   else
   {
      int nBlock0=-1;
      GMSPIPSBlockData_t** blocks = (GMSPIPSBlockData_t**) calloc(numBlocks,sizeof(GMSPIPSBlockData_t*));
      for (int blk=0; blk<numBlocks; blk++)
      {
         blocks[blk] = (GMSPIPSBlockData_t*) malloc(sizeof(GMSPIPSBlockData_t));
         rc = readOneBlock(numBlocks, blk, debug, offset, fType, pFileStem, pGAMSSysDir, blocks[blk]);
         if (rc)
         {
            printf("readOneBlock with blk=%d failed (rc=%d)\n", blk, rc);
            return rc;
         }
         if ( 0==blk )
            nBlock0 = blocks[blk]->n0;
         else
            if (nBlock0 != blocks[blk]->n0)
            {
               printf("(nBlock0 != blocks[blk]->n0 (%d != %d, blk=%d)\n", nBlock0, blocks[blk]->n0, blk);
               return -1;
            }
      }
      
      if ( printMat )
      {
         for (int blk=0; blk<numBlocks; blk++)
         {
            rc = writeBlock(NULL,blocks[blk],printMat);
            if (rc)
            {
               printf("writeBlock with blk=%d failed (rc=%d)\n", blk, rc);
               return rc;
            }
         }
         printf("All done\n");
      }
      else 
         printf("All done\n");

      for (int blk=0; blk<numBlocks; blk++)
      {
         freeBlock(blocks[blk]);
         free(blocks[blk]);
      }
      free(blocks);
   }      
   return rc;
}
