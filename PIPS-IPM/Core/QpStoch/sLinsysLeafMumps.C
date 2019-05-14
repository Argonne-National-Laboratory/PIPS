/*
 * sLinsysLeafMumps.C
 *
 *      Author: bzfrehfe
 */


#include "sLinsysLeafMumps.h"
#include "MumpsSolver.h"


sLinsysLeafMumps::~sLinsysLeafMumps()
{
   delete schurRightMatrix_csc;
   delete schurRightNzColId;
}

void sLinsysLeafMumps::addTermToSparseSchurCompl( sData *prob,
                      SparseSymMatrix& SC)
{
  this->addTermToSchurComplMumps(prob, true, SC);
}

void sLinsysLeafMumps::addTermToDenseSchurCompl( sData *prob,
                      DenseSymMatrix& SC)
{
  this->addTermToSchurComplMumps(prob, false, SC);
}


void sLinsysLeafMumps::buildSchurRightMatrix(sData *prob, SymMatrix& SC)
{
   assert(prob);
   assert(!schurRightMatrix_csc);
   assert(!schurRightNzColId);
   assert(nSC == -1);
   assert(mSchurRight == -1);


   // create new (Fortran indexed, implicit CSC) matrix
   //  (R^T A^T C^T)
   //  (0          )
   //  (F          )
   //  (G          )
   //
   // ... and remove empty rows

   assert(prob);

   SparseGenMatrix& A = prob->getLocalA();
   SparseGenMatrix& C = prob->getLocalC();
   SparseGenMatrix& F = prob->getLocalF();
   SparseGenMatrix& G = prob->getLocalG();
   SparseGenMatrix& R = prob->getLocalCrossHessian();

   int N, nxP;

   R.getSize(N, nxP);
   const bool withR = (nxP != -1);

   A.getSize(N, nxP);
   const bool withA = (nxP != -1);

   assert(N == locmy);
   assert(locmyl >= 0);
   assert(locmzl >= 0);

   nSC = SC.size();
   assert(nSC >= nxP);

   const int nxMyP = nSC - locmyl - locmzl;
   const int nxMyMzP = nSC - locmzl;

   if( nxP == -1 )
      C.getSize(N, nxP);

   int N2, nxP2;
   C.getSize(N2, nxP2);
   const bool withC = (nxP2 != -1);

   if( nxP == -1 )
      nxP = nSC;

   mSchurRight = locnx + locmy + locmz;

   SimpleVector nnzPerColRAC(nxP);

   if( withR )
      R.addNnzPerCol(nnzPerColRAC);

   if( withA )
      A.addNnzPerCol(nnzPerColRAC);

   if( withC )
      C.addNnzPerCol(nnzPerColRAC);

   const int withF = (locmyl > 0);
   const int withG = (locmzl > 0);
}


void sLinsysLeafMumps::addTermToSchurComplMumps(sData *prob, bool sparseSC,
          SymMatrix& SC)
{
   const int blocksizemax = 1000; // todo nThreads, or 64?

   if( !schurRightMatrix_csc )
      buildSchurRightMatrix(prob, SC);

   assert(nThreads >= 1);
   //std::cout << "blocksizemax " << blocksizemax << std::endl;

   // save columns in this array todo member variable, initialized at first call (if == NULL)
   double* colsBlockDense = new double[blocksizemax * mSchurRight];

   MumpsSolver* const solverMumps = dynamic_cast<MumpsSolver*>(solver);
   assert(solverMumps);


   //  do block-wise computation of
   //
   //     SC +=  B^T K^-1 B
   //
   //  ...and skip zero columns of B


   solverMumps->solve(*schurRightMatrix_csc, colsBlockDense);

   multLeftSchurComplBlocked(prob, colsBlockDense, schurRightNzColId, blocksizemax, sparseSC, SC);


#if 0
   // debug stuff
   int myrank;
   static int iteration = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   ofstream myfile;
   char filename[50];
   sprintf(filename, "../blocked_%d_%d.txt", myrank, iteration);
   myfile.open(filename);
   iteration++;
   SC.writeToStream(myfile); // todo write out in each iteration with global counter and MPI rank!
   myfile.close();

   assert(0);
#endif

   delete[] colsBlockDense;
}
