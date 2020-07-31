/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysRootAug.h"
#include "DeSymIndefSolver.h"
#include "DeSymIndefSolver2.h"
#include "DeSymPSDSolver.h"
#include "PardisoSolver.h"
#include "PardisoIndefSolver.h"
#include "sData.h"
#include "sTree.h"
#include <limits>
#ifdef WITH_MUMPS_ROOT
#include "MumpsSolverRoot.h"
#endif

//#define DUMPKKT

#ifdef DUMPKKT
#include <iostream>
#include <fstream>
#endif

#include "pipsport.h"
#include <unistd.h>
#include "math.h"

#ifdef STOCH_TESTING
extern double g_iterNumber;
#endif
extern int gInnerBiCGIter;
extern int gInnerBiCGFails;



static void biCGStabCommunicateStatus(int flag, int it)
{
   gInnerBiCGIter = it;

   if( flag != 0 )
      gInnerBiCGFails++;
}

sLinsysRootAug::sLinsysRootAug(sFactory * factory_, sData * prob_)
  : sLinsysRoot(factory_, prob_), CtDC(nullptr)
{ 
  prob_->getLocalSizes(locnx, locmy, locmz, locmyl, locmzl);
  kkt = createKKT(prob_);
  solver = createSolver(prob_, kkt);
  assert(locmyl >= 0 && locmzl >= 0);
  redRhs = new SimpleVector(locnx+locmy+locmz+locmyl+locmzl);
}

sLinsysRootAug::sLinsysRootAug(sFactory* factory_,
			       sData* prob_,
			       OoqpVector* dd_, 
			       OoqpVector* dq_,
			       OoqpVector* nomegaInv_,
			       OoqpVector* rhs_)
  : sLinsysRoot(factory_, prob_, dd_, dq_, nomegaInv_, rhs_), CtDC(nullptr)
{ 
  assert(locmyl == 0 && locmzl == 0);

  kkt = createKKT(prob_);
  solver = createSolver(prob_, kkt);
  redRhs = new SimpleVector(locnx+locmy+locmz);
}

sLinsysRootAug::~sLinsysRootAug()
{
  if(CtDC) delete CtDC;
  delete redRhs;
}

SymMatrix* 
sLinsysRootAug::createKKT(sData* prob)
{
   const int n = locnx + locmy + locmyl + locmzl;

   if( hasSparseKkt )
   {
      int myRank; MPI_Comm_rank(mpiComm, &myRank);
      SparseSymMatrix* sparsekkt;

      if( myRank == 0)
         std::cout << "getSchurCompMaxNnz " << prob->getSchurCompMaxNnz() << std::endl;

      if( usePrecondDist )
      {
         sparsekkt = prob->createSchurCompSymbSparseUpperDist(childrenProperStart, childrenProperEnd);
      }
      else
      {
         sparsekkt = prob->createSchurCompSymbSparseUpper();
      }

      assert(sparsekkt->size() == n);

      return sparsekkt;
   }
   else
   {
      return new DenseSymMatrix(n);
   }
}


DoubleLinearSolver*
sLinsysRootAug::createSolver(sData* prob, SymMatrix* kktmat_)
{
   int myRank; MPI_Comm_rank(mpiComm, &myRank);


#ifdef WITH_MUMPS_ROOT
   if( hasSparseKkt )
   {
      if( 0 == myRank )
         cout << "Using MUMPS for summed Schur complement - sLinsysRootAug" << endl;

      SparseSymMatrix* kktmat = dynamic_cast<SparseSymMatrix*>(kktmat_);

      return new MumpsSolverRoot(mpiComm, kktmat);
   }
   else
#elif defined(WITH_PARDISO)
   if( hasSparseKkt )
   {
      if( 0 == myRank )
         cout << "Using Pardiso for summed Schur complement - sLinsysRootAug" << endl;

      SparseSymMatrix* kktmat = dynamic_cast<SparseSymMatrix*>(kktmat_);

      return new PardisoIndefSolver(kktmat);
   }
   else
#else
   assert(!hasSparseKkt);
#endif
   {
      if( 0 == myRank )
         cout << "Using LAPACK dsytrf for summed Schur complement - sLinsysRootAug" << endl;

      DenseSymMatrix* kktmat = dynamic_cast<DenseSymMatrix*>(kktmat_);

      return new DeSymIndefSolver(kktmat);
      //return new DeSymIndefSolver2(kktmat, locnx); // saddle point solver
      //return new DeSymPSDSolver(kktmat);
      //return new PardisoSolver(kktmat);
   }
}

#ifdef TIMING
static double t_start, troot_total, taux, tchild_total, tcomm_total;
#endif


void sLinsysRootAug::finalizeKKT(sData* prob, Variables* vars)
{
  stochNode->resMon.recFactTmLocal_start();
  stochNode->resMon.recSchurMultLocal_start();

  if( usePrecondDist )
  {
     // don't do anything, already done previously
  }
  else
  {
     if( hasSparseKkt )
        finalizeKKTsparse(prob, vars);
     else
        finalizeKKTdense(prob, vars);
  }

  stochNode->resMon.recSchurMultLocal_stop();
  stochNode->resMon.recFactTmLocal_stop();
}


void sLinsysRootAug::finalizeKKTdist(sData* prob)
{
   assert(kkt && hasSparseKkt && prob);

   SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);

   int myRank; MPI_Comm_rank(mpiComm, &myRank);
   int mpiCommSize; MPI_Comm_size(mpiComm, &mpiCommSize);
   const bool iAmLastRank = (myRank == mpiCommSize - 1);
   const int childStart = childrenProperStart;
   const int childEnd = childrenProperEnd;

#ifndef NDEBUG
   int* const jcolKkt = kkts.jcolM();
#endif
   int* const krowKkt = kkts.krowM();
   double* const MKkt = kkts.M();

   assert(childStart >= 0 && childStart < childEnd);
   assert(kkts.size() == locnx + locmy + locmyl + locmzl);
   assert(!kkts.isLower);
   assert(locmyl >= 0 && locmzl >= 0);
   assert(prob->getLocalQ().krowM()[locnx] == 0 && "Q currently not supported for dist. sparse kkt");

   /////////////////////////////////////////////////////////////
   // update the KKT with the diagonals
   // xDiag is in fact diag(Q)+X^{-1}S
   /////////////////////////////////////////////////////////////

   if( xDiag && iAmLastRank )
   {
      const SimpleVector& sxDiag = dynamic_cast<const SimpleVector&>(*xDiag);

      for( int i = 0; i < locnx; i++ )
      {
         const int diagIdx = krowKkt[i];
         assert(jcolKkt[diagIdx] == i);

         MKkt[diagIdx] += sxDiag[i];
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with   - C' * diag(zDiag) *C
   /////////////////////////////////////////////////////////////
   if( locmz > 0 && iAmLastRank )
   {
      assert(zDiag);

      SparseGenMatrix& C = prob->getLocalD();
      C.matTransDinvMultMat(*zDiag, &CtDC);
      assert(CtDC->size() == locnx);

      //aliases for internal buffers of CtDC
      SparseSymMatrix* CtDCsp = dynamic_cast<SparseSymMatrix*>(CtDC);
      const int* krowCtDC = CtDCsp->krowM();
      const int* jcolCtDC = CtDCsp->jcolM();
      const double* dCtDC = CtDCsp->M();

      for( int i = 0; i < locnx; i++ )
      {
         const int pend = krowCtDC[i + 1];
         for( int p = krowCtDC[i]; p < pend; p++ )
         {
            const int col = jcolCtDC[p];

            if( col >= i )
            {
               // get start position of dense kkt block
               const int blockStart = krowKkt[i];
               assert(col < locnx && jcolKkt[blockStart + col - i] == col);

               MKkt[blockStart + col - i] -= dCtDC[p];
            }
         }
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with At (symmetric update forced)
   /////////////////////////////////////////////////////////////
   if( locmy > 0 && iAmLastRank )
   {
      SparseGenMatrix& At = prob->getLocalB().getTranspose(); // yes, B
      const double* MAt = At.M();
      const int* krowAt = At.krowM();

      for( int i = 0; i < locnx; ++i )
      {
         const int pstart = krowAt[i];
         const int pend = krowAt[i + 1];

         // get start position of sparse kkt block
         const int blockStart = krowKkt[i] + locnx - i;

         assert(blockStart <= krowKkt[i + 1]);

         for( int p = pstart; p < pend; ++p )
         {
            assert(At.jcolM()[p] < locmy);
            assert(blockStart + (p - pstart) <= krowKkt[i + 1]);
            assert(jcolKkt[blockStart + (p - pstart)] == (locnx + At.jcolM()[p]));

            MKkt[blockStart + (p - pstart)] += MAt[p];
         }
      }
   }

   int local2linksStartEq;
   int local2linksEndEq;
   int local2linksStartIneq;
   int local2linksEndIneq;

   prob->getSCrangeMarkersMy(childStart, childEnd, local2linksStartEq, local2linksEndEq,
         local2linksStartIneq, local2linksEndIneq);

   const int n2linksRowsLocalEq = local2linksEndEq - local2linksStartEq;

   PIPSdebugMessage("rank %d FT local columns: %d-%d \n", myRank, local2linksStartEq, local2linksEndEq);
   PIPSdebugMessage("rank %d GT local columns: %d-%d \n", myRank, local2linksStartIneq, local2linksEndIneq);
   PIPSdebugMessage("rank %d FT local columns: %d-%d \n", myRank, local2linksStartEq, local2linksEndEq);

   /////////////////////////////////////////////////////////////
   // update the KKT with Ft
   /////////////////////////////////////////////////////////////
   if( locmyl > 0 )
   {
      SparseGenMatrix& Ft = prob->getLocalF().getTranspose();

      // add locally owned sparse part of Ft
      addLinkConsBlock0Matrix(prob, Ft, locnx + locmy, 0, local2linksStartEq, local2linksEndEq);

      if( myRank == 0 )
      {
         const int n2linksRowsEq = prob->n2linkRowsEq();
         const int bordersizeEq = locmyl - n2linksRowsEq;
         const int borderstartEq = locnx + locmy + n2linksRowsEq;

         PIPSdebugMessage("rank %d FT border columns: %d-%d \n", myRank, borderstartEq, borderstartEq + bordersizeEq);

         // add (shared) border part of Ft
         addLinkConsBlock0Matrix(prob, Ft, locnx + locmy, n2linksRowsLocalEq, borderstartEq, borderstartEq + bordersizeEq);
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with Gt and add z diagonal
   /////////////////////////////////////////////////////////////
   if( locmzl > 0 )
   {
      SparseGenMatrix& Gt = prob->getLocalG().getTranspose();
      const int n2linksRowsIneq = prob->n2linkRowsIneq();
      const int bordersizeIneq = locmzl - n2linksRowsIneq;
      const int borderstartIneq = locnx + locmy + locmyl + n2linksRowsIneq;

      // add locally owned sparse part of Gt
      addLinkConsBlock0Matrix(prob, Gt, locnx + locmy + locmyl, n2linksRowsLocalEq, local2linksStartIneq, local2linksEndIneq);

      if( myRank == 0 )
      {
         const int n2linksRowsLocalIneq = local2linksEndIneq - local2linksStartIneq;
         PIPSdebugMessage("rank %d GT border columns: %d-%d\n", myRank, borderstartIneq, borderstartIneq + bordersizeIneq);

         // add (shared) border part of Gt
         addLinkConsBlock0Matrix(prob, Gt, locnx + locmy + locmyl, n2linksRowsLocalEq + n2linksRowsLocalIneq,
               borderstartIneq, borderstartIneq + bordersizeIneq);
      }

      assert(zDiagLinkCons);
      const SimpleVector& szDiagLinkCons = dynamic_cast<const SimpleVector&>(*zDiagLinkCons);

      assert(local2linksStartIneq >= locnx + locmy + locmyl);
      assert(local2linksEndIneq <= locnx + locmy + locmyl + locmzl);

      const int szDiagLocalStart = local2linksStartIneq - (locnx + locmy + locmyl);
      assert(szDiagLocalStart >= 0);
      assert(szDiagLocalStart < locmzl || (szDiagLocalStart == locmzl && local2linksStartIneq == local2linksEndIneq));

      // add locally owned part of z diagonal
      for( int i = szDiagLocalStart, iKkt = local2linksStartIneq; iKkt < local2linksEndIneq; ++i, ++iKkt )
      {
         const int idx = krowKkt[iKkt];
         assert(jcolKkt[idx] == iKkt);
         assert(i < locmzl);

         MKkt[idx] += szDiagLinkCons[i];
      }

      if( myRank == 0 )
      {
         const int szDiagBorderStart = borderstartIneq - (locnx + locmy + locmyl);

         assert(szDiagBorderStart >= 0 && szDiagBorderStart <= locmzl);
         assert(szDiagBorderStart + bordersizeIneq == locmzl);

         // add border part of diagonal
         for( int i = szDiagBorderStart, iKkt = borderstartIneq; iKkt < borderstartIneq + bordersizeIneq; ++i, ++iKkt )
         {
            const int idx = krowKkt[iKkt];
            assert(jcolKkt[idx] == iKkt);
            assert(i < locmzl);

            MKkt[idx] += szDiagLinkCons[i];
         }
      }
   }
}


extern int gLackOfAccuracy;
void sLinsysRootAug::solveReduced( sData *prob, SimpleVector& b)
{
  int myRank; MPI_Comm_rank(mpiComm, &myRank);

#ifdef TIMING
  t_start=MPI_Wtime();
  troot_total=tchild_total=tcomm_total=0.0; 
#endif

  assert(locnx+locmy+locmz == b.length());
  SimpleVector& r = (*redRhs);
  assert(r.length() <= b.length());
  SparseGenMatrix& C = prob->getLocalD();

  ///////////////////////////////////////////////////////////////////////
  // LOCAL SOLVE
  ///////////////////////////////////////////////////////////////////////
 
  ///////////////////////////////////////////////////////////////////////
  // b=[b1;b2;b3] is a locnx+locmy+locmz vector 
  // the new rhs should be 
  //           r = [b1-C^T*(zDiag)^{-1}*b3; b2]
  ///////////////////////////////////////////////////////////////////////

  r.copyFromArray(b.elements()); //will copy only as many elems as r has

  // aliases to parts (no mem allocations)
  SimpleVector r3(&r[locnx+locmy], locmz); //r3 is used as a temp
                                           //buffer for b3
  SimpleVector r2(&r[locnx],       locmy);
  SimpleVector r1(&r[0],           locnx);

  ///////////////////////////////////////////////////////////////////////
  // compute r1 = b1 - C^T*(zDiag)^{-1}*b3
  ///////////////////////////////////////////////////////////////////////
  if(locmz>0) {
    assert(r3.length() == zDiag->length());
    r3.componentDiv(*zDiag);//r3 is a copy of b3
    C.transMult(1.0, r1, -1.0, r3);
  }
  ///////////////////////////////////////////////////////////////////////
  // r contains all the stuff -> solve for it
  ///////////////////////////////////////////////////////////////////////
  if( innerSCSolve == 0 ) {
    // Option 1. - solve with the factors
    solver->Dsolve(r);
  } else if( innerSCSolve == 1 ) {
    // Option 2 - solve with the factors and perform iter. ref.
    solveWithIterRef(prob, r);
  } else {
    assert( innerSCSolve == 2 );
    // Option 3 - use the factors as preconditioner and apply BiCGStab
    solveWithBiCGStab(prob, r);
  }

  ///////////////////////////////////////////////////////////////////////
  // r is the sln to the reduced system
  // the sln to the aug system should be 
  //      x = [r1; r2;  (zDiag)^{-1} * (b3-C*r1);
  ///////////////////////////////////////////////////////////////////////
  SimpleVector b1(&b[0],           locnx);
  SimpleVector b2(&b[locnx],       locmy);
  SimpleVector b3(&b[locnx+locmy], locmz);
  b1.copyFrom(r1);
  b2.copyFrom(r2);

  if(locmz>0) {
    C.mult(1.0, b3, -1.0, r1);
    b3.componentDiv(*zDiag);
  }
#ifdef TIMING
  if( myRank == 0 && innerSCSolve >= 1 )
    cout << "Root - Refin times: child=" << tchild_total << " root=" << troot_total
	 << " comm=" << tcomm_total << " total=" << MPI_Wtime()-t_start << endl;
#endif
}

void sLinsysRootAug::solveReducedLinkCons( sData *prob, SimpleVector& b)
{
  int myRank; MPI_Comm_rank(mpiComm, &myRank);

#ifdef TIMING
  t_start=MPI_Wtime();
  troot_total=tchild_total=tcomm_total=0.0;
#endif
  assert(locmyl >= 0 && locmzl >= 0);

  assert(locnx+locmy+locmz+locmyl+locmzl == b.length());
  SimpleVector& r = (*redRhs);

  assert(r.length() == b.length());

  SparseGenMatrix& C = prob->getLocalD();

  ///////////////////////////////////////////////////////////////////////
  // LOCAL SOLVE
  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // b=[b1;b2;b3;b4;b5] is a locnx+locmy+locmz+locmyl+locmz vector
  // the new rhs should be
  //           r = [b1-C^T*(zDiag)^{-1}*b3; b2; b4; b5]
  ///////////////////////////////////////////////////////////////////////

  //copy all elements from b into r except for the the residual values corresponding to z0
  assert(r.n > 0 && sizeof( double ) == sizeof(r[0]));

  memcpy( &r[0], &b[0], (locnx+locmy) * sizeof( double ) );
  if( locmyl > 0 )
     memcpy( &r[locnx+locmy], &b[locnx+locmy+locmz], locmyl * sizeof( double ) );
  if( locmzl > 0 )
     memcpy( &r[locnx+locmy+locmyl], &b[locnx+locmy+locmz+locmyl], locmzl * sizeof( double ) );

  // aliases to parts (no mem allocations)
  SimpleVector r1(&r[0],           locnx);

  ///////////////////////////////////////////////////////////////////////
  // compute r1 = b1 - C^T*(zDiag)^{-1}*b3
  ///////////////////////////////////////////////////////////////////////
  if(locmz>0) {
    memcpy( &r[locnx+locmy+locmyl+locmzl], &b[locnx+locmy], locmz * sizeof( double ) );
	 SimpleVector b3copy(&r[locnx+locmy+locmyl+locmzl], locmz); // b3copy is used as a temp buffer for b3
    assert(b3copy.length() == zDiag->length());
    b3copy.componentDiv(*zDiag);
    C.transMult(1.0, r1, -1.0, b3copy);
  }
  ///////////////////////////////////////////////////////////////////////
  // r contains all the stuff -> solve for it
  ///////////////////////////////////////////////////////////////////////

  // we do not need the last locmz elements of r
  SimpleVector rshort(&r[0], locnx+locmy+locmyl+locmzl);

  if( innerSCSolve == 0 ) {
    // Option 1. - solve with the factors
    solver->Dsolve(rshort);
  } else if( innerSCSolve == 1 ) {
    // Option 2 - solve with the factors and perform iter. ref.
    solveWithIterRef(prob, rshort);
  } else {
    assert( innerSCSolve == 2 );
    // Option 3 - use the factors as preconditioner and apply BiCGStab
    solveWithBiCGStab(prob, rshort);
  }

  ///////////////////////////////////////////////////////////////////////
  // r is the sln to the reduced system
  // the sln to the aug system should be
  //      x = [r1; r2;  (zDiag)^{-1} * (b3-C*r1);
  ///////////////////////////////////////////////////////////////////////
  SimpleVector b1(&b[0],           locnx);

  b1.copyFrom(r1);

  if( locmy > 0 )
  {
     SimpleVector r2(&r[locnx],       locmy);
     SimpleVector b2(&b[locnx],       locmy);
     b2.copyFrom(r2);
  }

  if( locmz > 0 )
  {
    SimpleVector b3(&b[locnx+locmy], locmz);
    C.mult(1.0, b3, -1.0, r1);
    b3.componentDiv(*zDiag);
  }

  if( locmyl > 0 )
  {
    SimpleVector rmyl(&r[locnx+locmy], locmyl);
    SimpleVector b4(&b[locnx+locmy+locmz], locmyl);
    b4.copyFrom(rmyl);
  }

  if( locmzl > 0 )
  {
    SimpleVector rmzl(&r[locnx+locmy+locmyl], locmzl);
    SimpleVector b5(&b[locnx+locmy+locmz+locmyl], locmzl);
    b5.copyFrom(rmzl);
  }

#ifdef TIMING
  if( myRank == 0 && innerSCSolve >= 1 )
    cout << "Root - Refin times: child=" << tchild_total << " root=" << troot_total
	 << " comm=" << tcomm_total << " total=" << MPI_Wtime()-t_start << endl;
#endif
}

/** Ht should be either Ft or Gt */
void sLinsysRootAug::addLinkConsBlock0Matrix( sData *prob, SparseGenMatrix& Ht, int nHtOffsetCols,
      int nKktOffsetCols, int startCol, int endCol)
{
   assert(startCol >= 0 && startCol <= endCol && nKktOffsetCols >= 0 && nKktOffsetCols <= startCol);

   if( startCol == endCol )
      return;

   SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);

   int* const jcolKkt = kkts.jcolM();
   int* const krowKkt = kkts.krowM();
   double* const MKkt = kkts.M();
   const double* const MHt = Ht.M();
   const int* const krowHt = Ht.krowM();
   const int* const jcolHt = Ht.jcolM();
   const int n0Links = prob->getN0LinkVars();

   /* main loop going over all rows of Ht */
   for( int i = 0; i < locnx; ++i )
   {
      const bool sparseRow = (i >= locnx - n0Links);

      // note: upper left block ignores 0-link variables pattern, since CtC pattern is not implemented
      int pKkt = krowKkt[i] + locnx - i;

      if( !sparseRow )
         pKkt += nKktOffsetCols;

      assert(pKkt <= krowKkt[i + 1]);
      assert(sparseRow || pKkt == krowKkt[i + 1] || jcolKkt[pKkt] <= startCol);

      if( jcolKkt[pKkt] >= endCol )
      {
#ifndef NDEBUG
         // make sure that there is no entry of Ht in the given range
         int pHt;

         for( pHt = krowHt[i]; pHt < krowHt[i + 1]; pHt++ )
         {
            const int colHt = jcolHt[pHt] + nHtOffsetCols;
            if( colHt >= startCol && colHt < endCol )
               break;
         }

         assert(pHt == krowHt[i + 1]);
#endif
         return;
      }

      bool hit = false;

      // get first in-range entry of Kkt
      for( ; pKkt < krowKkt[i + 1]; pKkt++ )
      {
         const int colKkt = jcolKkt[pKkt];
         if( colKkt >= startCol && colKkt < endCol )
         {
            hit = true;
            break;
         }

         if( colKkt >= endCol )
            break;
      }

      // no entry of Kkt in range?
      if( !hit )
      {
         assert(startCol == endCol || sparseRow);
         continue;
      }

      assert(pKkt < krowKkt[i + 1]);

      int pHt;
      int colHt = -1;
      hit = false;

      // get first in-range entry of Ht
      for( pHt = krowHt[i]; pHt < krowHt[i + 1]; pHt++ )
      {
         colHt = jcolHt[pHt] + nHtOffsetCols;
         if( colHt >= startCol && colHt < endCol )
         {
            hit = true;
            break;
         }

         if( colHt >= endCol )
            break;
      }

      // no entry of Ht in range?
      if( !hit )
         continue;

      assert(colHt >= startCol && colHt < endCol);

      // add in-range entries of Ht to Kkt
      for( ; pKkt < krowKkt[i + 1]; pKkt++ )
      {
         const int colKkt = jcolKkt[pKkt];

         if( colKkt >= endCol )
            break;

         if( colKkt == colHt )
         {
            assert(pHt < krowHt[i + 1]);

            MKkt[pKkt] += MHt[pHt++];

            // end of Ht row reached?
            if( pHt == krowHt[i + 1] )
               break;

            colHt = jcolHt[pHt] + nHtOffsetCols;
         }
      }

      assert(pHt == krowHt[i + 1] || jcolHt[pHt] + nHtOffsetCols >= endCol); // asserts that no entry of Ht has been missed
   }
}


/** rxy = beta*rxy + alpha * SC * x */
void sLinsysRootAug::SCmult( double beta, SimpleVector& rxy,
              double alpha, SimpleVector& x,
              sData* prob)
{
  //if (iAmDistrib) {
  //only one process subtracts [ (Q+Dx0+C'*Dz0*C)*xx + A'*xy + F'*xxl + G'*xyl ] from r
  //                           [  A*xx                                         ]
  //                           [  F*xx                                         ]
  //                           [  G*xx                           + Omega * xyl ]

  int myRank; MPI_Comm_rank(mpiComm, &myRank);
  int mpiCommSize; MPI_Comm_size(mpiComm, &mpiCommSize);
  const bool iAmLastRank = (myRank == mpiCommSize - 1);
  assert(mpiCommSize >= 1);

  if( iAmLastRank ) { // needs to be the last rank because only this rank is guaranteed to have CtDC
    assert(rxy.length() == locnx + locmy + locmyl + locmzl);

    //only this proc subtracts from rxy
    rxy.scalarMult(beta);
    SparseSymMatrix& Q = prob->getLocalQ();
    Q.mult(1.0,&rxy[0],1, alpha,&x[0],1);

    if(locmz>0) {
      SparseSymMatrix* CtDC_sp = dynamic_cast<SparseSymMatrix*>(CtDC);
      assert(CtDC_sp);

      CtDC_sp->mult(1.0,&rxy[0],1, -alpha,&x[0],1);
    }

    SimpleVector& xDiagv = dynamic_cast<SimpleVector&>(*xDiag);
    assert(xDiagv.length() == locnx);
    for(int i=0; i<xDiagv.length(); i++)
      rxy[i] += alpha*xDiagv[i]*x[i];

    SparseGenMatrix& A=prob->getLocalB();
    A.transMult(1.0,&rxy[0],1, alpha,&x[locnx],1);
    A.mult(1.0,&rxy[locnx],1, alpha,&x[0],1);

    assert(locmyl >= 0 && locmzl >= 0);

    if( locmyl > 0 ) {
       SparseGenMatrix& F = prob->getLocalF();
       F.transMult(1.0,&rxy[0],1, alpha,&x[locnx+locmy],1);
       F.mult(1.0,&rxy[locnx+locmy],1, alpha,&x[0],1);
    }

    if( locmzl > 0 ) {
       SparseGenMatrix& G = prob->getLocalG();
       G.transMult(1.0,&rxy[0],1, alpha,&x[locnx+locmy+locmyl],1);
       G.mult(1.0,&rxy[locnx+locmy+locmyl],1, alpha,&x[0],1);

       SimpleVector& zDiagLinkConsv = dynamic_cast<SimpleVector&>(*zDiagLinkCons);
       assert(zDiagLinkConsv.length() == locmzl);
       const int shift = locnx+locmy+locmyl;
       for(int i=0; i<zDiagLinkConsv.length(); i++)
         rxy[i+shift] += alpha*zDiagLinkConsv[i]*x[i+shift];
    }
  } else {
    //other processes set r to zero since they will get this portion from process 0
    rxy.setToZero();
  }

#ifdef TIMING
    taux=MPI_Wtime();
#endif
    // now children add [0 A^T C^T ]*inv(KKT)*[0;A;C] x
    SimpleVector xx((locmyl || locmzl) ? (locnx + locmy + locmyl + locmzl) : locnx);
    xx.copyFromArray(x.elements());
    xx.scalarMult(-alpha);

    for(size_t it=0; it<children.size(); it++) {
      children[it]->addTermToSchurResidual(prob->children[it],rxy,xx);
    }

#ifdef TIMING
    tchild_total +=  (MPI_Wtime()-taux);
#endif
    //~done computing residual

#ifdef TIMING
    taux=MPI_Wtime();
#endif
    //all-reduce residual
    if(iAmDistrib) {
      SimpleVector buf(rxy.length());
      buf.setToZero(); //we use dx as the recv buffer
      MPI_Allreduce(rxy.elements(), buf.elements(), locnx+locmy+locmyl+locmzl, MPI_DOUBLE, MPI_SUM, mpiComm);
      rxy.copyFrom(buf);
    }
#ifdef TIMING
    tcomm_total += (MPI_Wtime()-taux);
#endif

}



void sLinsysRootAug::solveWithIterRef( sData *prob, SimpleVector& r)
{
  SimpleVector r2(&r[locnx],       locmy);
  SimpleVector r1(&r[0],           locnx);

  //SimpleVector realRhs(&r[0], locnx+locmy);
#ifdef TIMING
  taux=MPI_Wtime();
#endif

  double rhsNorm=r.twonorm(); //r== the initial rhs of the reduced system here

  int myRank; MPI_Comm_rank(mpiComm, &myRank);
  SimpleVector rxy(locnx+locmy); rxy.copyFrom(r);
  SimpleVector   x(locnx+locmy); x.setToZero(); //solution
  SimpleVector  dx(locnx+locmy);                //update from iter refinement
  SimpleVector x_prev(locnx+locmy);
  int refinSteps=0;
  std::vector<double> histResid;
  int maxRefinSteps=(gLackOfAccuracy>0?9:8);
  do { //iterative refinement
#ifdef TIMING
    taux=MPI_Wtime();
#endif

    x_prev.copyFrom(x);
    //dx = Ainv * r 
    dx.copyFrom(rxy);
    solver->Dsolve(dx);
    //update x
    x.axpy(1.0,dx);

#ifdef TIMING
    troot_total += (MPI_Wtime()-taux);
#endif  

    if(gLackOfAccuracy<0) break;
    if(refinSteps==maxRefinSteps) break;

    //////////////////////////////////////////////////////////////////////
    //iterative refinement
    //////////////////////////////////////////////////////////////////////
    //compute residual
    
    //if (iAmDistrib) {
    //only one process substracts [ (Q+Dx0+C'*Dz0*C)*xx + A'*xy ] from r
    //                            [  A*xx                       ]
    if(myRank==0) {
      rxy.copyFrom(r);
      if(locmz>0) {
	SparseSymMatrix* CtDC_sp = dynamic_cast<SparseSymMatrix*>(CtDC);
	CtDC_sp->mult(1.0,&rxy[0],1, 1.0,&x[0],1);
      }
      SparseSymMatrix& Q = prob->getLocalQ();
      Q.mult(1.0,&rxy[0],1, -1.0,&x[0],1);
      
      SimpleVector& xDiagv = dynamic_cast<SimpleVector&>(*xDiag);
      assert(xDiagv.length() == locnx);
      for(int i=0; i<xDiagv.length(); i++)
	rxy[i] -= xDiagv[i]*x[i];
      
      SparseGenMatrix& A=prob->getLocalB();
      A.transMult(1.0,&rxy[0],1, -1.0,&x[locnx],1);
      A.mult(1.0,&rxy[locnx],1, -1.0,&x[0],1);
    } else {
      //other processes set r to zero since they will get this portion from process 0
      rxy.setToZero();
    }

#ifdef TIMING
    taux=MPI_Wtime();
#endif  
    // now children add [0 A^T C^T ]*inv(KKT)*[0;A;C] x
    SimpleVector xx(&x[0], locnx);
    for(size_t it=0; it<children.size(); it++) {
      children[it]->addTermToSchurResidual(prob->children[it],rxy,xx);  
    }
#ifdef TIMING
    tchild_total +=  (MPI_Wtime()-taux);
#endif
    //~done computing residual 

#ifdef TIMING
    taux=MPI_Wtime();
#endif
    //all-reduce residual
    if(iAmDistrib) {
      dx.setToZero(); //we use dx as the recv buffer
      MPI_Allreduce(rxy.elements(), dx.elements(), locnx+locmy, MPI_DOUBLE, MPI_SUM, mpiComm);
      rxy.copyFrom(dx);
    }
#ifdef TIMING
    tcomm_total += (MPI_Wtime()-taux);
#endif

    double relResNorm=rxy.twonorm()/rhsNorm;
    
    if(relResNorm<1.0e-10) {
      break;
    } else {
      double prevRelResNorm=1.0e10;
      if(histResid.size()) 
	prevRelResNorm=histResid[histResid.size()-1];

      //check for stop, divergence or slow convergence conditions
      if(relResNorm>prevRelResNorm) {
	// diverging; restore iteration
	if(myRank==0) {
	  cout << "1st stg - iter refinement diverges relResNorm=" << relResNorm 
	       << "  before was " << prevRelResNorm << endl;
	  cout << "Restoring iterate." << endl;
	}
	x.copyFrom(x_prev);
	break;
      }else {
	//check slow convergence for the last xxx iterates.
	// xxx is 1 for now
	//if(relResNorm>0.*prevRelResNorm) {

	//  if(myRank==0) {
	//    cout << "1st stg - iter refinement stuck relResNorm=" << relResNorm 
	//	 << "  before was " << prevRelResNorm << endl;
	//    cout << "exiting refinement." << endl;
	//  }
	//  break;
	//
	//} else {
	//  //really nothing, continue
	//}
      }
      histResid.push_back(relResNorm);
      if(myRank==0)
	cout << "1st stg - sol does NOT  have enough accuracy (" << relResNorm << ") after " 
	     << refinSteps << " refinement steps" << endl;
    }
    refinSteps++;
  }while(refinSteps<=maxRefinSteps);

#ifdef TIMING
  taux = MPI_Wtime();
#endif

  r1.copyFrom(x);
  r2.copyFromArray(&x[locnx]);

#ifdef TIMING
  troot_total += (MPI_Wtime()-taux);
#endif  
}


void sLinsysRootAug::solveWithBiCGStab( sData *prob, SimpleVector& b)
{
  int n = b.length();

  const int maxit=75; //500
  const double tol=1e-10, EPS=1e-15; // EPS=2e-16

  int myRank; MPI_Comm_rank(mpiComm, &myRank);

  SimpleVector r(n);           //residual
  SimpleVector s(n);           //residual associated with half iterate
  SimpleVector rt(n);          //shadow residual
  SimpleVector xmin(n);        //minimal residual iterate
  SimpleVector x(n);           //iterate
  SimpleVector xhalf(n);       // half iterate of BiCG
  SimpleVector p(n),paux(n);
  SimpleVector v(n), t(n);
  int flag;
  double n2b;                  //norm of b 
  double normr, normrmin;      //norm of the residual and norm of residual at min-resid iterate
  double normr_act;
  double tolb;                 //relative tolerance
  double rho, omega, alpha;
  int stag, maxmsteps, maxstagsteps, moresteps;
  //double imin;
  //maxit = n/2+1;

  //////////////////////////////////////////////////////////////////
  //  Problem Setup and initialization
  //////////////////////////////////////////////////////////////////

  n2b = b.twonorm();
  tolb = n2b*tol;

  tolb = max(tolb, EPS);

#ifdef TIMING
  double relres;
  double iter=0.0;
  if( myRank == 0 )
     std::cout << "initial norm of b " << n2b << std::endl;
  taux = MPI_Wtime();
#endif
  //initial guess
  x.copyFrom(b);

  solver->Dsolve(x);
  //initial residual
  r.copyFrom(b);

#ifdef TIMING
  troot_total += (MPI_Wtime()-taux);
  taux = MPI_Wtime();
#endif 

  //applyA(1.0, r, -1.0, x);
  SCmult(1.0,r, -1.0,x, prob);

#ifdef TIMING
    tchild_total +=  (MPI_Wtime()-taux);
#endif

  normr = r.twonorm(); normr_act = normr;

  if( normr<=tolb ) {
    //initial guess is good enough
    b.copyFrom(x); flag=0; return;
  }

  if( myRank == 0 )
      std::cout << "innerBICG starts: " << normr << " > " << tolb << std::endl;

  rt.copyFrom(r); //Shadow residual
  double* resvec = new double[2*maxit+1];
  resvec[0] = normr; normrmin=normr;
  rho=1.0; omega=1.0;
  stag=0; maxmsteps=min(min(n/50, 5), n-maxit); 
  maxstagsteps=3; moresteps=0;

  //////////////////////////////////////////////////////////////////
  // loop over maxit iterations
  //////////////////////////////////////////////////////////////////
  int ii=0; while(ii<maxit) {
    //cout << ii << " ";
    flag=-1;
    ///////////////////////////////
    // First half of the iterate
    ///////////////////////////////
    double rho1=rho; double beta;
    rho = rt.dotProductWith(r); 
    //printf("rho=%g\n", rho);
    if(0.0==rho) { flag=4;  break; }

    if(ii==0) p.copyFrom(r);
    else {
      beta = (rho/rho1)*(alpha/omega);
      if(beta==0.0) { flag=4;  break; }

      //-------- p = r + beta*(p - omega*v) --------
      p.axpy(-omega, v); p.scale(beta); p.axpy(1.0, r);
    }

#ifdef TIMING
    taux = MPI_Wtime();
#endif
    //------ v = A*(M2inv*(M1inv*p)) and ph=M2inv*(M1inv*p)
    //first use v as temp storage
    //applyM1(0.0, v,    1.0, p);
    //applyM2(0.0, paux, 1.0, v);
    //applyA (0.0, v,    1.0, paux); 
    paux.copyFrom(p);
    solver->solve(paux);

#ifdef TIMING
  troot_total += (MPI_Wtime()-taux);
#endif 
    
    SCmult(0.0,v, 1.0,paux, prob);
    
    SimpleVector& ph = paux;

    double rtv = rt.dotProductWith(v);
    if(rtv==0.0) { flag=4; break; }

    alpha = rho/rtv;
    if(fabs(alpha)*ph.twonorm()<EPS*x.twonorm()) stag++;
    else                                         stag=0;

    // xhalf = x + alpha*ph and the associated residual
    xhalf.copyFrom(x); xhalf.axpy( alpha, ph);
    s.    copyFrom(r);     s.axpy(-alpha, v);
    normr = s.twonorm(); normr_act = normr;
    resvec[2*ii] = normr;

    //printf("iter %g normr=%g\n", ii+0.5, normr);
    //-------- check for convergence in the middle of the iterate.  -------- 
    if(normr<=tolb || stag>=maxstagsteps || moresteps) {
      s.copyFrom(b);
      //applyA(1.0, s, -1.0, xhalf); // s=b-Ax
      SCmult(1.0,s, -1.0,xhalf, prob);
      normr_act = s.twonorm();
      
      if(normr<=tolb) {
	//converged
	x.copyFrom(xhalf);	
	flag = 0;
#ifdef TIMING
	iter = 0.5+ii;
#endif
	break;
      } else {
	if(stag>=maxstagsteps && moresteps==0) {
	  stag=0;
	}
	moresteps++;
	if(moresteps>=maxmsteps) {
	  //method stagnated
	  flag=3; x.copyFrom(xhalf);
	  break;
	}
      }
    }
    if(stag>=maxstagsteps) { flag=3; break;} //stagnation

    //update quantities related to minimal norm iterate
    if(normr_act<normrmin) {
      xmin.copyFrom(xhalf); normrmin=normr_act;
      //imin=0.5+ii;
    }

#ifdef TIMING
    taux = MPI_Wtime();
#endif
    ///////////////////////////////
    // Second half of the iterate
    //////////////////////////////
    //applyM1(0.0, t,    1.0, s); //applyM1(s,     stemp);
    //applyM2(0.0, paux, 1.0, t); //applyM2(stemp, sh);
    //applyA (0.0, t,    1.0, paux); //applyA (sh, t);
    //kkt->mult(0.0,paux, 1.0,s);
    paux.copyFrom(s);
    solver->solve(paux);
#ifdef TIMING
    troot_total += (MPI_Wtime()-taux);
#endif

    SCmult(0.0,t, 1.0,paux, prob);

    SimpleVector& sh = paux; 
    double tt = t.dotProductWith(t);
    if(tt==0.0) { flag=4; break;}

    omega=t.dotProductWith(s); omega /= tt;

    if(fabs(omega)*sh.twonorm() < EPS*xhalf.twonorm()) stag++;
    else                                               stag=0;

    x.copyFrom(xhalf); x.axpy( omega, sh); // x=xhalf+omega*sh
    r.copyFrom(s);     r.axpy(-omega, t ); // r=s-omega*t

    normr = r.twonorm(); normr_act = normr;
    resvec[2*ii+1] = normr;

    //printf("stag=%d  maxstagsteps=%d moresteps=%d  normr=%g\n",
    //	   stag, maxstagsteps, moresteps, normr);    

    //-------- check for convergence at the end of the iterate.  -------- 
    if(normr<=tolb || stag>=maxstagsteps || moresteps) {
      r.copyFrom(b); 
      //applyA(1.0, r, -1.0, x); //r=b-Ax
      SCmult(1.0,r, -1.0,x, prob);
      normr_act=r.twonorm();

      if(normr<=tolb) {
         flag = 0;
#ifdef TIMING
         iter = 1.0+ii;
#endif
         break;
      }
      else {
	if(stag>=maxstagsteps && moresteps==0) {
	  stag = 0;
	}
	moresteps++;
	if(moresteps>=maxmsteps) {
	  //method stagnated
	  flag=3; break;
	}
      }
    } // end convergence check
    if(stag>=maxstagsteps) { flag=3; break;} //stagnation

    //update quantities related to minimal norm iterate
    if(normr_act<normrmin) {
      xmin.copyFrom(x); normrmin=normr_act;
      //imin=1.5+ii;
    }
    //printf("iter %g normr=%g\n", ii+1.0, normr);
    ///////////////////////////////
    // Next iterate
    ///////////////////////////////
    ii++;
    
  }//end while

  if(ii>=maxit) {
#ifdef TIMING
    iter=ii;
#endif
    flag=10;
  }
  
  if(flag==0 || flag==-1) {
#ifdef TIMING
    relres = normr_act/n2b;
    if(myRank==0) {
      printf("INNER BiCGStab converged: normResid=%g relResid=%g iter=%g\n",
	        normr_act, relres, iter);
    }
#endif
  } 
  else 
  {
    if(ii==maxit) flag=10;//aaa
    //FAILURE -> return minimum resid-norm iterate
    r.copyFrom(b); 
    //applyA(1.0, r, -1.0, xmin);
    SCmult(1.0,r, -1.0,xmin, prob);

    normr=r.twonorm();
    if(normr >= normr_act) {
      x.copyFrom(xmin);
      //iter=imin;
#ifdef TIMING
      relres=normr/n2b;
#endif
    } else {
#ifdef TIMING
      iter=1.0+ii;
      relres = normr/n2b;
#endif
    }

#ifdef TIMING
    if(myRank==0) {
      printf("INNERBiCGStab did not NOT converged after %g[%d] iterations.\n", iter,ii);
      printf("\t - Error code %d\n\t - Act res=%g\n\t - Rel res=%g %g\n\n", 
	     flag, normr, relres, normrmin);
    }
#endif
  }

  if( myRank == 0 )
     std::cout << "innerBICG: " << "ii=" << ii << " flag=" << flag << " normr=" << normr << " normr_act="
        << normr_act << " tolb=" << tolb << std::endl;

  biCGStabCommunicateStatus(flag, ii);

  b.copyFrom(x);
  delete[] resvec;
}

void sLinsysRootAug::finalizeKKTsparse(sData* prob, Variables* vars)
{
   SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);

#ifndef NDEBUG
   int* const jcolKkt = kkts.jcolM();
#endif
   int* const krowKkt = kkts.krowM();
   double* const MKkt = kkts.M();
   const int n0Links = prob->getN0LinkVars();

   assert(kkts.size() == locnx + locmy + locmyl + locmzl);
   assert(!kkts.isLower);
   assert(locmyl >= 0 && locmzl >= 0);

#if 0
   int myRank; MPI_Comm_rank(mpiComm, &myRank);

   if( myRank == 0)
   {
      xDiag->writefToStreamStats(std::cout, "xDiag");
      zDiag->writefToStreamStats(std::cout, "zDiag");

      const SimpleVector& szDiagLinkCons = dynamic_cast<const SimpleVector&>(*zDiagLinkCons);
      assert(szDiagLinkCons.length() == locmzl);
      int zerocount = 0;
      for( int i = 0; i < locmzl; i++ )
      {
         if( szDiagLinkCons[i] == 0)
            zerocount++;
      }

      zDiagLinkCons->writefToStreamStats(std::cout, "zDiagLinkCons");

      std::cout << "zDiagLinkCons zeroes: " << zerocount << std::endl;
   }
#endif

   //////////////////////////////////////////////////////
   // compute Q+diag(xdiag) - C' * diag(zDiag) * C
   // and update the KKT
   //////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////////
   // update the KKT with Q (DO NOT PUT DIAG)
   /////////////////////////////////////////////////////////////
   assert(prob->getLocalQ().krowM()[locnx] == 0 && "Q currently not supported for sparse kkt");

   /////////////////////////////////////////////////////////////
   // update the KKT with the diagonals
   // xDiag is in fact diag(Q)+X^{-1}S
   /////////////////////////////////////////////////////////////
   if( xDiag )
   {
      const SimpleVector& sxDiag = dynamic_cast<const SimpleVector&>(*xDiag);

      for( int i = 0; i < locnx; i++ )
      {
         const int diagIdx = krowKkt[i];
         assert(jcolKkt[diagIdx] == i);

         MKkt[diagIdx] += sxDiag[i];
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with   - C' * diag(zDiag) *C
   /////////////////////////////////////////////////////////////
   if( locmz > 0 )
   {
      assert(zDiag);

      SparseGenMatrix& C = prob->getLocalD();
      C.matTransDinvMultMat(*zDiag, &CtDC);
      assert(CtDC->size() == locnx);

      //aliases for internal buffers of CtDC
      SparseSymMatrix* CtDCsp = dynamic_cast<SparseSymMatrix*>(CtDC);
      const int* krowCtDC = CtDCsp->krowM();
      const int* jcolCtDC = CtDCsp->jcolM();
      const double* dCtDC = CtDCsp->M();

      for( int i = 0; i < locnx; i++ )
      {
         const int pend = krowCtDC[i + 1];
         for( int p = krowCtDC[i]; p < pend; p++ )
         {
            const int col = jcolCtDC[p];

            if( col >= i )
            {
               // get start position of dense kkt block
               const int blockStart = krowKkt[i];
               assert(col < locnx && jcolKkt[blockStart + col - i] == col);

               MKkt[blockStart + col - i] -= dCtDC[p];
            }
         }
      }
   } //~end if locmz>0

   /////////////////////////////////////////////////////////////
   // update the KKT with At (symmetric update forced)
   /////////////////////////////////////////////////////////////
   if( locmy > 0 )
   {
      SparseGenMatrix& At = prob->getLocalB().getTranspose(); // yes, B
      const double* MAt = At.M();
      const int* krowAt = At.krowM();

      for( int i = 0; i < locnx; ++i )
      {
         const int pstart = krowAt[i];
         const int pend = krowAt[i + 1];

         // get start position of sparse kkt block
         const int blockStart = krowKkt[i] + locnx - i;

         assert(blockStart <= krowKkt[i + 1]);

         for( int p = pstart; p < pend; ++p )
         {
            assert(At.jcolM()[p] < locmy);
            assert(blockStart + (p - pstart) <= krowKkt[i + 1]);
            assert(jcolKkt[blockStart + (p - pstart)] == (locnx + At.jcolM()[p]));

            MKkt[blockStart + (p - pstart)] += MAt[p];
         }
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with Ft
   /////////////////////////////////////////////////////////////
   if( locmyl > 0 )
   {
      SparseGenMatrix& Ft = prob->getLocalF().getTranspose();
      const double* MFt = Ft.M();
      const int* krowFt = Ft.krowM();
      const int* jcolFt = Ft.jcolM();

      int* krowGt = nullptr;

      if( locmzl > 0 )
      {
         SparseGenMatrix& Gt = prob->getLocalG().getTranspose();
         krowGt = Gt.krowM();
      }

      for( int i = 0; i < locnx; ++i )
      {
         const bool sparseRow = (i >= locnx - n0Links);
         const int pend = krowFt[i + 1];

         if( sparseRow )
         {
            int blockStart =  krowKkt[i + 1] - (krowFt[i + 1] - krowFt[i]);

            if( locmzl > 0 )
               blockStart -= (krowGt[i + 1] - krowGt[i]);

            assert(blockStart >= krowKkt[i]);

            for( int p = krowFt[i], shift = 0; p < pend; ++p, ++shift )
            {
               assert(jcolFt[p] < locmyl && jcolKkt[blockStart + shift] == (locnx + locmy + jcolFt[p]));

               MKkt[blockStart + shift] += MFt[p];
            }
         }
         else
         {
            const int blockStart = krowKkt[i + 1] - locmyl - locmzl;
            assert(blockStart >= krowKkt[i]);

            for( int p = krowFt[i]; p < pend; ++p )
            {
               const int col = jcolFt[p];
               assert(col < locmyl && jcolKkt[blockStart + col] == (locnx + locmy + col));

               MKkt[blockStart + col] += MFt[p];
            }
         }
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with Gt and add z diagonal
   /////////////////////////////////////////////////////////////
   if( locmzl > 0 )
   {
      SparseGenMatrix& Gt = prob->getLocalG().getTranspose();
      const double* MGt = Gt.M();
      const int* krowGt = Gt.krowM();
      const int* jcolGt = Gt.jcolM();

      for( int i = 0; i < locnx; ++i )
      {
         const int pend = krowGt[i + 1];
         const bool sparseRow = (i >= locnx - n0Links);

         if( sparseRow )
         {
            const int blockStart = krowKkt[i + 1] - (krowGt[i + 1] - krowGt[i]);

            assert(blockStart >= krowKkt[i]);

            for( int p = krowGt[i], shift = 0; p < pend; ++p, ++shift )
            {
               assert(jcolGt[p] < locmzl && jcolKkt[blockStart + shift] == (locnx + locmy + locmyl + jcolGt[p]));

               MKkt[blockStart + shift] += MGt[p];
            }
         }
         else
         {
            const int blockStart = krowKkt[i + 1] - locmzl;

            assert(blockStart >= krowKkt[i]);

            for( int p = krowGt[i]; p < pend; ++p )
            {
               const int col = jcolGt[p];
               assert(col < locmzl && jcolKkt[blockStart + col] == (locnx + locmy + locmyl + col));

               MKkt[blockStart + col] += MGt[p];
            }
         }
      }

      assert(zDiagLinkCons);

      const SimpleVector& szDiagLinkCons = dynamic_cast<const SimpleVector&>(*zDiagLinkCons);

      for( int i = 0, iKkt = locnx + locmy + locmyl; i < locmzl; ++i, ++iKkt )
      {
         const int idx = krowKkt[iKkt];
         assert(jcolKkt[idx] == iKkt);

         MKkt[idx] += szDiagLinkCons[i];
      }
   }

#ifdef DUMPKKT
   ofstream myfile;
   myfile.open("../sparsekkt");

   int zerocount = 0;
   const int sizeKkt = locnx + locmy + locmyl + locmzl;

   for( int r = 0; r < sizeKkt; r++ )
   {
      for( int i = krowKkt[r]; i < krowKkt[r + 1]; i++ )
      {
         const double val = MKkt[i];
         const double col = jcolKkt[i];
         if( val != 0.0 )
            myfile << r << " " << col << " " << val << std::endl;
         else
            zerocount++;
      }
   }

   std::cout << "zero-count " << zerocount << " of " << krowKkt[sizeKkt] << std::endl;

   myfile.close();

   assert(0);

#endif

}


void sLinsysRootAug::finalizeKKTdense(sData* prob, Variables* vars)
{
   int j, p, pend;

   DenseSymMatrix * const kktd = dynamic_cast<DenseSymMatrix*>(kkt);

   //alias for internal buffer of kkt
   double** const dKkt = kktd->Mat();

   //////////////////////////////////////////////////////
   // compute Q+diag(xdiag) - C' * diag(zDiag) * C
   // and update the KKT
   //////////////////////////////////////////////////////


   /////////////////////////////////////////////////////////////
   // update the KKT with Q (DO NOT PUT DIAG)
   /////////////////////////////////////////////////////////////
   SparseSymMatrix& Q = prob->getLocalQ();
   int* krowQ=Q.krowM(); int* jcolQ=Q.jcolM(); double* dQ=Q.M();
   for(int i=0; i<locnx; i++) {
     pend = krowQ[i+1];
     for(p=krowQ[i]; p<pend; p++) {
       j = jcolQ[p];
       if(i==j) continue;
       double val = dQ[p];
       dKkt[i][j] += val;
       dKkt[j][i] += val;

       assert(0 && "non-empty Q currently not supported");
     }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with the diagonals
   // xDiag is in fact diag(Q)+X^{-1}S
   /////////////////////////////////////////////////////////////
   //kktd->atPutDiagonal( 0, *xDiag );
   if( xDiag )
   {
      SimpleVector& sxDiag = dynamic_cast<SimpleVector&>(*xDiag);
      for(int i=0; i<locnx; i++) dKkt[i][i] += sxDiag[i];
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with   - C' * diag(zDiag) *C
   /////////////////////////////////////////////////////////////
   if(locmz>0) {
     assert(zDiag);

     SparseGenMatrix& C = prob->getLocalD();
     C.matTransDinvMultMat(*zDiag, &CtDC);

     assert(CtDC->size() == locnx);

     //aliases for internal buffers of CtDC
     SparseSymMatrix* CtDCsp = reinterpret_cast<SparseSymMatrix*>(CtDC);
     int* krowCtDC=CtDCsp->krowM(); int* jcolCtDC=CtDCsp->jcolM(); double* dCtDC=CtDCsp->M();

     for( int i = 0; i < locnx; i++ )
     {
       pend = krowCtDC[i + 1];
       for( p = krowCtDC[i]; p < pend; p++ )
       {
          j = jcolCtDC[p];

          if( j <= i )
             dKkt[i][j] -= dCtDC[p];
        }
     }
   } //~end if locmz>0
   /////////////////////////////////////////////////////////////
   // update the KKT with A (symmetric update forced)
   /////////////////////////////////////////////////////////////
   if(locmy>0)
   {
     SparseGenMatrix& A = prob->getLocalB(); // yes, B
     const double* dA = A.M();
     const int* krowA = A.krowM();
     const int* jcolA = A.jcolM();

     int iKkt = locnx;
     for( int i = 0; i < locmy; ++i, ++iKkt ) {

       for( p = krowA[i], pend = krowA[i + 1]; p < pend; ++p ) {
         j = jcolA[p];
         assert(j < locnx);

         dKkt[iKkt][j] += dA[p];
       }
     }
   }
   //prob->getLocalB().getStorageRef().dump("stage1eqmat2.dump");


   /////////////////////////////////////////////////////////////
   // update the KKT with F
   /////////////////////////////////////////////////////////////
   if( locmyl > 0 )
   {
     SparseGenMatrix& F = prob->getLocalF();
     const double* dF = F.M();
     const int* krowF = F.krowM();
     const int* jcolF = F.jcolM();

     int iKkt = locnx + locmy;
     for( int i = 0; i < locmyl; ++i, ++iKkt ) {
       for( p = krowF[i], pend = krowF[i+1]; p < pend; ++p ) {
         j = jcolF[p];
         assert(j < locnx);

         const double val = dF[p];
         dKkt[iKkt][j] += val;
       }
     }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with G and put z diagonal
   /////////////////////////////////////////////////////////////
   if( locmzl > 0 )
   {
     SparseGenMatrix& G = prob->getLocalG();
     assert(zDiagLinkCons);
     SimpleVector& szDiagLinkCons = dynamic_cast<SimpleVector&>(*zDiagLinkCons);

     const double* dG = G.M();
     const int* krowG = G.krowM();
     const int* jcolG = G.jcolM();

     int iKkt = locnx + locmy + locmyl;
     for( int i = 0; i < locmzl; ++i, ++iKkt ) {

       dKkt[iKkt][iKkt] += szDiagLinkCons[i];
       for( p = krowG[i], pend = krowG[i+1]; p < pend; ++p ) {
         j = jcolG[p];
         assert(j < locnx);

         const double val = dG[p];
         dKkt[iKkt][j] += val;
       }
     }
   }

#ifdef DUMPKKT
   const int msize = locnx + locmy + locmyl + locmzl;

   ofstream myfile;
   myfile.open("../densekkt");

   for( int col = 0; col < msize; col++ )
      for( int row = col; row < msize; row++ )
         if( dKkt[row][col] != 0.0 )
            myfile << col << " " << row << " " << dKkt[row][col] << std::endl;

   myfile.close();

   assert(0);
#endif

   /////////////////////////////////////////////////////////////
   // update the KKT zeros for the lower right block
   /////////////////////////////////////////////////////////////
   //kktd->storage().atPutZeros(locnx, locnx, locmy+locmz, locmy+locmz);
   //myAtPutZeros(kktd, locnx, locnx, locmy, locmy);
}

