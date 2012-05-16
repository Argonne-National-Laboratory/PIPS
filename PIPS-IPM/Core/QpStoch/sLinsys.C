/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#include "sLinsys.h"
#include "StochTree.h"
#include "sFactory.h"
#include "sData.h"
#include "SparseLinearAlgebraPackage.h"

#ifndef MIN
#define MIN(a,b) ((a > b) ? b : a)
#endif

sLinsys::sLinsys(sFactory* factory_, sData* prob)
  : QpGenLinsys(), kkt(NULL), solver(NULL)
{
  factory = factory_;

  nx = prob->nx; my = prob->my; mz = prob->mz;
  ixlow = prob->ixlow;
  ixupp = prob->ixupp;
  iclow = prob->iclow;
  icupp = prob->icupp;

  nxlow = prob->nxlow;
  nxupp = prob->nxupp;
  mclow = prob->mclow;
  mcupp = prob->mcupp;

  if( nxupp + nxlow > 0 ) {
    dd      = factory_->tree->newPrimalVector();
    dq      = factory_->tree->newPrimalVector();
    prob->getDiagonalOfQ( *dq );
  }
  nomegaInv   = factory_->tree->newDualZVector();
  rhs         = factory_->tree->newRhs();

  useRefs=0;
  data = prob;
  stochNode = prob->stochNode;
}

sLinsys::sLinsys(sFactory* factory_,
				   sData* prob,				    
				   OoqpVector* dd_, 
				   OoqpVector* dq_,
				   OoqpVector* nomegaInv_,
				   OoqpVector* rhs_)
  : QpGenLinsys(), kkt(NULL), solver(NULL)
{
  factory = factory_;


  nx = prob->nx; my = prob->my; mz = prob->mz;
  ixlow = prob->ixlow;
  ixupp = prob->ixupp;
  iclow = prob->iclow;
  icupp = prob->icupp;

  nxlow = prob->nxlow;
  nxupp = prob->nxupp;
  mclow = prob->mclow;
  mcupp = prob->mcupp;

  if( nxupp + nxlow > 0 ) {
    //dd      = OoqpVectorHandle(dd_);
    dd= dd_;
    //dq      = OoqpVectorHandle(dq_);
    dq = dq_;
  }
  //nomegaInv   = OoqpVectorHandle(nomegaInv_);
  //rhs         = OoqpVectorHandle(rhs_);
  nomegaInv = nomegaInv_;
  rhs = rhs_;

  useRefs=1;
  data = prob;
  stochNode = prob->stochNode;
}


sLinsys::~sLinsys()
{
  if (solver) delete solver;
  if (kkt)    delete kkt;
}

void sLinsys::joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
				OoqpVector& rhs2_in, OoqpVector& rhs3_in )
{
  StochVector& rhs  = dynamic_cast<StochVector&>(rhs_in);
  StochVector& rhs1 = dynamic_cast<StochVector&>(rhs1_in);
  StochVector& rhs2 = dynamic_cast<StochVector&>(rhs2_in);
  StochVector& rhs3 = dynamic_cast<StochVector&>(rhs3_in);

  rhs.jointCopyFrom(rhs1, rhs2, rhs3);
}

void sLinsys::separateVars( OoqpVector& x_in, OoqpVector& y_in,
				     OoqpVector& z_in, OoqpVector& vars_in )
{
  StochVector& x    = dynamic_cast<StochVector&>(x_in);
  StochVector& y    = dynamic_cast<StochVector&>(y_in);
  StochVector& z    = dynamic_cast<StochVector&>(z_in);
  StochVector& vars = dynamic_cast<StochVector&>(vars_in);

  vars.jointCopyTo(x, y, z);
}




void sLinsys::factor(Data *prob_, Variables *vars)
{
  // the call to the the parent's method takes care of all necessary updates
  // to the KKT system (updating diagonals mainly). This is done reccursevely,
  // we don't have to worry about it anymore. 
  QpGenLinsys::factor(prob_, vars);

  // now DO THE LINEAR ALGEBRA!
  
  sData* prob = dynamic_cast<sData*>(prob_);
  // in order to avoid a call to QpGenLinsys::factor, call factor2 method.
  factor2(prob, vars);
}
 

/**
 * Computes U = Li\Gi^T.
 *        [ 0 0 0 ]
 * Gi^T = [ A 0 0 ]
 *        [ C 0 0]
 *
 * We have special structure here:
 *             [ 0 ]
 *   U   = Li\ [ A ] ,   U is (nx+my+mz)-by-(np)
 *             [ C ]
 *
 *   V = Di\U
 */
void sLinsys::computeU_V(sData *prob, 
			 DenseGenMatrix* U, DenseGenMatrix* V)
{
  U->scalarMult(0.0);
  V->scalarMult(0.0);

  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();

  int N, nxP;
  A.getSize(N, nxP); assert(N==locmy);

  N = locnx+locmy+locmz;
  SimpleVector uCol(N);

  for(int it=0; it<nxP; it++) {
    
    double* p = &uCol[0];
    for(int it1=0; it1<locnx; it1++) p[it1]=0.0;

    A.fromGetDense(0, it, &uCol[locnx],  1, locmy, 1);    
    C.fromGetDense(0, it, &uCol[locnx+locmy], 1, locmz, 1);

    solver->Lsolve(uCol);
    U->atPutDense(0, it, &uCol[0], 1, N, 1);

    solver->Dsolve(uCol);
    V->atPutDense(0, it, &uCol[0], 1, N, 1);

  }
}

void sLinsys::allocU(DenseGenMatrix ** U, int n0)
{
  int lines,cols;
  if(*U==NULL) {
    *U = new DenseGenMatrix(locnx+locmy+locmz, n0);
  } else {
    (*U)->getSize(lines,cols);
    
    if(lines!=locnx+locmy+locmz || cols != n0) {

      delete (*U);
      *U = new DenseGenMatrix(locnx+locmy+locmz, n0);
    }
  }
}

void sLinsys::allocV(DenseGenMatrix ** V, int n0)
{
  int lines,cols;
  if(*V==NULL)
    *V = new DenseGenMatrix(locnx+locmy+locmz, n0);
  else {
    (*V)->getSize(lines,cols);
    
    if(lines!=locnx+locmy+locmz || cols != n0) {
      
      delete (*V);
      *V = new DenseGenMatrix(locnx+locmy+locmz, n0);
    }
  }
}



/**
 *       [ 0 Ai^T Ci^T ]          [    ]
 * b0 -= [ 0   0   0   ] * Li\Di\ [ zi ]
 *       [ 0   0   0   ]          [    ]
 *
 * Changes only the first n0 entries of b0
 */
void sLinsys::addLnizi(sData *prob, OoqpVector& z0_, OoqpVector& zi_)
{
  SimpleVector& z0 = dynamic_cast<SimpleVector&>(z0_);
  SimpleVector& zi = dynamic_cast<SimpleVector&>(zi_);
  /*
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  static int c = 0;
  if (mype == 0 && c == 0) {
    printf("RHS IN\n");
    for (int i = 0; i < zi.length(); i++) {
      printf("%d %.10E\n",i,zi[i]);
    }
  }*/

  solver->Dsolve (zi);
  

  
  solver->Ltsolve(zi);


  /*if (mype == 0 && c == 0) {
    printf("SOL OUT\n");
    for (int i = 0; i < zi.length(); i++) {
      printf("%d %.10E\n",i,zi[i]);
    }
  }
  c++;*/

  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();

  //get n0= nx(parent)= #cols of A or C
  int dummy, n0;
  A.getSize(dummy, n0);

  // zi2 and zi3 are just references to fragements of zi
  SimpleVector zi2 (&zi[locnx], locmy);
  SimpleVector zi3 (&zi[locnx+locmy], locmz);
  // same for z01
  SimpleVector z01 (&z0[0], n0);
  
  A.transMult(1.0, z01, -1.0, zi2);
  C.transMult(1.0, z01, -1.0, zi3);


  if(stochNode->id()==-12) {
      printf("sleeping (why?)\n");
      //sleep(10);
      //usleep(10000000);
  }
}



void sLinsys::solveCompressed( OoqpVector& rhs_ )
{
  StochVector& rhs = dynamic_cast<StochVector&>(rhs_);
  //!log
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
  

  //!opt - mix execution of the following calls

  Lsolve (data, rhs); 
  //sleep(rank);
  //if(rank>=0)
  //{printf("-----------------\n");rhs.writeToStream(cout);printf("~~~~~~~~~~~~~~~~~~~~~~~~~~\n");}
  //!log

  Dsolve(data,rhs);
  //sleep(rank);
  //if(rank>=0)
  //{printf("-----------------\n");rhs.writeToStream(cout);printf("~~~~~~~~~~~~~~~~~~~~~~~~~~\n");}

 
  Ltsolve(data,rhs);
  //if(rank==0) {
    //printf("Lsolve\n");   rhs.writeToStream(cout);
  //}
}


/*
 *  y = alpha*Lni^T x + beta*y
 *
 *                       ( [ 0 0 0 ]     )
 *  y = beta*y + Di\Li\ (  [ A 0 0 ] * x )
 *                      (  [ C 0 0 ]    )
 */
void sLinsys::LniTransMult(sData *prob, 
				    SimpleVector& y, 
				    double alpha, SimpleVector& x)
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();

  int N, nx0;

  //get nx(parent) from the number of cols of A (or C). Use N as dummy
  A.getSize(N,nx0);
  // a mild assert
  assert(nx0 <= x.length());

  N = locnx+locmy+locmz;
  assert( y.length() == N);
  
  //!memopt
  SimpleVector LniTx(N);

  // shortcuts
  SimpleVector x1(&x[0], nx0);
  SimpleVector LniTx1(&LniTx[0], locnx);
  SimpleVector LniTx2(&LniTx[locnx], locmy);
  SimpleVector LniTx3(&LniTx[locnx+locmy], locmz);
  
  LniTx1.setToZero();
  A.mult(0.0, LniTx2, 1.0, x1);
  C.mult(0.0, LniTx3, 1.0, x1);
 
  solver->Lsolve(LniTx); 
  solver->Dsolve(LniTx);

  //!execopt : does LniTx have the first nx entries zero ?
  y.axpy(alpha,LniTx); 
 
}

/**
 * Computes U = Gi * inv(H_i) * Gi^T.
 *        [ 0 0 0 ]
 * Gi^T = [ A 0 0 ]
 *        [ C 0 0]
 *
 */

void sLinsys::addTermToDenseSchurCompl(sData *prob, 
				       DenseSymMatrix& SC) 
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();

  int N, nxP, NP;
  A.getSize(N, nxP); assert(N==locmy);
  NP = SC.size(); assert(NP>=nxP);

  if(nxP==-1) C.getSize(N,nxP);
  if(nxP==-1) nxP = NP;
  N = locnx+locmy+locmz;
	
/*
	const int blocksize = 20;
	for (int it=0; it<nxP; it += blocksize) {
		DenseGenMatrix out(SC[it], blocksize, NP);
		addColsToDenseSchurCompl(prob, out, it, MIN(it+blocksize,nxP));
	}
*/
		
  SimpleVector col(N);

  for(int it=0; it<nxP; it++) {
    
    double* pcol = &col[0];
    for(int it1=0; it1<locnx; it1++) pcol[it1]=0.0;

    A.fromGetDense(0, it, &col[locnx],  1, locmy, 1);    
    C.fromGetDense(0, it, &col[locnx+locmy], 1, locmz, 1);

    solver->solve(col);

    //here we have colGi = inv(H_i)* it-th col of Gi^t
    
    //now do colSC = Gi * inv(H_i)* it-th col of Gi^t
    //for(int it1=0; it1<locnx; it1++) pcol[it1]=0.0;
    
    // SC+=At*y
    //A.transMult( 1.0, &SC[0][it],       NP,  
		//-1.0, &col[locnx],       1);
    // SC+=Ct*z
    //C.transMult( 1.0, &SC[0][it], NP,
		//-1.0, &col[locnx+locmy], 1);

    // SC+=At*y
    A.transMult( 1.0, &SC[it][0],       1,  
		-1.0, &col[locnx],       1);
    // SC+=Ct*z
    C.transMult( 1.0, &SC[it][0], 1,
		-1.0, &col[locnx+locmy], 1);

    //SimpleVector a(&col[locnx+locmy], locmz);
    //a.writeToStream(cout);assert(false);
  }

}

#include <set>
#include <algorithm>

// we load the calculated columns into rows of out
// to match the column-major scalapack format

void sLinsys::addColsToDenseSchurCompl(sData *prob, 
				       DenseGenMatrix& out, 
				       int startcol, int endcol) 
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();


	int ncols = endcol-startcol;
  int N, nxP, NP, ncols_t, N_out;
  A.getSize(N, nxP); assert(N==locmy);
  out.getSize(ncols_t, N_out); 
  assert(N_out == nxP);
  assert(endcol <= nxP &&  ncols_t >= ncols);

  if(nxP==-1) C.getSize(N,nxP);
  if(nxP==-1) nxP = NP;

  N = locnx+locmy+locmz;
	DenseGenMatrix cols(ncols,N);
  bool allzero = true;
  memset(cols[0],0,N*ncols*sizeof(double));
  A.getStorage()->fromGetColBlock(startcol, &cols[0][locnx], N, endcol-startcol, allzero);
  C.getStorage()->fromGetColBlock(startcol, &cols[0][locnx+locmy], N, endcol-startcol, allzero);



	//int mype; MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	//printf("solving with multiple RHS %d \n", mype);	
	solver->solve(cols);
	//printf("done solving %d \n", mype);

	
	const int blocksize = 20;

	for (int it=0; it < ncols; it += blocksize) {
		int end = MIN(it+blocksize,ncols);
		int numcols = end-it;
		
		// SC-=At*y
		A.getStorage()->transMultMat( 1.0, out[it], numcols, N_out,  
		-1.0, &cols[it][locnx], N);
		// SC-=Ct*z
		C.getStorage()->transMultMat( 1.0, out[it], numcols, N_out,
		-1.0, &cols[it][locnx+locmy], N);
	}
	
//for (int it=startcol; it<endcol; it++) {
//  double* pcol = cols[it-startcol];
// 	A.transMult( 1.0, out[it-startcol], 1,  
//		-1.0, &pcol[locnx],       1);
//    // SC-=Ct*z
//   C.transMult( 1.0, out[it-startcol], 1,
//		-1.0, &pcol[locnx+locmy], 1);
//
//	}

	//printf("did transmult %d \n", mype);

}

/*
void sLinsys::addColsToDenseSchurCompl(sData *prob, 
				       DenseGenMatrix& out, 
				       int startcol, int endcol) 
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();

  int N, nxP, NP, ncols;
  A.getSize(N, nxP); assert(N==locmy);
  out.getSize(ncols, N); assert(N == nxP);
  assert(endcol <= nxP &&  ncols >= endcol - startcol);

  if(nxP==-1) C.getSize(N,nxP);
  if(nxP==-1) nxP = NP;

  N = locnx+locmy+locmz;
  SimpleVector col(N);

  for(int it=startcol; it<endcol; it++) {
    
    double* pcol = &col[0];
    for(int it1=0; it1<locnx; it1++) pcol[it1]=0.0;

    A.fromGetDense(0, it, &col[locnx],  1, locmy, 1);    
    C.fromGetDense(0, it, &col[locnx+locmy], 1, locmz, 1);

    solver->solve(col);

    //here we have colGi = inv(H_i)* it-th col of Gi^t
    
    //now do colSC = Gi * inv(H_i)* it-th col of Gi^t
    //for(int it1=0; it1<locnx; it1++) pcol[it1]=0.0;

    // SC-=At*y
    A.transMult( 1.0, out[it-startcol], 1,  
		-1.0, &col[locnx],       1);
    // SC-=Ct*z
    C.transMult( 1.0, out[it-startcol], 1,
		-1.0, &col[locnx+locmy], 1);

    //SimpleVector a(&col[locnx+locmy], locmz);
    //a.writeToStream(cout);assert(false);
  }
}*/


// adds only lower triangular elements to out

void sLinsys::symAddColsToDenseSchurCompl(sData *prob, 
				       double *out, 
				       int startcol, int endcol) 
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();

  int N, nxP, NP;
  A.getSize(N, nxP); assert(N==locmy);
  //out.getSize(ncols, N); assert(N == nxP);
  assert(endcol <= nxP);

  if(nxP==-1) C.getSize(N,nxP);
  if(nxP==-1) {assert(false); nxP = NP;} //petra - found that NP may be unitialized; initialized NP (to remove the compile warning) but added an assert

  N = locnx+locmy+locmz;

  const int BLOCKSIZE = 40;
  
  DenseGenMatrix cols(BLOCKSIZE,N);
  int outi = 0;
  
  
  // get list of completely zero columns
  /*
  static int go = true;
    if (go) {
    set<int> zeros;
    for (int i = 0; i < nxP; i++) {
      zeros.insert(i);
    }
    for (int i = 0; i < A.getStorage()->len; i++) {
      zeros.erase(A.getStorage()->jcolM[i]);
    }
    for (int i = 0; i < C.getStorage()->len; i++) {
      zeros.erase(C.getStorage()->jcolM[i]);
    }
    cout << "Zero columns: (of " << nxP << ")\n";
    for (set<int>::iterator it=zeros.begin(); it != zeros.end(); it++) {
      cout << *it << " ";
    }
    cout << endl;
    go = false;
  }*/

  for (int col = startcol; col < endcol; col += BLOCKSIZE) {
    int ecol = MIN(col+BLOCKSIZE,endcol);
    int nbcols = ecol-col;
    
    
    memset(cols[0],0,BLOCKSIZE*N*sizeof(double));
    
    //for (int c = col; c < ecol; c++) {
    //    A.fromGetDense(0,c,&cols[c-col][locnx],1,locmy,1);
    //    C.fromGetDense(0,c,&cols[c-col][locnx+locmy],1,locmz,1);
    //}
    bool allzero = true;
    
    A.getStorage()->fromGetColBlock(col, &cols[0][locnx], N, nbcols, allzero);
    C.getStorage()->fromGetColBlock(col, &cols[0][locnx+locmy], N, nbcols, allzero);
    
    if (!allzero) {
      solver->solve(cols);
      
      
      A.getStorage()->transMultMatLower(out+outi, nbcols, col,
					-1.0, &cols[0][locnx], N);
      C.getStorage()->transMultMatLower(out+outi, nbcols, col,
					-1.0, &cols[0][locnx+locmy], N);
    }
    for (int c = col; c < ecol; c++) {
      outi += nxP-c;
    }
    //
    //for (int c = col; c < ecol; c++) {
    //  A.getStorage()->transMultLower(1.0, out+outi,
    //       -1.0, &cols[c-col][locnx],c);
    //  C.getStorage()->transMultLower(1.0, out+outi,
    //       -1.0, &cols[c-col][locnx+locmy],c);
    //  outi += nxP - c;
    //}
  }
}

/*

void sLinsys::symAddColsToDenseSchurCompl(sData *prob, 
				       double *out, 
				       int startcol, int endcol) 
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();

  int N, nxP, NP, ncols;
  A.getSize(N, nxP); assert(N==locmy);
  N = nxP;
  //out.getSize(ncols, N); assert(N == nxP);
  assert(endcol <= nxP);

  if(nxP==-1) C.getSize(N,nxP);
  if(nxP==-1) nxP = NP;

  N = locnx+locmy+locmz;
  SimpleVector col(N);
  int outi = 0;

  for(int it=startcol; it<endcol; it++) {
    
    double* pcol = &col[0];
    for(int it1=0; it1<locnx; it1++) pcol[it1]=0.0;

    A.fromGetDense(0, it, &col[locnx],  1, locmy, 1);    
    C.fromGetDense(0, it, &col[locnx+locmy], 1, locmz, 1);

    solver->solve(col);

    //here we have colGi = inv(H_i)* it-th col of Gi^t
    
    //now do colSC = Gi * inv(H_i)* it-th col of Gi^t
    //for(int it1=0; it1<locnx; it1++) pcol[it1]=0.0;
    
    // SC-=At*y
    A.getStorage()->transMultLower( 1.0, out+outi,
                        -1.0, &col[locnx], it);
    // SC-=Ct*z
    C.getStorage()->transMultLower( 1.0, out+outi, 
                        -1.0, &col[locnx+locmy], it);

    int nelts = nxP-it;
    outi += nelts;
  }
}*/
