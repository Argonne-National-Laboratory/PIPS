
/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "sLinsys.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
#include "sVars.h"

#include "SparseLinearAlgebraPackage.h"

#include "NlpGenLinsys.h"

#include "NlpGenSparseLinsys.h"

#include "RegularizationAlg.h"

#ifndef MIN
#define MIN(a,b) ((a > b) ? b : a)
#endif

extern int gOuterIterRefin;

extern int gOuterSolve;
extern int separateHandDiag;

sLinsys::sLinsys(sFactory* factory_, sData* prob)
  : NlpGenLinsys(), kkt(NULL), solver(NULL), isActive(true)
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

  //cout << "sLinsys1: nxupp=" << nxupp << " nxlow=" << nxlow << "  sum=" << (nxlow+nxupp) << endl;

  //if( nxupp + nxlow > 0 ) {
  dd      	 = factory_->tree->newPrimalVector();
  assert(dd!=NULL);
  
  dq      	 = factory_->tree->newPrimalVector();
  temp_diagX = factory_->tree->newPrimalVector();
  
  prob->getDiagonalOfQ( *dq );
    //}
  nomegaInv   = factory_->tree->newDualZVector();

  if(gOuterSolve<3) {
  	rhs         =  factory_->tree->newRhs(); ;
  }
  else if(gOuterSolve>=3 && separateHandDiag==1){
  	rhs         =  factory_->tree->newRhsXSYZ();
	additiveDiag =  factory_->tree->newRhsXSYZ();
  }  
  else if(gOuterSolve>=3 && separateHandDiag==0){
  	rhs         =  factory_->tree->newRhsXSYZ();
  }else{
	assert(0);
  }

  temp_diagS =  factory_->tree->newDualZVector() ;
  temp_diagZ =	factory_->tree->newDualZVector() ;

  temp_diagY =  factory_->tree->newDualYVector() ;



  useRefs=0;
  data = prob;
  stochNode = prob->stochNode;
}

sLinsys::sLinsys(sFactory* factory_,
		 sData* prob,				    
		 OoqpVector* dd_, 
		 OoqpVector* dq_,
		 OoqpVector* nomegaInv_,
		 OoqpVector* rhs_,
		 OoqpVector* additiveDiag_)
  : NlpGenLinsys(), kkt(NULL), solver(NULL), isActive(true)
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

  //cout << "sLinsys2: nxupp=" << nxupp << " nxlow=" << nxlow << endl;

  //if( nxupp + nxlow > 0 ) {
  dd= dd_;
  dq = dq_;
    //}
  nomegaInv = nomegaInv_;
  rhs = rhs_;
  additiveDiag = additiveDiag_;

//  temp_diagX =	factory_->tree->newPrimalVector();
//  temp_diagS =  factory_->tree->newDualZVector() ;
//  temp_diagZ =  factory_->tree->newDualZVector() ;
//  temp_diagY =  factory_->tree->newDualYVector() ;

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

void sLinsys::joinRHSXSYZ( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
		  OoqpVector& rhs2_in, OoqpVector& rhs3_in, OoqpVector& rhs4_in )
{
  StochVector& rhs  = dynamic_cast<StochVector&>(rhs_in);
  StochVector& rhs1 = dynamic_cast<StochVector&>(rhs1_in);
  StochVector& rhs2 = dynamic_cast<StochVector&>(rhs2_in);
  StochVector& rhs3 = dynamic_cast<StochVector&>(rhs3_in);
  StochVector& rhs4 = dynamic_cast<StochVector&>(rhs4_in);

  rhs.jointCopyFromXSYZ(rhs1, rhs2, rhs3,rhs4);
}

void sLinsys::separateVarsXSYZ( OoqpVector& x_in, OoqpVector& s_in, 
		  OoqpVector& y_in, OoqpVector& z_in, OoqpVector& vars_in)
{
  StochVector& x    = dynamic_cast<StochVector&>(x_in);
  StochVector& s    = dynamic_cast<StochVector&>(s_in);
  StochVector& y    = dynamic_cast<StochVector&>(y_in);
  StochVector& z    = dynamic_cast<StochVector&>(z_in);
  StochVector& vars = dynamic_cast<StochVector&>(vars_in);

  vars.jointCopyToXSYZ(x, s, y, z);
}


void sLinsys::factor(Data *prob_, Variables *vars_in,RegularizationAlg *RegInfo)
{
  sData* prob = dynamic_cast<sData*>(prob_);
  sVars* vars = dynamic_cast<sVars*>(vars_in);
  
  bool skipUpdateReg=false;
  long long Num_NegEVal=-1;
  double priReg=0.0,dualReg=0.0;

#ifdef TIMING
  double tTot=MPI_Wtime();
#endif
  // DoEvalReg =   1 -> when factorizing the matrix, add regularizations to correct inertia and singularity 
  //				  	(we ONLY calling this routine once in IBR)
  //			  2 -> when factorizing the matrix, add regularizations to correct singularity only 
  //					(this is always the 1st call of this routine when IFR is used)
  //			  0 -> when factorizing the matrix, force to use primal regularizaion. called iff xWx tests fail  
  //					(the other calls of this routine when IFR is used, now matrix is nonsingular for sure)
  if(RegInfo->DoEvalReg >= 1){
  	RegInfo->newLinearSystem();

	if(RegInfo->ForceReg)
	  NlpGenLinsys::factorNoMatChange(prob_, vars, RegInfo);
	else
      NlpGenLinsys::factorNoMatChange(prob_, vars, NULL);	

    // now DO THE LINEAR ALGEBRA!
    // in order to avoid a call to NlpGenLinsys::factor, call factor2 method.
    Num_NegEVal = factor2(prob, vars);

	long long gbMy = prob->getGlobalMy();
	long long gbMz = prob->getGlobalMz();

	// check if matrix is singular
	if(Num_NegEVal < 0 || (Num_NegEVal < gbMy + gbMz && RegInfo->DoEvalReg == 1) )
	  RegInfo->MatrixSingular = 1;
	else
	  RegInfo->MatrixSingular = 0;	

	// skip update regularizaion if: 	1) have correct inertia and mat is nonsingular 
	//						OR 	2) mat is nonsingular and we will do inertia-free test later	
	if( (RegInfo->DoEvalReg==1 && Num_NegEVal == gbMy + gbMz) || (RegInfo->DoEvalReg == 2 && Num_NegEVal != -1)){
	  skipUpdateReg=true;
	}
  }

  // update regularization
  while( !skipUpdateReg ){
	RegInfo->computeRegularization(priReg,dualReg,prob->currMu);
	
    NlpGenLinsys::factorNoMatChange(prob_, vars, RegInfo);
	Num_NegEVal=(long long)factor2(prob, vars);

	// check if matrix is singular
	if(Num_NegEVal < 0)
	  RegInfo->MatrixSingular = 1;
	else
	  RegInfo->MatrixSingular = 0;

	// skip update regularizaion if: 	1) have correct inertia and mat is nonsingular 
	//						OR 	2) mat is nonsingular and we will do inertia-free test later	
	//						OR  	3) we are doing inertia-free test now (mat is defenitly nonsingular)
	if( (RegInfo->DoEvalReg == 1 && Num_NegEVal == prob->getGlobalMy() + prob->getGlobalMz())
		|| (RegInfo->DoEvalReg == 2 && Num_NegEVal != -1) || RegInfo->DoEvalReg == 0 )
	{
	  skipUpdateReg = true;
	}  	  
  }  

#ifdef TIMING
  tTot = MPI_Wtime()-tTot;
  MPI_Barrier(MPI_COMM_WORLD);
  int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if(0==myRank)
      cout << "Outer fact. total time " << tTot << endl;
#endif
}


void sLinsys::factor(Data *prob_, Variables *vars)
{
#ifdef TIMING
  double tTot=MPI_Wtime();
#endif
  // the call to the the parent's method takes care of all necessary updates
  // to the KKT system (updating diagonals mainly). This is done reccursevely,
  // we don't have to worry about it anymore. 
  NlpGenLinsys::factorNoMatChange(prob_, vars, NULL);

  // now DO THE LINEAR ALGEBRA!
  
  sData* prob = dynamic_cast<sData*>(prob_);
  // in order to avoid a call to NlpGenLinsys::factor, call factor2 method.
  factor2(prob, vars);

#ifdef TIMING
  tTot = MPI_Wtime()-tTot;
  MPI_Barrier(MPI_COMM_WORLD);
  int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  //if(128*(myRank/128)==0)
  if(0==myRank)
      cout << "Outer fact. total time " << tTot << endl;
#endif
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
  assert(gOuterSolve<3);

  U->scalarMult(0.0);
  V->scalarMult(0.0);
  assert(false); //need code to deal with cross Hessian term
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
	assert(gOuterSolve<3 || separateHandDiag==0);

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
	assert(gOuterSolve<3 );

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
 *       [ R^i^T Ai^T Ci^T ]          [    ]
 * z0 -= [ 0      0   0   ] * Li\Di\ [ zi ]
 *       [ 0      0   0   ]          [    ]
 *
 * 
 * if gOuterSolve>=3 && separateHandDiag==0
 *       [ R^i^T  0  Ai^T Ci^T ]          [    ]
 * z0 -= [ 0      0  0   0   ] * Li\Di\ [ zi ]
 *       [ 0         0  0   0   ]          [    ]
 *
 * 
 */
void sLinsys::addLnizi(sData *prob, OoqpVector& z0_, OoqpVector& zi_)
{
  SimpleVector& z0 = dynamic_cast<SimpleVector&>(z0_);
  SimpleVector& zi = dynamic_cast<SimpleVector&>(zi_);

  solver->Dsolve (zi);  
  solver->Ltsolve(zi);

  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& R = prob->getLocalCrossHessian();
  SparseGenMatrix& E = prob->getLocalE();
  SparseGenMatrix& F = prob->getLocalF();
  int mle = prob->getmle();
  int mli = prob->getmli();

  //get n0= nx(parent)= #cols of A or C
  int dummy, n0;
  A.getSize(dummy, n0);

  if(gOuterSolve>=3 ){
  	int locns = locmz;
    // zi2 and zi3 are just references to fragments of zi
    SimpleVector zi1 (&zi[0],           locnx);
    SimpleVector zi2 (&zi[locnx],       locns);
    SimpleVector zi3 (&zi[locnx+locns], locmy);
    SimpleVector zi4 (&zi[locnx+locns+locmy], locmz);
    // same for z01 (only the first n0 entries in the output z0 are computed)
    SimpleVector z01 (&z0[0], n0);

    R.transMult(1.0, z01, -1.0, zi1);
    A.transMult(1.0, z01, -1.0, zi3);
    C.transMult(1.0, z01, -1.0, zi4);
    
    if (mle > 0)
      {
	int nx0, my0, mz0;
	stochNode->get_FistStageSize(nx0, my0,mz0);
	SimpleVector z0E (&z0[nx0+mz0+my0-mle], mle);
	E.mult(1.0, z0E, -1.0, zi1);
      }
    if (mli > 0)
      {
	int nx0, my0, mz0;
	stochNode->get_FistStageSize(nx0, my0,mz0);
	SimpleVector z0F (&z0[nx0+mz0+my0+mz0-mli], mli);
	F.mult(1.0, z0F, -1.0, zi1);
      }
#if 0
	A.printMatrixInMatlab("bord_A");
	C.printMatrixInMatlab("bord_C");
	R.printMatrixInMatlab("bord_R");
#endif
  }
  else{
    // zi2 and zi3 are just references to fragments of zi
    SimpleVector zi1 (&zi[0],           locnx);
    SimpleVector zi2 (&zi[locnx],       locmy);
    SimpleVector zi3 (&zi[locnx+locmy], locmz);
    // same for z01 (only the first n0 entries in the output z0 are computed)
    SimpleVector z01 (&z0[0], n0);

    R.transMult(1.0, z01, -1.0, zi1);
    A.transMult(1.0, z01, -1.0, zi2);
    C.transMult(1.0, z01, -1.0, zi3);
  }
}



void sLinsys::solveCompressed( OoqpVector& rhs_ )
{
  StochVector& rhs = dynamic_cast<StochVector&>(rhs_);
#ifdef TIMING
  //double tTot=MPI_Wtime();
#endif
  Lsolve (data,rhs); 
  Dsolve (data,rhs);
  Ltsolve(data,rhs);
#ifdef TIMING
  //cout << "SolveCompressed took: " << (MPI_Wtime()-tTot) << endl;
#endif
}


/*
 *  y = alpha*Lni^T x + beta*y
 *
 *                       ( [ R 0 0 ]     )
 *  y = beta*y + Di\Li\ (  [ A 0 0 ] * x )
 *                      (  [ C 0 0 ]    )
*
* if gOuterSolve>=3 && separateHandDiag==0 :
*
*
*						( [ R 0 0 ] 	)
*						( [ 0 0 0 ] 	)
*  y = beta*y + Di\Li\ 		( [ A 0 0 ] * x	)
*					  	( [ C 0 0 ]     	)
* 
*/

void sLinsys::LniTransMult(sData *prob, 
			   OoqpVector& y_, 
			   double alpha, OoqpVector& x_)
{
  SimpleVector& y = dynamic_cast<SimpleVector&>(y_);
  SimpleVector& x = dynamic_cast<SimpleVector&>(x_);


  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& R = prob->getLocalCrossHessian();
  SparseGenMatrix& E = prob->getLocalE();
  SparseGenMatrix& F = prob->getLocalF();
  int mle = prob->getmle();
  int mli = prob->getmli();

  int N, nx0;
  int locns = locmz;

  //get nx(parent) from the number of cols of A (or C). Use N as dummy
  A.getSize(N,nx0);
  // a mild assert
  assert(nx0 <= x.length());

  if(gOuterSolve>=3 ) 
	N = locnx+locns+locmy+locmz;
  else
    N = locnx+locmy+locmz;

  assert( y.length() == N);

  //!memopt
  SimpleVector LniTx(N);

  SimpleVector x1(&x[0], nx0);

  if(gOuterSolve>=3 ) {
    // shortcuts
    SimpleVector LniTx1(&LniTx[0], locnx);
    SimpleVector LniTx2(&LniTx[locnx], locns);
    SimpleVector LniTx3(&LniTx[locnx+locns], locmy);
    SimpleVector LniTx4(&LniTx[locnx+locns+locmy], locmz);
	
    LniTx1.setToZero();
	LniTx2.setToZero();
    R.mult(0.0, LniTx1, 1.0, x1);
    A.mult(0.0, LniTx3, 1.0, x1);
    C.mult(0.0, LniTx4, 1.0, x1);
    if (mle > 0)
      {
	int nx0, my0, mz0;
        stochNode->get_FistStageSize(nx0, my0,mz0);
        SimpleVector xE (&x[nx0+mz0+my0-mle], mle);
        E.transMult(1.0, LniTx1, 1.0, xE);
      }
    if (mli > 0)
      {
        int nx0, my0, mz0;
        stochNode->get_FistStageSize(nx0, my0,mz0);
        SimpleVector xF (&x[nx0+mz0+my0+mz0-mli], mli);
        F.transMult(1.0, LniTx1, 1.0, xF);
      }

  }else{
    // shortcuts
//    SimpleVector x1(&x[0], nx0);
    SimpleVector LniTx1(&LniTx[0], locnx);
    SimpleVector LniTx2(&LniTx[locnx], locmy);
    SimpleVector LniTx3(&LniTx[locnx+locmy], locmz);
  
    LniTx1.setToZero();
    R.mult(0.0, LniTx1, 1.0, x1);
    A.mult(0.0, LniTx2, 1.0, x1);
    C.mult(0.0, LniTx3, 1.0, x1);
  }
  
  solver->Lsolve(LniTx); 
  solver->Dsolve(LniTx);
  solver->Ltsolve(LniTx);
  
  y.axpy(alpha,LniTx); 
 
}
/*
 * Computes res += [0 A^T C^T ]*inv(KKT)*[0;A;C] x
 */

void sLinsys::addTermToSchurResidual(sData* prob, 
				     SimpleVector& res, 
				     SimpleVector& x)
{
	assert(gOuterSolve<3 );

  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& R = prob->getLocalCrossHessian();

  int nxP, aux;
  A.getSize(aux,nxP); assert(aux==locmy);
  C.getSize(aux,nxP); assert(aux==locmz);
  R.getSize(aux,nxP); assert(aux==locnx);
  assert(nxP==x.length());
  int N=locnx+locmy+locmz;
  SimpleVector y(N);
  //y.setToZero();

  R.mult( 0.0,&y[0],1,           1.0,&x[0],1);
  A.mult( 0.0,&y[locnx],1,       1.0,&x[0],1);
  C.mult( 0.0,&y[locnx+locmy],1, 1.0,&x[0],1);
  //cout << "4 - y norm:" << y.twonorm() << endl;
  //printf("%g  %g  %g  %g\n", y[locnx+locmy+0], y[locnx+locmy+1], y[locnx+locmy+2], y[locnx+locmy+3]);
  solver->solve(y);

  R.transMult(1.0,&res[0],1, 1.0,&y[0],1);
  A.transMult(1.0,&res[0],1, 1.0,&y[locnx],1);
  C.transMult(1.0,&res[0],1, 1.0,&y[locnx+locmy],1);
}

#include "PardisoSolver.h"
/**
 * Computes U = Gi * inv(H_i) * Gi^T.
 *        [ R 0 0 ]
 * Gi^T = [ A 0 0 ]
 *        [ C 0 0]
 *
 * A and C are the recourse eq. and ineq. matrices, R is the cross
 * Hessian term.
 *
 * if gOuterSolve>=3 && separateHandDiag==0
 *        		[ R 0 0 ]
 *        		[ 0 0 0 ]   <--- for the Ineq Cons slack var
 * Gi^T = 	[ A 0 0 ]
 *        		[ C 0 0]
 *
 * A and C are the recourse eq. and ineq. matrices, R is the cross
 * Hessian term.
 */

void sLinsys::addTermToDenseSchurCompl(sData *prob, 
				       DenseSymMatrix& SC) 
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& R = prob->getLocalCrossHessian();
  

  int N, nxP, NP,mR,nR;
  int locns = locmz;
  
  int mle = prob->getmle();
  int mli = prob->getmli();
  SparseGenMatrix& E = prob->getLocalE();
  SparseGenMatrix& F = prob->getLocalF();

  SparseGenMatrix ET;
  SparseGenMatrix FT;
  ET.transCopyof(E);
  FT.transCopyof(F);

  int nx0, my0, mz0;
  stochNode->get_FistStageSize(nx0, my0,mz0);
   
  A.getSize(N, nxP); assert(N==locmy);
  NP = SC.size(); assert(NP>=nxP);

  if(nxP==-1) C.getSize(N,nxP);
  if(nxP==-1) nxP = NP;
  if(gOuterSolve>=3 ) 
	N = locnx+locns+locmy+locmz;
  else
    N = locnx+locmy+locmz;


  int blocksize = 64;
  DenseGenMatrix cols(blocksize,N);

  bool ispardiso=false;
  PardisoSolver* pardisoSlv=NULL;
//  pardisoSlv = dynamic_cast<PardisoSolver*>(solver);
  int* colSparsity=NULL;
  if(pardisoSlv) {
    ispardiso=true;
    colSparsity=new int[N];
    //blocksize=32;
  }
  
  for (int it=0; it < nxP; it += blocksize) {
    int start=it;
    int end = MIN(it+blocksize,nxP);
    int numcols = end-start;
    cols.getStorageRef().m = numcols; // avoid extra solves

    bool allzero = true;
    memset(&cols[0][0],0,N*blocksize*sizeof(double));

    if(ispardiso) {
      for(int i=0; i<N; i++) colSparsity[i]=0;
	  if(gOuterSolve>=3 ) {
		R.getStorageRef().fromGetColBlock(start, &cols[0][0], N, numcols, colSparsity, allzero);
		A.getStorageRef().fromGetColBlock(start, &cols[0][locnx+locns], N, numcols, colSparsity, allzero);
		C.getStorageRef().fromGetColBlock(start, &cols[0][locnx+locns+locmy], N, numcols, colSparsity, allzero);
	  }
	  else{
		R.getStorageRef().fromGetColBlock(start, &cols[0][0], N, numcols, colSparsity, allzero);
	    A.getStorageRef().fromGetColBlock(start, &cols[0][locnx], N, numcols, colSparsity, allzero);
		C.getStorageRef().fromGetColBlock(start, &cols[0][locnx+locmy], N, numcols, colSparsity, allzero);
	  }
    } else {
	  if(gOuterSolve>=3 ) {
		R.getStorageRef().fromGetColBlock(start, &cols[0][0], N, numcols, allzero);
		A.getStorageRef().fromGetColBlock(start, &cols[0][locnx+locns], N, numcols, allzero);
		C.getStorageRef().fromGetColBlock(start, &cols[0][locnx+locns+locmy], N, numcols, allzero);
	  }
	  else{
		R.getStorageRef().fromGetColBlock(start, &cols[0][0], N, numcols, allzero);
		A.getStorageRef().fromGetColBlock(start, &cols[0][locnx], N, numcols, allzero);
		C.getStorageRef().fromGetColBlock(start, &cols[0][locnx+locmy], N, numcols, allzero);
	  }    
    }

    if(!allzero) {
      
      if(ispardiso)
		pardisoSlv->solve(cols,colSparsity);
      else 
		solver->solve(cols);

      if(gOuterSolve>=3 ) {
        R.getStorageRef().transMultMat( 1.0, &(SC[0][start]), numcols, NP,  
       				      -1.0, &cols[0][0], N);
        A.getStorageRef().transMultMat( 1.0, &(SC[0][start]), numcols, NP,  
       				      -1.0, &cols[0][locnx+locns], N);	
        C.getStorageRef().transMultMat( 1.0, &(SC[0][start]), numcols, NP,
       				      -1.0, &cols[0][locnx+locns+locmy], N);
	if(mle>0)
          ET.getStorageRef().transMultMat( 1.0,  &(SC.getStorageRef().M[nx0+mz0+my0-mle][start]), numcols, NP, -1.0, &cols[0][0], N);
	if(mli>0)
	    FT.getStorageRef().transMultMat( 1.0, &(SC.getStorageRef().M[nx0+mz0+my0+mz0-mli][start]), numcols, NP,
                                   -1.0, &cols[0][0], N);
	  }
	  else{
	    R.getStorageRef().transMultMat( 1.0, &(SC[0][start]), numcols, NP,  
       				      -1.0, &cols[0][0], N);
	    A.getStorageRef().transMultMat( 1.0, &(SC[0][start]), numcols, NP,  
       				      -1.0, &cols[0][locnx], N);
	    C.getStorageRef().transMultMat( 1.0, &(SC[0][start]), numcols, NP,
       				      -1.0, &cols[0][locnx+locmy], N);
	  }
    } //end !allzero
  }

  for (int it=0; it < mle; it += blocksize)
  {
    int start=it;
    int end = MIN(it+blocksize,mle);
    int numcols = end-start;
    cols.getStorageRef().m = numcols; // avoid extra solves                                                                                                                          
    bool allzero = true;
    memset(&cols[0][0],0,N*blocksize*sizeof(double));

    if(ispardiso) 
    {
        for(int i=0; i<N; i++) colSparsity[i]=0;
	ET.getStorageRef().fromGetColBlock(start, &cols[0][0], N, numcols, colSparsity, allzero);
    }
    else 
    {
	ET.getStorageRef().fromGetColBlock(start, &cols[0][0], N, numcols, allzero);
    }
    if(!allzero) {
      if(ispardiso)
	pardisoSlv->solve(cols,colSparsity);
      else
	solver->solve(cols);
      if(gOuterSolve>=3 ) {
	R.getStorageRef().transMultMat( 1.0, &(SC[0][nx0+mz0+my0-mle+start]), numcols, NP,
					-1.0, &cols[0][0], N);
	A.getStorageRef().transMultMat( 1.0, &(SC[0][nx0+mz0+my0-mle+start]), numcols, NP,
					-1.0, &cols[0][locnx+locns], N);
	C.getStorageRef().transMultMat( 1.0, &(SC[0][nx0+mz0+my0-mle+start]), numcols, NP,
					-1.0, &cols[0][locnx+locns+locmy], N);
	if(mle>0)
	  ET.getStorageRef().transMultMat( 1.0,  &(SC[nx0+mz0+my0-mle][nx0+mz0+my0-mle+start]), numcols, NP, -1.0, &cols[0][0], N);
	if(mli>0)
	  FT.getStorageRef().transMultMat( 1.0,  &(SC[nx0+mz0+my0+mz0-mli][nx0+mz0+my0-mle+start]), numcols, NP, -1.0, &cols[0][0], N);

	/*
	std::cout<<"cols:  "<<std::endl;
        for(int i=0; i<numcols;i++)
          for(int j=0; j<N;j++)
	    std::cout<<"col "<<i<<"row "<<j<<"elt "<<cols[i][j]<<std::endl;

	std::cout<<"SC:  "<<std::endl;
        for(int i=0; i<NP;i++)
          for(int j=0; j<NP;j++)
	    std::cout<<"row "<<i<<"col "<<j<<"elt "<<SC.getStorageRef().M[i][j]<<std::endl;
	*/
      }
      else{
	assert(false && "not implemented");
      }
    } //end !allzero 
  }

  for (int it=0; it < mli; it += blocksize)
  {
      int start=it;
      int end = MIN(it+blocksize,mle);
      int numcols = end-start;
      cols.getStorageRef().m = numcols; // avoid extra solves
      bool allzero = true;
      memset(&cols[0][0],0,N*blocksize*sizeof(double));

      if(ispardiso)
      {
	  for(int i=0; i<N; i++) colSparsity[i]=0;
	  FT.getStorageRef().fromGetColBlock(start, &cols[0][0], N, numcols, colSparsity, allzero);
      }
      else
      {
	  FT.getStorageRef().fromGetColBlock(start, &cols[0][0], N, numcols, allzero);
      }
      if(!allzero) {
	if(ispardiso)
	  pardisoSlv->solve(cols,colSparsity);
	else
	  solver->solve(cols);
	if(gOuterSolve>=3 ) {
	  R.getStorageRef().transMultMat( 1.0, &(SC[0][nx0+mz0+my0+mz0-mli+start]), numcols, NP,
					  -1.0, &cols[0][0], N);
	  A.getStorageRef().transMultMat( 1.0, &(SC[0][nx0+mz0+my0+mz0-mli+start]), numcols, NP,
					  -1.0, &cols[0][locnx+locns], N);
	  C.getStorageRef().transMultMat( 1.0, &(SC[0][nx0+mz0+my0+mz0-mli+start]), numcols, NP,
					  -1.0, &cols[0][locnx+locns+locmy], N);
	  if(mle>0)
	    ET.getStorageRef().transMultMat( 1.0,  &(SC[nx0+mz0+my0-mle][nx0+mz0+my0+mz0-mli+start]), numcols, NP, -1.0, &cols[0][0], N);
	  if(mli>0)
	    FT.getStorageRef().transMultMat( 1.0,  &(SC[nx0+mz0+my0+mz0-mli][nx0+mz0+my0+mz0-mli+start]), numcols, NP, -1.0, &cols[0][0], N);
	}
	else{
	  assert(false && "not implemented");
	}
      } //end !allzero
    }

  if(ispardiso) delete[] colSparsity;
}


/* this is the original code that was doing one column at a time. */
/* 
void sLinsys::addTermToDenseSchurCompl(sData *prob, 
				       DenseSymMatrix& SC) 
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& R = prob->getLocalCrossHessian();

  int N, nxP, NP;
  A.getSize(N, nxP); assert(N==locmy);
  NP = SC.size(); assert(NP>=nxP);

  if(nxP==-1) C.getSize(N,nxP);
  if(nxP==-1) nxP = NP;
  N = locnx+locmy+locmz;

  SimpleVector col(N);

  for(int it=0; it<nxP; it++) {
    
    double* pcol = &col[0];
    for(int it1=0; it1<locnx; it1++) pcol[it1]=0.0;

    R.fromGetDense(0, it, &col[0],           1, locnx, 1);
    A.fromGetDense(0, it, &col[locnx],       1, locmy, 1);    
    C.fromGetDense(0, it, &col[locnx+locmy], 1, locmz, 1);

    solver->solve(col);

    //here we have colGi = inv(H_i)* it-th col of Gi^t
     //now do colSC = Gi * inv(H_i)* it-th col of Gi^t
 
    // SC+=R*x
    R.transMult( 1.0, &SC[it][0],     1,
		 -1.0, &col[0],      1);

    // SC+=At*y
    A.transMult( 1.0, &SC[it][0],   1,  
		-1.0, &col[locnx],  1);
    // SC+=Ct*z
    C.transMult( 1.0, &SC[it][0],   1,
		-1.0, &col[locnx+locmy], 1);

  }
}*/
 
#include <set>
#include <algorithm>

// we load the calculated columns into rows of out
// to match the column-major scalapack format

void sLinsys::addColsToDenseSchurCompl(sData *prob, 
				       DenseGenMatrix& out, 
				       int startcol, int endcol) 
{
	assert(gOuterSolve<3 );

  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& R = prob->getLocalCrossHessian();

  int ncols = endcol-startcol;
  int N, nxP, ncols_t, N_out;
  A.getSize(N, nxP); assert(N==locmy);
  out.getSize(ncols_t, N_out); 
  assert(N_out == nxP);
  assert(endcol <= nxP &&  ncols_t >= ncols);

  if(nxP==-1) C.getSize(N,nxP);
  //if(nxP==-1) nxP = NP;

  N = locnx+locmy+locmz;
  DenseGenMatrix cols(ncols,N);
  bool allzero = true;
  memset(cols[0],0,N*ncols*sizeof(double));

  R.getStorageRef().fromGetColBlock(startcol, &cols[0][0],
				    N, endcol-startcol, allzero);
  A.getStorageRef().fromGetColBlock(startcol, &cols[0][locnx], 
				    N, endcol-startcol, allzero);
  C.getStorageRef().fromGetColBlock(startcol, &cols[0][locnx+locmy], 
				    N, endcol-startcol, allzero);

  //int mype; MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  //printf("solving with multiple RHS %d \n", mype);	
  solver->solve(cols);
  //printf("done solving %d \n", mype);
  
  
  const int blocksize = 20;
  
  for (int it=0; it < ncols; it += blocksize) {
    int end = MIN(it+blocksize,ncols);
    int numcols = end-it;
    assert(false); //add Rt*x -- and test the code
    // SC-=At*y
    A.getStorageRef().transMultMat( 1.0, out[it], numcols, N_out,  
				  -1.0, &cols[it][locnx], N);
    // SC-=Ct*z
    C.getStorageRef().transMultMat( 1.0, out[it], numcols, N_out,
				  -1.0, &cols[it][locnx+locmy], N);
  }
  

}



// adds only lower triangular elements to out

void sLinsys::symAddColsToDenseSchurCompl(sData *prob, 
				       double *out, 
				       int startcol, int endcol) 
{
	assert(gOuterSolve<3 );

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
  

  for (int col = startcol; col < endcol; col += BLOCKSIZE) {
    int ecol = MIN(col+BLOCKSIZE,endcol);
    int nbcols = ecol-col;
    
    
    memset(cols[0],0,BLOCKSIZE*N*sizeof(double));
    
 
    bool allzero = true;
    assert(false); //! code needs to consider the cross Hessian; fixme
    A.getStorageRef().fromGetColBlock(col, &cols[0][locnx], N, nbcols, allzero);
    C.getStorageRef().fromGetColBlock(col, &cols[0][locnx+locmy], N, nbcols, allzero);
    
    if (!allzero) {
      solver->solve(cols);
      
      A.getStorageRef().transMultMatLower(out+outi, nbcols, col,
					-1.0, &cols[0][locnx], N);
      C.getStorageRef().transMultMatLower(out+outi, nbcols, col,
					  -1.0, &cols[0][locnx+locmy], N); 
    }
    for (int c = col; c < ecol; c++) {
      outi += nxP-c;
    }
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
    A.getStorageRef().transMultLower( 1.0, out+outi,
                        -1.0, &col[locnx], it);
    // SC-=Ct*z
    C.getStorageRef().transMultLower( 1.0, out+outi, 
                        -1.0, &col[locnx+locmy], it);

    int nelts = nxP-it;
    outi += nelts;
  }
}*/

void sLinsys::_backSolve(sData *prob, OoqpVector& ParSol_ , OoqpVector& Vec_, StochVector* End_Par_Pos_)
{
  assert("not done" &&0);
}

void sLinsys::_addTargetParsLnizi(sData *prob, OoqpVector& ParSol_ , OoqpVector& Vec_, OoqpVector* goal_Par)
{
  assert("not done" &&0);
}

void sLinsys::_setupColOfBordMat( sData *prob, OoqpVector* rhs_St, const int ColIDX, bool &allzero, 
									  const int aimlevel)
{
  assert("not done" &&0);
}

void sLinsys::_assembleSC( sData *prob, OoqpVector* rhs_St_in, const int ColIDX, 
										const int aimlevel, DenseSymMatrix& SC)
{
  assert("not done" &&0);
}


void sLinsys:: _setupYaddLniTx( sData *prob, OoqpVector& y_, 
			   double alpha, SimpleVector& x, const int aimlevel)
{
	assert("not done" &&0);
}

