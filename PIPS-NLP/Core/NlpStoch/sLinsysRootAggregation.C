/* PIPS-NLP                                                           	*
 * Author:  Nai-Yuan Chiang                                       	*
 * (C) 2015 Argonne National Laboratory. 				*/

#include "sLinsysRootAggregation.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
#include "sDummyLinsys.h"
#include "sLinsysLeaf.h"

#include "Ma57Solver.h"

#include "StochVector.h"
#include <cassert>

#include "NlpGenVars.h"


#ifdef STOCH_TESTING
extern double g_iterNumber;
extern double g_scenNum;
#endif

#include "constants.h"
extern PreCondInfo *preCond;

/*********************************************************************/
/************************** ROOT *************************************/
/*********************************************************************/

sLinsysRootAggregation::sLinsysRootAggregation(sFactory * factory_, sData * prob_)
  : sLinsysRoot(factory_, prob_)
{
  prob_->getLocalSizes(locnx, locmy, locmz);
  assert(locmz==0);

  n_col = preCond->n_col;  
  n_row = preCond->n_row;

  JNnz  = preCond->jac_nnz;
  QNnz  = preCond->hes_nnz;

  redDim = n_col + n_row;
  redNnz = JNnz + QNnz + n_col + n_row; // nnz for jac, hes and diag part
  
  // vectors for the reduced system
  redRhs 		= new SimpleVector(redDim);
  temp_ncol 	= new SimpleVector(n_col);
  temp_nrow 	= new SimpleVector(n_row);

  temp_ixlow 	= new SimpleVector(n_col);
  temp_ixupp 	= new SimpleVector(n_col);

  temp_x 		= new SimpleVector(n_col);
  temp_y 		= new SimpleVector(n_row);
  temp_gamma 	= new SimpleVector(n_col);
  temp_phi 		= new SimpleVector(n_col);
  temp_v 		= new SimpleVector(n_col);
  temp_w 		= new SimpleVector(n_col);

  rQ 			= new SimpleVector(n_col);
  rA 			= new SimpleVector(n_row);
  rv 			= new SimpleVector(n_col);
  rw 			= new SimpleVector(n_col);
  rgamma 		= new SimpleVector(n_col);
  rphi 			= new SimpleVector(n_col);

  temp_rhs_x 	= new SimpleVector(n_col);
  temp_rhs_y 	= new SimpleVector(n_row); 

  temp_xl 		= new SimpleVector(n_col);
  temp_xu		= new SimpleVector(n_col);


  kkt 	 = createKKT(prob_);
  solver = createSolver(prob_, kkt);

  // this is Schur Complement matrix, set it to dummy matrix
  CtDC =  new DenseSymMatrix(0);

  
  setIX = false;
}

sLinsysRootAggregation::sLinsysRootAggregation(sFactory* factory_,
			 sData* prob_,
			 OoqpVector* dd_, 
			 OoqpVector* dq_,
			 OoqpVector* nomegaInv_,
			 OoqpVector* rhs_, OoqpVector* additiveDiag)
  : sLinsysRoot(factory_, prob_, dd_, dq_, nomegaInv_, rhs_,additiveDiag)
{
  n_col = preCond->n_col;  
  n_row = preCond->n_row;

  JNnz = preCond->jac_nnz;
  QNnz = preCond->hes_nnz;

  redDim = n_col + n_row;
  redNnz = JNnz + QNnz + n_col + n_row; // nnz for jac, hes and diag part
  
  redRhs = new SimpleVector(redDim);
  kkt = createKKT(prob_);
  solver = createSolver(prob_, kkt);

  temp_ncol = new SimpleVector(n_col);
  temp_nrow = new SimpleVector(n_row);

  // this is Schur Complement matrix, set it to dummy matrix
  CtDC =  new DenseSymMatrix(0);
}

sLinsysRootAggregation::~sLinsysRootAggregation()
{ 
	if(MatJ) delete MatJ;
	if(MatQ) delete MatQ;

	delete CtDC;
	delete redRhs;
	delete temp_ncol;
	delete temp_nrow;

	delete temp_ixlow;
	delete temp_ixupp;	

	delete temp_x;
	delete temp_y;	
	delete temp_gamma;
	delete temp_phi;	
	delete temp_v;
	delete temp_w;	

	delete rQ;
	delete rA;
	delete rv;
	delete rw;	
	delete rgamma;
	delete rphi;	

	delete temp_rhs_x;
	delete temp_rhs_y;
	
	delete temp_xl;
	delete temp_xu;		
}


SymMatrix* 
sLinsysRootAggregation::createKKT(sData* prob)
{
  int info;

  // creat aggregation matrix
  SparseSymMatrix *Mat = new SparseSymMatrix( redDim, redNnz) ;

  MatJ = new SparseGenMatrix( n_row, n_col, JNnz );
  MatJ->putSparseTriple(preCond->jac_irow, JNnz, preCond->jac_jcol, preCond->jac_value, info );

  MatQ = new SparseSymMatrix( n_col, QNnz );
  MatQ->putSparseTriple(preCond->hes_irow, QNnz, preCond->hes_jcol, preCond->hes_value, info );

  
  SimpleVectorHandle v( new SimpleVector(redDim) );
  v->setToZero();
  Mat->setToDiagonal(*v);

  Mat->symAtPutSubmatrix( 0, 0, *MatQ, 0, 0, n_col, n_col );
  Mat->symAtPutSubmatrix( n_col, 0, *MatJ, 0, 0, n_row, n_col );

  redNnz = Mat->numberOfNonZeros();
  
  return Mat;
}


DoubleLinearSolver*
sLinsysRootAggregation::createSolver(sData* prob, SymMatrix* kktmat_)
{
  SparseSymMatrix* kktmat = dynamic_cast<SparseSymMatrix*>(kktmat_);
  return new Ma57Solver(kktmat);
}

// note that diagonal of the original full prob has been computed
int sLinsysRootAggregation::factor2(sData *prob, Variables *vars)
{
  int negEVal=0, tempNegEVal=0;
  int return_NegEval;  
  int matIsSingular=0,matIsSingularAllReduce;
  int mype; MPI_Comm_rank(mpiComm, &mype);

  SparseSymMatrix& kktsp = dynamic_cast<SparseSymMatrix&>(*kkt);
  
  DenseSymMatrix& kktd = dynamic_cast<DenseSymMatrix&>(*CtDC);

  // First tell children to factorize. 
  for(int c=0; c<children.size(); c++) {
    tempNegEVal = children[c]->factor2(prob->children[c], vars);
	if(tempNegEVal<0){
	  matIsSingular = 1; 
	  break;
	}else{
	  negEVal += tempNegEVal;
	}	
  }
  MPI_Allreduce(&matIsSingular, &matIsSingularAllReduce, 1, MPI_INT, MPI_SUM, mpiComm);
  MPI_Allreduce(&negEVal, &return_NegEval, 1, MPI_INT, MPI_SUM, mpiComm);

  // all the diag mat is nonsingular
  if(0==matIsSingularAllReduce){

	// generate aggregation matrix
	reduceKKT();
	
	// factorize aggregation matrix
	tempNegEVal = factorizeKKT();
	MPI_Bcast(&return_NegEval, 1, MPI_INT, 0, mpiComm);	 
	
	if(tempNegEVal<0){ 
	  return_NegEval = -1;
	}else{
	  return_NegEval += negEVal;
	}
  }  
 
  return return_NegEval;
}



// put diag into the aggregation sys
void sLinsysRootAggregation::reduceKKT()
{
  int myRank; MPI_Comm_rank(mpiComm,&myRank);
  if(0==myRank){
    SparseSymMatrix* kktsp = dynamic_cast<SparseSymMatrix*>(kkt); 
    StochVector* ddsp = dynamic_cast<StochVector*>(dd); 

    assert(n_col == dd->length());
    temp_ncol->setToZero();

    ddsp->copyIntoArrayWithIndex_AggVarCon(temp_ncol->elements(),NULL,preCond->n_col,true);
    kktsp->atPutDiagonal( 0, *temp_ncol ); 
  }
}


int sLinsysRootAggregation::factorizeKKT()
{
  int negEValTemp = 0;
  int myRank; MPI_Comm_rank(mpiComm,&myRank);
  if(0==myRank)
    negEValTemp = solver->matrixChanged();
  return negEValTemp;
}

void sLinsysRootAggregation::Lsolve(sData *prob, OoqpVector& rhs)
{
  StochVector& b = dynamic_cast<StochVector&>(rhs);
  assert(children.size() == b.children.size() );
  int myRank; MPI_Comm_rank(mpiComm, &myRank);

  // children compute their part
  for(size_t it=0; it<children.size(); it++) {
    children[it]->Lsolve(prob->children[it], *b.children[it]);  
  }

  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.vec);

  //this code actually works on a single CPU too :)
  if (iAmDistrib) {
    //only one process add b0
    if(myRank>0) {
      b0.setToZero();
    }
  }
 
}

void sLinsysRootAggregation::Dsolve( sData *prob, OoqpVector& rhs )
{
  StochVector& b = dynamic_cast<StochVector&>(rhs);
  int myRank; MPI_Comm_rank(mpiComm,&myRank);

  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.vec);

  if(0==myRank)
    solveReduced(prob, b0);
  
  // broadcast b0 from rank 0 to others
  MPI_Bcast(b0.elements(), b0.length(), MPI_DOUBLE, 0, mpiComm);
}

void sLinsysRootAggregation::Ltsolve( sData *prob, OoqpVector& rhs )
{
  StochVector& b   = dynamic_cast<StochVector&>(rhs);
  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.vec);

  SimpleVector& x0 = b0; //just another name, for clarity
  
  // Li^T\bi for each child i. The backsolve needs z0
  for(size_t it=0; it<children.size(); it++) {
    children[it]->Ltsolve2(prob->children[it], *b.children[it], x0);
  }
}


// solve aggregation system
void sLinsysRootAggregation::solveReduced( sData *prob, SimpleVector& sol_0)
{
  int myRank; MPI_Comm_rank(mpiComm, &myRank);
  SimpleVector& aggRHS = (*redRhs);

  assert(locmz==0); // no inequalities
  assert(locnx+locmy==sol_0.length());  
  assert(aggRHS.length() >= sol_0.length() && aggRHS.length()>0);

  // build aggregation rhs
  calcPreCondKKTResids(prob);
  computeReducedRhs();
  _joinRedRHS( aggRHS, *temp_rhs_x, *temp_rhs_y);

  // solve aggregation system
  solver->Dsolve(aggRHS);

  // reset 1st stage var
  _separateVars( aggRHS, *temp_rhs_x, *temp_rhs_y);
  _set1stStVar(aggRHS, sol_0);

}


void sLinsysRootAggregation::calcPreCondKKTResids(Data *prob_in)
{
  NlpGenVars * vars = (NlpGenVars *) preCond->curr_Iter;
  NlpGenData * prob = (NlpGenData *) prob_in;

  double componentNorm, norm=0.0, gap=0.0;

  //assemble x,y  
  assert(vars->x->length() >= preCond->n_col);
  assert(vars->y->length() >= preCond->n_row);
  assert(locmz==0);

  temp_x->setToZero();
  vars->x->copyIntoArrayWithIndex_AggVarCon(temp_x->elements(),NULL,preCond->n_col,true);  
  temp_y->setToZero();
  vars->y->copyIntoArrayWithIndex_AggVarCon(temp_y->elements(),NULL,preCond->n_row,false);
  
  // compute rQ
  rQ->copyFromArray(preCond->grad_x_val);
  MatQ->mult(1.0, *rQ,  1.0, *temp_x);
  MatJ->transMult( 1.0, *rQ, -1.0, *temp_y);

  if(!setIX){
  	setIX = true;
	temp_ixlow->setToZero();
	ixlow->copyIntoArrayWithIndex_AggVarCon(temp_ixlow->elements(),NULL,preCond->n_col,true);  
	temp_ixupp->setToZero();
	ixupp->copyIntoArrayWithIndex_AggVarCon(temp_ixupp->elements(),NULL,preCond->n_col,true); 	
	
	if( nxlow > 0 ){ 
	  temp_xl->setToZero();
	  prob->xlowerBound().copyIntoArrayWithIndex_AggVarCon(temp_xl->elements(),NULL,preCond->n_col,true);  
	}	
	if( nxupp > 0 ){
	  temp_xu->setToZero();
	  prob->xupperBound().copyIntoArrayWithIndex_AggVarCon(temp_xu->elements(),NULL,preCond->n_col,true);	
	}	
  }
	
  vars->gamma->selectNonZeros(*ixlow);
  vars->phi->selectNonZeros( *ixupp );

  if( nxlow > 0 ){ 
	temp_v->setToZero();
	vars->v->copyIntoArrayWithIndex_AggVarCon(temp_v->elements(),NULL,preCond->n_col,true);   	
	temp_gamma->setToZero();
	vars->gamma->copyIntoArrayWithIndex_AggVarCon(temp_gamma->elements(),NULL,preCond->n_col,true);  
	
    rQ->axpy( -1.0, *temp_gamma );
  }  
  if( nxupp > 0 ){ 
	temp_w->setToZero();
	vars->w->copyIntoArrayWithIndex_AggVarCon(temp_w->elements(),NULL,preCond->n_col,true); 
	temp_phi->setToZero();
	vars->phi->copyIntoArrayWithIndex_AggVarCon(temp_phi->elements(),NULL,preCond->n_col,true);	  
    	
    rQ->axpy( 1.0, *temp_phi );
  }

  // compute rA
  rA->copyFromArray(preCond->cons_b_val);
  MatJ->mult( -1.0, *rA, 1.0, *temp_x );


  if( nxlow > 0 ) {
	rv->copyFrom( *temp_x );
    rv->axpy( -1.0, *temp_xl ); 	
    rv->axpy( -1.0, *temp_v );

    rv->selectNonZeros( *temp_ixlow );

	rgamma->setToZero();
	rgamma->axzpy( 1.0, *temp_v, *temp_gamma );
  }
  
  if( nxupp > 0 ) {
	rw->copyFrom( *temp_x );
    rw->axpy( -1.0, *temp_xu );
    rw->axpy(  1.0, *temp_w );
		
    rw->selectNonZeros( *temp_ixupp );

	rphi->setToZero();
	rphi->axzpy( 1.0, *temp_w, *temp_phi );	
  }

}


void sLinsysRootAggregation::computeReducedRhs()
{

  NlpGenVars	   * vars  = (NlpGenVars *) (preCond->curr_Iter);
  NlpGenVars	   * step  = (NlpGenVars *) (preCond->curr_Step);

  assert( vars->validNonZeroPattern() );
  assert(locmz==0);

  temp_rhs_x->copyFrom( *rQ );
  if( nxlow > 0 ) {
  	// step now contains the rhs for the full sys. It has been computed in QpGenLinsys::solve
    OoqpVector & vInvGamma = *step->v;
  	temp_ncol->setToZero();
	vInvGamma.copyIntoArrayWithIndex_AggVarCon(temp_ncol->elements(),NULL,preCond->n_col,true);  
		
    temp_rhs_x->axzpy ( 1.0, *temp_ncol, *rv );
    temp_rhs_x->axdzpy( 1.0, *rgamma, *temp_v, *temp_ixlow );
  }

  if( nxupp > 0 ) {
    OoqpVector & wInvPhi	 = *step->w;
  	temp_ncol->setToZero();
	wInvPhi.copyIntoArrayWithIndex_AggVarCon(temp_ncol->elements(),NULL,preCond->n_col,true);  

    temp_rhs_x->axzpy (	1.0, *temp_ncol,	*rw );
    temp_rhs_x->axdzpy( -1.0, *rphi, *temp_w, *temp_ixupp );
  }

  temp_rhs_y->copyFrom( *rA );

}

void sLinsysRootAggregation::_set1stStVar(SimpleVector& redSol,SimpleVector& sol0)
{
  int startVarID = n_col - locnx;
  int startConID = n_row - locmy;
  
  memcpy( &sol0[0], &temp_rhs_x->elements()[startVarID], locnx * sizeof( double ) );
  if ( n_row > 0 ) memcpy( &sol0[locnx], &temp_rhs_y->elements()[startConID], locmy * sizeof( double ) );
}


void sLinsysRootAggregation::_joinRedRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			  OoqpVector& rhs2_in)
{
  SimpleVector & rhs  = (SimpleVector &) rhs_in;
  SimpleVector & rhs1 = (SimpleVector &) rhs1_in;
  SimpleVector & rhs2 = (SimpleVector &) rhs2_in;

  memcpy( &rhs[0], &rhs1[0], n_col * sizeof( double ) );
  if( n_row > 0 ) memcpy( &rhs[n_col],      &rhs2[0], n_row * sizeof( double ) );
}

void
sLinsysRootAggregation::_separateVars( OoqpVector& vars_in, OoqpVector& x_in, OoqpVector& y_in )
{
  SimpleVector & vars  = (SimpleVector &) vars_in;
  SimpleVector & x = (SimpleVector &) x_in;
  SimpleVector & y = (SimpleVector &) y_in;

  memcpy( &x[0], &vars[0], n_col * sizeof( double ) );
  if ( n_row > 0 ) memcpy( &y[0], &vars[n_col],      n_row * sizeof( double ) );
}





