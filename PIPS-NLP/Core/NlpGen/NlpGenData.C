/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */
 
 /* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "NlpGenData.h"
#include "NlpGenVars.h"
#include "NlpGenResiduals.h"

#include "DoubleMatrix.h"
#include "OoqpVector.h"
#include <cmath>

#include "SimpleVector.h"
#include "LinearAlgebraPackage.h"
#include "MpsReader.h"

#include "NlpInfo.h"

#include <stdlib.h>

using namespace std;
NlpGenData::NlpGenData()
{
 KryIter=0;
}

NlpGenData::NlpGenData(LinearAlgebraPackage * la_,
		     int nx_, int my_, int mz_,
		     int nnzQ_, int nnzA_, int nnzC_)
{
	la = la_;
	
	nx = nx_;
	my = my_;
	mz = mz_;
	
	H = SymMatrixHandle( la->newSymMatrix( nx,	   nnzQ_) );
	Jeq = GenMatrixHandle( la->newGenMatrix( my, nx, nnzA_ ) );
	Jineq= GenMatrixHandle( la->newGenMatrix( mz, nx, nnzC_ ) );
	
	grad  = OoqpVectorHandle( la->newVector( nx ) );
	blx   = OoqpVectorHandle( la->newVector( nx ) );
	ixlow = OoqpVectorHandle( la->newVector( nx ) );
	bux   = OoqpVectorHandle( la->newVector( nx ) );
	ixupp = OoqpVectorHandle( la->newVector( nx ) );

	CeqBody = OoqpVectorHandle( la->newVector( my ) );	
	bA	  = OoqpVectorHandle( la->newVector( my ) );

	CIneqBody = OoqpVectorHandle( la->newVector( mz ) );		
	bl	  = OoqpVectorHandle( la->newVector( mz ) );
	iclow = OoqpVectorHandle( la->newVector( mz ) );
	bu	  = OoqpVectorHandle( la->newVector( mz ) );
	icupp = OoqpVectorHandle( la->newVector( mz ) );
	sc	  = OoqpVectorHandle( la->newVector( nx  ) );

//	inputNlp = new NlpInfoAMPL(nx,my,mz,nnzQ_,nnzA_,nnzC_);
//	inputNlp = new NlpInfo(nx,my,mz,nnzQ_,nnzA_,nnzC_);


	trialBarrGrad_x  = OoqpVectorHandle( la->newVector( nx ) );	
	trialBarrGrad_s  = OoqpVectorHandle( la->newVector( mz ) );	


	trialCeqBody = OoqpVectorHandle( la->newVector( my ) );	
	trialCIneqBody = OoqpVectorHandle( la->newVector( mz ) );	

	dampind_xL_v  = OoqpVectorHandle( la->newVector( nx ) );
	dampind_xU_w  = OoqpVectorHandle( la->newVector( nx ) );	
	dampind_sL_t  = OoqpVectorHandle( la->newVector( mz ) );	
	dampind_sU_u  = OoqpVectorHandle( la->newVector( mz ) );		

	currMu = 0;
	linsysRes = 0;
	linsysRes_Full=0;
	KryIter=0;

	inputNlp = NULL;
	schurVarConID = NULL;
	schurSize = -1;
	setDampingVarMap();


	var_Part_idx_in = NULL;	
	con_Part_idx_in = NULL;	
}





NlpGenData::NlpGenData(LinearAlgebraPackage * la_,
		     long long  nx_, long long my_, long long mz_,
		     long long nnzQ_, long long nnzA_, long long nnzC_,
		     long long nxL_in,long long nxU_in,long long nsL_in,long long nsU_in)
{
	la = la_;
	
	nx = nx_;
	my = my_;
	mz = mz_;
	
	H = SymMatrixHandle( la->newSymMatrix( nx,	   nnzQ_) );
	Jeq = GenMatrixHandle( la->newGenMatrix( my, nx, nnzA_ ) );
	Jineq= GenMatrixHandle( la->newGenMatrix( mz, nx, nnzC_ ) );
	
	grad  = OoqpVectorHandle( la->newVector( nx ) );
	blx   = OoqpVectorHandle( la->newVector( nx ) );
	ixlow = OoqpVectorHandle( la->newVector( nx ) );
	bux   = OoqpVectorHandle( la->newVector( nx ) );
	ixupp = OoqpVectorHandle( la->newVector( nx ) );

	CeqBody = OoqpVectorHandle( la->newVector( my ) );	
	trialCeqBody = OoqpVectorHandle( la->newVector( my ) );	
	bA	  = OoqpVectorHandle( la->newVector( my ) );

	CIneqBody = OoqpVectorHandle( la->newVector( mz ) );	
	trialCIneqBody = OoqpVectorHandle( la->newVector( mz ) );	
	bl	  = OoqpVectorHandle( la->newVector( mz ) );
	iclow = OoqpVectorHandle( la->newVector( mz ) );
	bu	  = OoqpVectorHandle( la->newVector( mz ) );
	icupp = OoqpVectorHandle( la->newVector( mz ) );
	sc	  = OoqpVectorHandle( la->newVector( nx  ) );


	trialBarrGrad_x  = OoqpVectorHandle( la->newVector( nx ) );	
	trialBarrGrad_s  = OoqpVectorHandle( la->newVector( mz ) );	


	dampind_xL_v  = OoqpVectorHandle( la->newVector( nx ) );
	dampind_xU_w  = OoqpVectorHandle( la->newVector( nx ) );	
	dampind_sL_t  = OoqpVectorHandle( la->newVector( mz ) );	
	dampind_sU_u  = OoqpVectorHandle( la->newVector( mz ) );	


	currMu = 0;

    inputNlp = NULL;
	linsysRes = 0;
	linsysRes_Full=0;
	KryIter=0;

	
	nxlow = nxL_in; nxupp = nxU_in;
	mclow = nsL_in; mcupp = nsU_in;

	schurVarConID = NULL;
	schurSize = -1;
	
	setDampingVarMap();

	var_Part_idx_in = NULL;	
	con_Part_idx_in = NULL;		
}



NlpGenData::NlpGenData( LinearAlgebraPackage * la_in,
		      OoqpVector * grad_in, SymMatrix * H_in,
		      OoqpVector * xlow_in, OoqpVector * ixlow_in, long long nxlow_,
		      OoqpVector * xupp_in, OoqpVector * ixupp_in, long long nxupp_,
		      GenMatrix  * A_in, OoqpVector * bA_in,
		      GenMatrix  * C_in,
		      OoqpVector * clow_in, OoqpVector * iclow_in, long long mclow_,
		      OoqpVector * cupp_in, OoqpVector * icupp_in, long long mcupp_,
		      OoqpVector * CeqBody_in, OoqpVector * CIneqBody_in,
		      OoqpVector * trialBarrGrad_x_in,OoqpVector * trialBarrGrad_s_in,
	     	  OoqpVector * trialCeqBody_in, OoqpVector *trialCIneqBody_in, 
	     	  OoqpVector * dampind_xL_v_in,OoqpVector * dampind_xU_w_in,
	     	  OoqpVector * dampind_sL_t_in, OoqpVector *dampind_sU_u_in)
{

  SpReferTo( grad,     grad_in  );
  SpReferTo( bA,    bA_in );
  SpReferTo( blx,   xlow_in  );
  SpReferTo( ixlow, ixlow_in );
  SpReferTo( bux,   xupp_in  );
  SpReferTo( ixupp, ixupp_in );
  SpReferTo( bl,    clow_in  );
  SpReferTo( iclow, iclow_in );
  SpReferTo( bu,    cupp_in  );
  SpReferTo( icupp, icupp_in );
  SpReferTo( CeqBody,    CeqBody_in  );
  SpReferTo( CIneqBody, CIneqBody_in );


  long long dummy;
  la = la_in;

  nx = grad->length();
  SpReferTo( H, H_in );

  SpReferTo( Jeq, A_in );
  Jeq->getSize( my, dummy );
  
  SpReferTo( Jineq, C_in );
  Jineq->getSize( mz, dummy );
  
  currMu = 0;
  inputNlp = NULL;

  nxlow = ixlow->numberOfNonzeros(); 
  nxupp = ixupp->numberOfNonzeros(); 
  mclow = iclow->numberOfNonzeros(); 
  mcupp = icupp->numberOfNonzeros(); 


  SpReferTo( trialBarrGrad_x,    trialBarrGrad_x_in  );
  SpReferTo( trialBarrGrad_s, trialBarrGrad_s_in );
    
  SpReferTo( trialCeqBody,    trialCeqBody_in  );
  SpReferTo( trialCIneqBody, trialCIneqBody_in );


  SpReferTo( dampind_xL_v,    dampind_xL_v_in  );
  SpReferTo( dampind_xU_w, dampind_xU_w_in );
  SpReferTo( dampind_sL_t,    dampind_sL_t_in  );
  SpReferTo( dampind_sU_u, dampind_sU_u_in );

  linsysRes = 0;
  linsysRes_Full=0;
  KryIter=0;

  setDampingVarMap();
  schurVarConID = NULL;
  schurSize = -1;

  var_Part_idx_in = NULL;	
  con_Part_idx_in = NULL;	  
}

void
NlpGenData::setDampingVarMap()
{
    if( nxlow > 0 ) {
      //set value -1
      dampind_xL_v->setToConstant(-1);
      dampind_xL_v->selectNonZeros( *ixlow );

      if(nxupp>0){
	    //find index where var have both bounds, val = -1
	    dampind_xL_v->selectNonZeros( *ixupp);
	    //find index where var only have lower bound, now val=1 where only lb, otherwise 0
	    dampind_xL_v->addSomeConstants(1, *ixlow);
      }else{
	    dampind_xL_v->negate();
	  }
    }
  
    if( nxupp > 0 ) {
      //set value -1
      dampind_xU_w->setToConstant(-1);
      dampind_xU_w->selectNonZeros( *ixupp );

      if(nxlow>0){
	    //find index where var have both bounds, val = -1
	    dampind_xU_w->selectNonZeros( *ixlow);
	    //find index where var only have upper bound, now val=1 where only ub, otherwise 0
	    dampind_xU_w->addSomeConstants(1, *ixupp);
      }else{
	    dampind_xU_w->negate();
	  }
	}

    if( mclow > 0 ) {
      //set value -1
      dampind_sL_t->setToConstant(-1);
      dampind_sL_t->selectNonZeros( *iclow );

      if(mcupp>0){
	    //find index where var have both bounds, val = -1
	    dampind_sL_t->selectNonZeros( *icupp);
	    //find index where var only have lower bound, now val=1 where only lb, otherwise 0
	    dampind_sL_t->addSomeConstants(1, *iclow);
      }else{
	    dampind_sL_t->negate();
	  }
    }

    if( mcupp > 0 ) {
      //set value -1
      dampind_sU_u->setToConstant(-1);
      dampind_sU_u->selectNonZeros( *icupp );

      if(mclow>0){
	    //find index where var have both bounds, val = -1
	    dampind_sU_u->selectNonZeros( *iclow);
	    //find index where var only have upper bound, now val=1 where only ub, otherwise 0
	    dampind_sU_u->addSomeConstants(1, *icupp);
      }else{
	    dampind_sU_u->negate();
	  }
    }

}

NlpGenData::NlpGenData( LinearAlgebraPackage * la_in,
		      OoqpVector * grad_in, SymMatrix * H_in,
		      OoqpVector * xlow_in, OoqpVector * ixlow_in,
		      OoqpVector * xupp_in, OoqpVector * ixupp_in,
		      GenMatrix  * A_in, OoqpVector * bA_in,
		      GenMatrix  * C_in,
		      OoqpVector * clow_in, OoqpVector * iclow_in,
		      OoqpVector * cupp_in, OoqpVector * icupp_in)
{

  SpReferTo( grad,     grad_in  );
  SpReferTo( bA,    bA_in );
  SpReferTo( blx,   xlow_in  );
  SpReferTo( ixlow, ixlow_in );
  SpReferTo( bux,   xupp_in  );
  SpReferTo( ixupp, ixupp_in );
  SpReferTo( bl,    clow_in  );
  SpReferTo( iclow, iclow_in );
  SpReferTo( bu,    cupp_in  );
  SpReferTo( icupp, icupp_in );

  long long dummy;
  la = la_in;

  nx = grad->length();
  SpReferTo( H, H_in );

  SpReferTo( Jeq, A_in );
  Jeq->getSize( my, dummy );
  
  SpReferTo( Jineq, C_in );
  Jineq->getSize( mz, dummy );   

  currMu = 0;
  inputNlp = NULL;
}


void NlpGenData::Qmult( double beta,  OoqpVector& y,
		       double alpha, OoqpVector& x )
{
  H->mult( beta, y, alpha, x );
}

void NlpGenData::Amult( double beta,  OoqpVector& y,
		       double alpha, OoqpVector& x)
{
  Jeq->mult( beta, y, alpha, x );
}

void NlpGenData::Cmult( double beta,  OoqpVector& y,
		       double alpha, OoqpVector& x )
{
  Jineq->mult( beta, y, alpha, x );
}

void NlpGenData::ATransmult( double beta,  OoqpVector& y,
			    double alpha, OoqpVector& x )
{
  Jeq->transMult( beta, y, alpha, x );
}

void NlpGenData::CTransmult( double beta,  OoqpVector& y,
			    double alpha, OoqpVector& x )
{
  Jineq->transMult( beta, y, alpha, x );
}


void NlpGenData::getg( OoqpVector& myG )
{
  myG.copyFrom( *grad );
}

void NlpGenData::getbA( OoqpVector& bout )
{
  bout.copyFrom( *bA );
}

void NlpGenData::getInEqCons( OoqpVector& InEqCon_out )
{
  InEqCon_out.copyFrom( *CIneqBody);
}

double NlpGenData::datanorm()
{
  double norm = 0.0;
  double componentNorm;

  componentNorm = grad->infnorm();
  if( componentNorm > norm ) norm = componentNorm;
  
  componentNorm = H->abmaxnorm();
  if( componentNorm > norm ) norm = componentNorm;

  componentNorm = bA->infnorm();
  if( componentNorm > norm ) norm = componentNorm;

  componentNorm = Jeq->abmaxnorm();
  if( componentNorm > norm ) norm = componentNorm;

  componentNorm = Jineq->abmaxnorm();
  if( componentNorm > norm ) norm = componentNorm;

  assert( blx->matchesNonZeroPattern( *ixlow ) );
  componentNorm = blx->infnorm();
  if( componentNorm > norm ) norm = componentNorm;

  assert( bux->matchesNonZeroPattern( *ixupp ) );
  componentNorm = bux->infnorm();
  if( componentNorm > norm ) norm = componentNorm;

  assert( bl->matchesNonZeroPattern( *iclow ) );
  componentNorm = bl->infnorm();
  if( componentNorm > norm ) norm = componentNorm;

  assert( bu->matchesNonZeroPattern( *icupp ) );
  componentNorm = bu->infnorm();
  if( componentNorm > norm ) norm = componentNorm;

  return norm;
}

void NlpGenData::datainput( MpsReader * reader, int& iErr ) 
{
    reader->readQpGen( *grad, *H, *blx, *ixlow, *bux, *ixupp,
		     *Jeq, *bA,
		     *Jineq, *bl, *iclow, *bu, *icupp, iErr );

    if( reader->scalingOption == 1){
        // Create the scaling vector
        this->createScaleFromQ();

        //Scale the variables
        this->scaleQ();
        this->scaleA();
        this->scaleC();
        this->scaleg();
        this->scalexlow();
        this->scalexupp();
        }

    /* If objective sense is "MAX", flip the C and Q matrices */
    if( !strncmp( reader->objectiveSense, "MAX", 3)){
        this->flipg();
        this->flipQ();
        }  
}

void NlpGenData::print()
{
  cout << "begin Q\n";
  H->writeToStream( cout );
  cout << "end Q\n";
  cout << "begin c\n";
  grad->writeToStream( cout );
  cout << "end c\n";

  cout << "begin xlow\n";
  blx->writeToStream( cout );
  cout << "end xlow\n";
  cout << "begin ixlow\n";
  ixlow->writeToStream( cout );
  cout << "end ixlow\n";

  cout << "begin xupp\n";
  bux->writeToStream( cout );
  cout << "end xupp\n";  
  cout << "begin ixupp\n";
  ixupp->writeToStream( cout );
  cout << "end ixupp\n";
  cout << "begin A\n";

  Jeq->writeToStream( cout );
  cout << "end A\n";
  cout << "begin b\n";
  bA->writeToStream( cout );
  cout << "end b\n";
  cout << "begin C\n";
  Jineq->writeToStream( cout );
  cout << "end C\n";
  
  cout << "begin clow\n";
  bl->writeToStream( cout );
  cout << "end clow\n";
  cout << "begin iclow\n";
  iclow->writeToStream( cout );
  cout << "end iclow\n";

  cout << "begin cupp\n";
  bu->writeToStream( cout );
  cout << "end cupp\n";
  cout << "begin icupp\n";
  icupp->writeToStream( cout );
  cout << "end icupp\n";

}

#include <fstream>


void NlpGenData::putQIntoAt( GenMatrix& M, int row, int col )
{
  M.atPutSubmatrix( row, col, *H, 0, 0, nx, nx );
}

void NlpGenData::putQIntoAt( SymMatrix& M, int row, int col )
{
  M.symAtPutSubmatrix( row, col, *H, 0, 0, nx, nx );
}

void NlpGenData::putAIntoAt( GenMatrix& M, int row, int col )
{
  M.atPutSubmatrix( row, col, *Jeq, 0, 0, my, nx );
}

void NlpGenData::putAIntoAt( SymMatrix& M, int row, int col )
{
  M.symAtPutSubmatrix( row, col, *Jeq, 0, 0, my, nx );
}

void NlpGenData::putCIntoAt( GenMatrix& M, int row, int col )
{
  M.atPutSubmatrix( row, col, *Jineq, 0, 0, mz, nx );
}

void NlpGenData::putCIntoAt( SymMatrix& M, int row, int col )
{
  M.symAtPutSubmatrix( row, col, *Jineq, 0, 0, mz, nx );
}

void NlpGenData::setQIntoAt( SymMatrix& M, int row, int col, bool firstCall, std::map<int,int> &ValIdxMap )
{
  M.symAtSetSubmatrix( row, col, *H, 0, 0, nx, nx, firstCall, ValIdxMap);
}

void NlpGenData::setAIntoAt( SymMatrix& M, int row, int col, bool firstCall, std::map<int,int> &ValIdxMap )
{
  M.symAtSetSubmatrix( row, col, *Jeq, 0, 0, my, nx, firstCall, ValIdxMap);
}


void NlpGenData::setCIntoAt( SymMatrix& M, int row, int col, bool firstCall, std::map<int,int> &ValIdxMap )
{
  M.symAtSetSubmatrix( row, col, *Jineq, 0, 0, mz, nx, firstCall, ValIdxMap);
}






void NlpGenData::getDiagonalOfQ( OoqpVector& dq )
{
  H->fromGetDiagonal(0, dq);
}

/*
double NlpGenData::objectiveValue( NlpGenVars * vars )
{
  OoqpVectorHandle temp( la->newVector( nx ) );
  this->getg( *temp );
  this->Qmult( 1.0, *temp, 0.5, *vars->x );

  return temp->dotProductWith( *vars->x );
}
*/



void NlpGenData::createScaleFromQ()
{
  // Stuff the diagonal elements of Q into the vector "sc"
  this->getDiagonalOfQ( *sc);

  // Modifying scVector is equivalent to modifying sc
  SimpleVector & scVector = dynamic_cast<SimpleVector &>(*sc);

  int scLength = scVector.length();

  for( int i = 0; i < scLength; i++){
    if( scVector[i] > 1)
        scVector[i] = 1.0/sqrt( scVector[i]);
    else
        scVector[i] = 1.0;
    }
}

void NlpGenData::scaleQ()
{
    H->SymmetricScale( *sc);

}


void NlpGenData::scaleA()
{
    Jeq->ColumnScale( *sc);

}

void NlpGenData::scaleC()
{

    Jineq->ColumnScale( *sc);

}

void NlpGenData::scaleg()
{
  SimpleVector & scVector = dynamic_cast<SimpleVector &>(*sc);

  assert ( scVector.length() == grad->length());

  // D * g
  grad->componentMult( scVector);
}

void NlpGenData::scalexupp()
{
  SimpleVector & scVector = dynamic_cast<SimpleVector &>(*sc);

  assert ( scVector.length() == bux->length());

  // inverse(D) * bux
  bux->componentDiv( scVector);

}


void NlpGenData::scalexlow()
{
  SimpleVector & scVector = dynamic_cast<SimpleVector &>(*sc);

  assert ( scVector.length() == blx->length());

  // inverse(D) * blx
  blx->componentDiv( scVector);

}

void NlpGenData::flipg()
{
  // Multiply C matrix by -1
  grad->scalarMult( -1.0);
}

void NlpGenData::flipQ()
{
  // Multiply Q matrix by -1
  H->scalarMult( -1.0);
}


NlpGenData::~NlpGenData()
{
}


double 
NlpGenData::objectiveValue( Variables* vars_in )
{ 
  NlpGenVars *vars = (NlpGenVars *) vars_in;
  	
  if(vars)
    return inputNlp->ObjValue( vars );
  else 
  	return PriObj;
}


void 
NlpGenData::evalData(Variables* vars_in)
{
  NlpGenVars * vars = (NlpGenVars *) vars_in;
  int ifIncludeQx=0;
 
  PriObj  = objectiveValue(vars);
  BarrObj = BarrObjValue(vars,PriObj);

  ifIncludeQx = inputNlp->ObjGrad(vars,grad);
  if(!ifIncludeQx)
  	Qmult( 1.0, *grad,  1.0, *vars->x );

  inputNlp->ConstraintBody(vars, CeqBody, CIneqBody);
  inputNlp->JacFull(vars,Jeq,Jineq);
  inputNlp->Hessian(vars,H);
}




// evaluate barrier objective with damping factor
double
NlpGenData::BarrObjValue( NlpGenVars * vars, double PriObj_in, const double dampingFact)
{  
  double PriObjWrk,mu,BarrObjWrk=0;
  mu = currMu;

  if(PriObj_in!=0.0)
  	PriObjWrk = PriObj_in;
  else
  	PriObjWrk = objectiveValue( vars );
  

  // we can add   if( vars->t->numberOfNonzeros()>0)  before each call to skip calls;
  // Howevr, I think we would test the length of each vector when doing the evaluation (e.g addLog). 
  // hence we can skip the if routines here.) 
  if(vars->mclow>0) 
  	BarrObjWrk += vars->t->sumLog(iclow); 
  if(vars->mcupp>0) 
  	BarrObjWrk += vars->u->sumLog(icupp);
  if(vars->nxlow>0) 
  	BarrObjWrk += vars->v->sumLog(ixlow);
  if(vars->nxupp>0) 
  	BarrObjWrk += vars->w->sumLog(ixupp);

  BarrObjWrk *= -mu;
  BarrObjWrk += PriObjWrk;


  // Include the linear damping term if kappa_d is nonzero.
  if (dampingFact > 0) {
  	if(vars->mclow>0) 
	  BarrObjWrk += dampingFact * mu * vars->t->dotProductWith(*dampind_sL_t); 
    if(vars->mcupp>0) 
	  BarrObjWrk += dampingFact * mu * vars->u->dotProductWith(*dampind_sU_u);
	if(vars->nxlow>0) 
	  BarrObjWrk += dampingFact * mu * vars->v->dotProductWith(*dampind_xL_v);
	if(vars->nxupp>0)
	  BarrObjWrk += dampingFact * mu * vars->w->dotProductWith(*dampind_xU_w);
  }

  return(BarrObjWrk);

}




// evaluate merit function  with 2 norm
// FIXME_NY: with 2 norm only. I think this routine should be moved to the upper layer (alg level?)
// 
double 
NlpGenData::evalMeritFunc( NlpGenVars * vars, double penalty, double BarrObj_in)
{  
  double MeritFuncWrk=0,sum=0;


  if(BarrObj_in != 0.0)
  	MeritFuncWrk = BarrObj_in;
  else
  	MeritFuncWrk = BarrObjValue( vars, 0.0 );

  sum += vars->y->sumPowElt(2);
  sum += vars->z->sumPowElt(2);
  
  sum += vars->lambda->sumPowElt(2);
  sum += vars->pi->sumPowElt(2);
  
  sum += vars->gamma->sumPowElt(2);
  sum += vars->phi->sumPowElt(2);

  sum = sqrt(sum);

  MeritFuncWrk += -penalty*sum;

  
  return(MeritFuncWrk);

}


// evaluate merit function  with 2 norm
// 
double 
NlpGenData::evalMeritFunc( double priErr, double penalty, double BarrObj_in)
{  
  double MeritFuncWrk = 0;

  MeritFuncWrk = BarrObj_in;
  MeritFuncWrk += penalty*priErr;
  
  return(MeritFuncWrk);

}

// evaluate merit function  
double 
NlpGenData::evalMeritFunc( double BarrObj_in, Variables *iterate_in, Residuals *resid_in)
{  
  NlpGenVars *vars = dynamic_cast<NlpGenVars *> (iterate_in);
  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);  

  double MeritFuncWrk = BarrObj_in;

  // kkt.metric = BarrObj + lam'*kkt.lam; 
  MeritFuncWrk += (vars->y)->dotProductWith(*resid->rA);
  MeritFuncWrk += (vars->z)->dotProductWith(*resid->rC);

  // FIXME_NY: do we need to include the additional constraints due to the slacks?
/*
  if( mclow > 0)
  	MeritFuncWrk+= (vars->lambda)->dotProductWith(*resid->rt);
  if( mcupp > 0)
  	MeritFuncWrk += (vars->pi)->dotProductWith(*resid->ru);
  if( nxlow > 0)
  	MeritFuncWrk += (vars->gamma)->dotProductWith(*resid->rv);
  if( nxupp > 0)
  	MeritFuncWrk += (vars->phi)->dotProductWith(*resid->rw);  
*/
  
  return(MeritFuncWrk);
}

// evaluate merit function  
double 
NlpGenData::evalScaledConstraintNorm( Variables *iterate_in, Residuals *resid_in, const int isTrialStep )
{  
  NlpGenVars *vars = dynamic_cast<NlpGenVars *> (iterate_in);
  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);  

  double kktThWrk;

  // kkt.th = 0.5*norm(kkt.lam)^2;
  kktThWrk = 0.5 * resid->getKKTRhsNorm_Primal(this,vars,PIPS_NORM_TWO, isTrialStep);

  // FIXME_NY: do we need to include the additional constraints due to the slacks?
  // kktThWrk = 0.5 * resid->priErr();
  
  return(kktThWrk);
}




// evaluate constraint body
void 
NlpGenData::evalConstraintBody(Variables* vars_in, const int IfTrialStep)
{
  NlpGenVars * vars = (NlpGenVars *) vars_in;
  
  if(IfTrialStep==0){
    inputNlp->ConstraintBody(vars, CeqBody, CIneqBody);
  }else{
    inputNlp->ConstraintBody(vars, trialCeqBody, trialCIneqBody);
  }
}



// evaluate BarrGrad times vector(direction), use trialBarrGrad_x as temp vec
double 
NlpGenData::getBarrGradTimesD(Variables* vars_in, Variables* step_in, 
					const int IfTrialStep, const double dampingFact)
{
  double result = 0.0;

  NlpGenVars * vars = (NlpGenVars *) vars_in;
  NlpGenVars * sstep = (NlpGenVars *) step_in;

  trialBarrGrad_x->copyFrom(*grad);
  result += trialBarrGrad_x->dotProductWith( *sstep->x );

  if( nxlow > 0 ) {
    trialBarrGrad_x->setToConstant(-currMu);
	trialBarrGrad_x->divideSome( *vars->v, *ixlow );
	trialBarrGrad_x->selectNonZeros( *ixlow );
	
	result += trialBarrGrad_x->dotProductWith( *sstep->x );
  }
  if( nxupp > 0 ) {
    trialBarrGrad_x->setToConstant(currMu);
	trialBarrGrad_x->divideSome( *vars->w, *ixupp );
	trialBarrGrad_x->selectNonZeros( *ixupp );
	
	result += trialBarrGrad_x->dotProductWith( *sstep->x );
  }

  if( mclow > 0 ) {
    trialBarrGrad_s->setToConstant(-currMu);
	trialBarrGrad_s->divideSome( *vars->t, *iclow );
	trialBarrGrad_s->selectNonZeros( *iclow );
	
	result += trialBarrGrad_s->dotProductWith( *sstep->s );
  }
  if( mcupp > 0 ) {
    trialBarrGrad_s->setToConstant(currMu);
	trialBarrGrad_s->divideSome( *vars->u, *icupp );
	trialBarrGrad_s->selectNonZeros( *icupp );
	
	result += trialBarrGrad_s->dotProductWith( *sstep->s );
  }

  // Include the linear damping term if kappa_d is nonzero.
  if (dampingFact > 0.0) {

	double dampingTimeMu = dampingFact*currMu;

	if(nxlow>0) 
	  result += dampingTimeMu * dampind_xL_v->dotProductWith(*sstep->x);
	if(nxupp>0)
	  result -= dampingTimeMu  * dampind_xU_w->dotProductWith(*sstep->x);	
	
  	if(mclow>0) 
	  result += dampingTimeMu * dampind_sL_t->dotProductWith(*sstep->s); 
    if(mcupp>0) 
	  result -= dampingTimeMu  * dampind_sU_u->dotProductWith(*sstep->s);
  }

  return result;
}



void
NlpGenData::moveBounds(OoqpVector *priWrk_X, OoqpVector *priWrk_S, const double move_bound_scalar)
{
  if( mclow > 0 ) {
    // Calculate priWrk_S = max(|lb|,1)
    priWrk_S->absVal(bl);
	priWrk_S->SetComponentFromMaxXorConstant(priWrk_S,1.0);
	
	bl->axpy(-move_bound_scalar,*priWrk_S);
	bl->selectNonZeros(*iclow);
  }
  if( mcupp > 0 ) {
    // Calculate priWrk_S = max(|ub|,1)
    priWrk_S->absVal(bu);
	priWrk_S->SetComponentFromMaxXorConstant(priWrk_S,1.0);
	
	bu->axpy(move_bound_scalar,*priWrk_S);
	bu->selectNonZeros(*icupp);
  }
  if( nxlow > 0 ) {
    // Calculate priWrk_X = max(|lb|,1)
    priWrk_X->absVal(blx);
	priWrk_X->SetComponentFromMaxXorConstant(priWrk_X,1.0);
	
	blx->axpy(-move_bound_scalar,*priWrk_X);
	blx->selectNonZeros(*ixlow);
  }
  if( nxupp > 0 ) {
    // Calculate priWrk_X = max(|ub|,1)
    priWrk_X->absVal(bux);
	priWrk_X->SetComponentFromMaxXorConstant(priWrk_X,1.0);
	
	bux->axpy(move_bound_scalar,*priWrk_X);
	bux->selectNonZeros(*ixupp);
  }

}


void NlpGenData::getInitX(OoqpVector *initVecX)
{
  inputNlp->get_InitX0(initVecX);
}



// evaluate cons times vector(direction), use trialBarrGrad_y as temp vec
double 
NlpGenData::getConTimesD(Variables* vars_in, Variables* step_in, Residuals *resid_in)
{
  NlpGenVars *steps = dynamic_cast<NlpGenVars *> (step_in);
  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);  

  double ConsResidTd = 0.0;

  //  ConsResidTd = d'*kkt.lam; 
  ConsResidTd += (steps->y)->dotProductWith(*resid->rA);
  ConsResidTd += (steps->z)->dotProductWith(*resid->rC);
  
  return ConsResidTd;
}


