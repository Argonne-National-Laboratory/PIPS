/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QpGenLinsys.h"

#include "QpGenData.h"
#include "QpGenResiduals.h"
#include "QpGenVars.h"

#include "OoqpVector.h"
#include "DoubleMatrix.h"
#include "DoubleLinearSolver.h"
#include "SimpleVector.h"
#include "LinearAlgebraPackage.h"
#include "QpGen.h"

#ifdef TIMING
#include "mpi.h"
#include <vector>
#endif


#include <fstream>
using namespace std;

extern int gOuterIterRefin;

QpGenLinsys::QpGenLinsys( QpGen * factory_,
			  QpGenData * prob,
			  LinearAlgebraPackage * la ) 
: factory( factory_), rhs(NULL), dd(NULL), dq(NULL), useRefs(0)
{

  nx = prob->nx; my = prob->my; mz = prob->mz;
  ixlow = prob->ixlow;
  ixupp = prob->ixupp;
  iclow = prob->iclow;
  icupp = prob->icupp;

  nxlow = ixlow->numberOfNonzeros();
  nxupp = ixupp->numberOfNonzeros();
  mclow = iclow->numberOfNonzeros();
  mcupp = icupp->numberOfNonzeros();

  if( nxupp + nxlow > 0 ) {
    dd      = la->newVector( nx );
    dq      = la->newVector( nx );
    prob->getDiagonalOfQ( *dq );
  }
  nomegaInv   = la->newVector( mz );
  rhs         = la->newVector( nx + my + mz );

  if(gOuterIterRefin) {
    //for iterative refinement
    sol  = la->newVector( nx + my + mz );
    res  = la->newVector( nx + my + mz );
    resx = la->newVector( nx );
    resy = la->newVector( my );
    resz = la->newVector( mz );
  } else {
    sol = res = resx = resy = resz = NULL;
  }
}

QpGenLinsys::QpGenLinsys()
 : factory( NULL), rhs(NULL), dd(NULL), dq(NULL), useRefs(0),
   sol(NULL), res(NULL), resx(NULL), resy(NULL), resz(NULL)
{
}

QpGenLinsys::~QpGenLinsys()
{
  if(!useRefs) {
    delete dd; delete dq;
    delete rhs;
    delete nomegaInv;
  }

  if(sol)  delete sol;
  if(res)  delete res;
  if(resx) delete resx;
  if(resy) delete resy;
  if(resz) delete resz;
  
}
void QpGenLinsys::factor(Data * /* prob_in */, Variables *vars_in)
{
  QpGenVars * vars = (QpGenVars *) vars_in;

  assert( vars->validNonZeroPattern() );

  if( nxlow + nxupp > 0 ) dd->copyFrom(*dq);
  this->computeDiagonals( *dd, *nomegaInv,
			  *vars->t, *vars->lambda,
			  *vars->u, *vars->pi,
			  *vars->v, *vars->gamma,
			  *vars->w, *vars->phi );
  if( nxlow + nxupp > 0 ) this->putXDiagonal( *dd );

  nomegaInv->invert();
  nomegaInv->negate();

  if( mclow + mcupp > 0 ) this->putZDiagonal( *nomegaInv );
 
}


void QpGenLinsys::computeDiagonals( OoqpVector& dd_, OoqpVector& omega,
				    OoqpVector& t,  OoqpVector& lambda,
				    OoqpVector& u,  OoqpVector& pi,
				    OoqpVector& v,  OoqpVector& gamma,
				    OoqpVector& w,  OoqpVector& phi )
{
  if( nxupp + nxlow > 0 ) {
    if( nxlow > 0 ) dd_.axdzpy( 1.0, gamma, v, *ixlow );
    if( nxupp > 0 ) dd_.axdzpy( 1.0, phi  , w, *ixupp );
  }
  omega.setToZero();
  if ( mclow > 0 ) omega.axdzpy( 1.0, lambda, t, *iclow );
  if ( mcupp > 0 ) omega.axdzpy( 1.0, pi,     u, *icupp );
  //assert( omega->allPositive() );
}

void QpGenLinsys::solve(Data * prob_in, Variables *vars_in,
			Residuals *res_in, Variables *step_in)
{
  QpGenData      * prob  = (QpGenData *) prob_in;
  QpGenVars      * vars  = (QpGenVars *) vars_in;
  QpGenVars      * step  = (QpGenVars *) step_in;
  QpGenResiduals * res   = (QpGenResiduals *) res_in;


  assert( vars->validNonZeroPattern() );
  assert( res ->validNonZeroPattern() );
  
  step->x->copyFrom( *res->rQ );
  if( nxlow > 0 ) {
    OoqpVector & vInvGamma = *step->v;
    vInvGamma.copyFrom( *vars->gamma );
    vInvGamma.divideSome( *vars->v, *ixlow );
	
    step->x->axzpy ( 1.0, vInvGamma, *res->rv );
    step->x->axdzpy( 1.0, *res->rgamma, *vars->v, *ixlow );
  }
  if( nxupp > 0 ) {
    OoqpVector & wInvPhi   = *step->w;
    wInvPhi.copyFrom( *vars->phi );
    wInvPhi.divideSome( *vars->w, *ixupp );
	  
    step->x->axzpy (  1.0, wInvPhi,   *res->rw );
    step->x->axdzpy( -1.0, *res->rphi, *vars->w, *ixupp );
  }
  // start by partially computing step->s
  step->s->copyFrom( *res->rz );
  if( mclow > 0 ) {
    OoqpVector & tInvLambda = *step->t;
	
    tInvLambda.copyFrom( *vars->lambda );
    tInvLambda.divideSome( *vars->t, *iclow );

    step->s-> axzpy( 1.0, tInvLambda, *res->rt );
    step->s->axdzpy( 1.0, *res->rlambda, *vars->t, *iclow );
  }

  if( mcupp > 0 ) {
    OoqpVector & uInvPi = *step->u;
	
    uInvPi.copyFrom( *vars->pi );
    uInvPi.divideSome( *vars->u, *icupp );

    step->s-> axzpy(  1.0, uInvPi, *res->ru );
    step->s->axdzpy( -1.0, *res->rpi, *vars->u, *icupp );
  }

  step->y->copyFrom( *res->rA );
  step->z->copyFrom( *res->rC );

  {
    // Unfortunately, we need a temporary  OoqpVector for the solve,
    // Use step->lambda or step->pi
    OoqpVectorHandle ztemp;
    if( mclow > 0 ) {
      ztemp = step->lambda;
    } else {
      ztemp = step->pi;
    }

    this->solveXYZS( *step->x, *step->y, *step->z, *step->s,
		     *ztemp, prob );
  }

  if( mclow > 0 ) {
    step->t->copyFrom( *step->s );
    step->t->axpy( -1.0, *res->rt );
    step->t->selectNonZeros( *iclow );

    step->lambda->copyFrom( *res->rlambda );
    step->lambda->axzpy( -1.0, *vars->lambda, *step->t );
    step->lambda->divideSome( *vars->t, *iclow );
  }
  if( mcupp > 0 ) {
    step->u->copyFrom( *res->ru );
    step->u->axpy( -1.0, *step->s );
    step->u->selectNonZeros( *icupp );

    step->pi->copyFrom( *res->rpi );
    step->pi->axzpy( -1.0, *vars->pi, *step->u );
    step->pi->divideSome( *vars->u, *icupp );
  }
  if( nxlow > 0 ) {
    step->v->copyFrom( *step->x );
    step->v->axpy( -1.0, *res->rv );
    step->v->selectNonZeros( *ixlow );
	
    step->gamma->copyFrom( *res->rgamma );
    step->gamma->axzpy( -1.0, *vars->gamma, *step->v );
    step->gamma->divideSome( *vars->v, *ixlow );
  }
  if( nxupp > 0 ) {
    step->w->copyFrom( *res->rw );
    step->w->axpy( -1.0, *step->x );
    step->w->selectNonZeros( *ixupp );
	
    step->phi->copyFrom( *res->rphi );
    step->phi->axzpy( -1.0, *vars->phi, *step->w );
    step->phi->divideSome( *vars->w, *ixupp );
  }
  assert( step->validNonZeroPattern() );

}


void QpGenLinsys::solveXYZS( OoqpVector& stepx, OoqpVector& stepy,
			       OoqpVector& stepz, OoqpVector& steps,
			       OoqpVector& /* ztemp */,
			       QpGenData* prob )
{

#ifdef TIMING
    int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    vector<double> histRelResid;
    //vector<double> histRelResidInf;
    if(0==myRank) cout << "solveXYZS solve ---" << endl;
#endif
  stepz.axzpy( -1.0, *nomegaInv, steps );
  this->joinRHS( *rhs, stepx, stepy, stepz );

  if(gOuterIterRefin) {
    
#ifdef TIMING
      double tTot=MPI_Wtime(), tSlv=0., tResid=0., tTmp;
#endif
    res->copyFrom(*rhs);
    sol->setToZero();
    double bnorm=rhs->twonorm();
    int refinSteps=-1;  double resNorm;

    do{
#ifdef TIMING
      tTmp=MPI_Wtime();
#endif
      this->solveCompressed( *res );
#ifdef TIMING
      tSlv += (MPI_Wtime()-tTmp);
#endif

      //x = x+dx
      sol->axpy(1.0, *res);
      refinSteps++;
      

#ifdef TIMING
      tTmp=MPI_Wtime();
#endif
      res->copyFrom(*rhs);    
      //  stepx, stepy, stepz are used as temporary buffers
      computeResidualXYZ( *sol, *res, stepx, stepy, stepz, prob );
#ifdef TIMING
      tResid += (MPI_Wtime()-tTmp);
#endif 

      resNorm=res->twonorm(); 
#ifdef TIMING
      histRelResid.push_back(resNorm/bnorm);
      //histRelResidInf.push_back(res->infnorm()/bnorm);
      if(0==myRank) cout << "resid.nrm xyz: " << resNorm << "   "
			 << "rhs.nrm xyz: " << bnorm << endl;
#endif      
      
      if(resNorm/bnorm<1e-12)
	break;

    } while(true);
#ifdef TIMING
    tTot = MPI_Wtime() - tTot;
    if(0==myRank) {// && refinSteps>0)  {
      cout << "Outer " << refinSteps << " iterative steps. Rel resids: "; //Norm rel res:" 
      for(size_t it=0; it<histRelResid.size(); it++) 
	  cout << histRelResid[it] << " | ";
      cout << endl;
      cout << "Outer times:  solve=" << tSlv 
	   << "  resid=" << tResid 
	   << "  total=" << tTot << endl; 
    }
#endif
    this->separateVars( stepx, stepy, stepz, *sol );
  } else {
    this->solveCompressed( *rhs );
    this->separateVars( stepx, stepy, stepz, *rhs );
  }
  stepy.negate();
  stepz.negate();
	
  steps.axpy( -1.0, stepz );
  steps.componentMult( *nomegaInv );
  steps.negate();
}

/**
 * res = res - mat*sol
 * stepx, stepy, stepz are used as temporary buffers
 */
void QpGenLinsys::computeResidualXYZ(OoqpVector& sol, 
				     OoqpVector& res, 
				     OoqpVector& solx, 
				     OoqpVector& soly, 
				     OoqpVector& solz, 
				     QpGenData* data)
{
  this->separateVars( solx, soly, solz, sol );
  this->separateVars( *resx, *resy, *resz, res);

  data->Qmult(1.0, *resx, -1.0, solx);
  resx->axzpy(-1.0, *dd, solx);
  data->ATransmult(1.0, *resx, -1.0, soly);
  data->CTransmult(1.0, *resx, -1.0, solz);
  //cout << "resx norm: " << resx->twonorm() << endl;
  
  data->Amult(1.0, *resy, -1.0, solx);
  //cout << "resy norm: " << resy->twonorm() << endl;
  data->Cmult(1.0, *resz, -1.0, solx);
  resz->axzpy(-1.0, *nomegaInv, solz);
  //cout << "resz norm: " << resz->twonorm() << endl;
  this->joinRHS( res, *resx, *resy, *resz );
}

void QpGenLinsys::joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			     OoqpVector& rhs2_in, OoqpVector& rhs3_in )
{
  // joinRHS has to be delegated to the factory. This is true because
  // the rhs may be distributed across processors, so the factory is the
  // only object that knows with certainly how to scatter the elements.
  factory->joinRHS( rhs_in, rhs1_in, rhs2_in, rhs3_in );
}

void QpGenLinsys::separateVars( OoqpVector& x_in, OoqpVector& y_in,
				  OoqpVector& z_in, OoqpVector& vars_in )
{
  // separateVars has to be delegated to the factory. This is true because
  // the rhs may be distributed across processors, so the factory is the
  // only object that knows with certainly how to scatter the elements.
  factory->separateVars( x_in, y_in, z_in, vars_in );
}

