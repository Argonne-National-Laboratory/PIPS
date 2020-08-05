/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QpGenResiduals.h"
#include "QpGenVars.h"
#include "QpGenData.h"

#include "OoqpVector.h"
#include "LinearAlgebraPackage.h"

#include "pipsdef.h"

#include <iostream>
#include <fstream>

#include "mpi.h"

QpGenResiduals::QpGenResiduals( LinearAlgebraPackage * la,
				long long nx_, long long my_, long long mz_,
				OoqpVector * ixlow_in, OoqpVector * ixupp_in,
				OoqpVector * iclow_in, OoqpVector * icupp_in )
{
  /* constructor not in use in PIPS-IPM */
  assert(false); 
  nx = nx_;
  my = my_;
  mz = mz_;

  SpReferTo( ixlow, ixlow_in );
  nxlow = ixlow->numberOfNonzeros();

  SpReferTo( ixupp, ixupp_in );
  nxupp = ixupp->numberOfNonzeros();

  SpReferTo( iclow, iclow_in );
  mclow = iclow->numberOfNonzeros();

  SpReferTo( icupp, icupp_in );
  mcupp = icupp->numberOfNonzeros();

  rQ = OoqpVectorHandle( la->newVector( nx ) );
  rA = OoqpVectorHandle( la->newVector( my ) );
  rC = OoqpVectorHandle( la->newVector( mz ) );

  rz = OoqpVectorHandle( la->newVector( mz ) );
  if ( mclow > 0 ) {
    rt      = OoqpVectorHandle( la->newVector( mz ) );
    rlambda = OoqpVectorHandle( la->newVector( mz ) );
  }
  if ( mcupp > 0 ) {
    ru     = OoqpVectorHandle( la->newVector( mz ) );
    rpi    = OoqpVectorHandle( la->newVector( mz ) );
  }
  if( nxlow > 0 ) {
    rv     = OoqpVectorHandle( la->newVector( nx ) );
    rgamma = OoqpVectorHandle( la->newVector( nx ) );
  }
  if( nxupp > 0 ) {
    rw   = OoqpVectorHandle( la->newVector( nx ) );
    rphi = OoqpVectorHandle( la->newVector( nx ) );
  }
}

QpGenResiduals::QpGenResiduals( const QpGenResiduals& res) : Residuals(res)
{
  nx = res.nx;
  my = res.my;
  mz = res.mz;

  ixlow = OoqpVectorHandle( res.ixlow->cloneFull() );
  nxlow = res.nxlow;

  ixupp = OoqpVectorHandle( res.ixupp->cloneFull() );
  nxupp = res.nxupp;

  iclow = OoqpVectorHandle( res.iclow->cloneFull() );
  mclow = res.mclow;

  icupp = OoqpVectorHandle( res.icupp->cloneFull() );
  mcupp = res.mcupp;

  rQ = OoqpVectorHandle(res.rQ->cloneFull());
  rA = OoqpVectorHandle(res.rA->cloneFull());
  rC = OoqpVectorHandle(res.rC->cloneFull());

  rz = OoqpVectorHandle(res.rz->cloneFull());

  rt = OoqpVectorHandle(res.rt->cloneFull());
  rlambda = OoqpVectorHandle(res.rlambda->cloneFull());

  ru = OoqpVectorHandle(res.ru->cloneFull());
  rpi = OoqpVectorHandle(res.rpi->cloneFull());

  rv = OoqpVectorHandle(res.rv->cloneFull());
  rgamma = OoqpVectorHandle(res.rgamma->cloneFull());
  
  rw = OoqpVectorHandle(res.rw->cloneFull());
  rphi = OoqpVectorHandle(res.rphi->cloneFull());
}


void QpGenResiduals::calcresids(Data *prob_in, Variables *vars_in, bool print_resids)
{
#ifdef TIMING
  print_resids = true;
#endif
  const int myRank = PIPS_MPIgetRank(MPI_COMM_WORLD);
  QpGenVars * vars = (QpGenVars *) vars_in;
  QpGenData * prob = (QpGenData *) prob_in;

  double componentNorm, norm = 0.0, gap = 0.0;
 
  prob->getg( *rQ );
  prob->Qmult( 1.0, *rQ,  1.0, *vars->x );

  // calculate x^T (g+Qx) - contribution to the duality gap
  gap += rQ->dotProductWith(*vars->x); 

  prob->ATransmult( 1.0, *rQ, -1.0, *vars->y );
  prob->CTransmult( 1.0, *rQ, -1.0, *vars->z );

  vars->gamma->selectNonZeros(*ixlow);
  vars->phi->selectNonZeros( *ixupp );
  if( nxlow > 0 ) rQ->axpy( -1.0, *vars->gamma );
  if( nxupp > 0 ) rQ->axpy(  1.0, *vars->phi );

  componentNorm = rQ->infnorm();

  if( print_resids )
  {
     const double rQtwonorm=rQ->twonorm();
     if( 0 == myRank )
     {
//        std::cout << gap << std::endl;
        std::cout << " rQ infnorm = " << componentNorm << " | twonorm = " << rQtwonorm << std::endl;
     }
  }

  if( componentNorm > norm ) norm = componentNorm;

  prob->getbA( *rA );
  prob->Amult( -1.0, *rA, 1.0, *vars->x );

  //printf("gap2=%20.16f\n", gap);
  
  // contribution -d^T y to duality gap
  gap -= prob->bA->dotProductWith(*vars->y);

  componentNorm = rA->infnorm();
  if( print_resids )
  {
     if( 0 == myRank )
     {
//        std::cout << gap << std::endl;
        std::cout << " rA norm = " << componentNorm << std::endl;
     }
  }
  if( componentNorm > norm ) norm = componentNorm;

  rC->copyFrom( *vars->s );
  prob->Cmult( -1.0, *rC, 1.0, *vars->x );

  componentNorm = rC->infnorm();
  if( print_resids )
  {
     if( 0 == myRank )
        std::cout << " rC norm = " << componentNorm << std::endl;
  }
  if( componentNorm > norm ) norm = componentNorm;

  rz->copyFrom( *vars->z );

  if( mclow > 0 ) {
    rz->axpy( -1.0, *vars->lambda );

#ifdef TIMING
    componentNorm =  (vars->lambda)->infnorm();
    if( 0 == myRank )
       std::cout << "lambda norm " << componentNorm << std::endl;
#endif


    rt->copyFrom( *vars->s );
    rt->axpy( -1.0, prob->slowerBound() );
    rt->selectNonZeros( *iclow );
    rt->axpy( -1.0, *vars->t );
    gap -= prob->bl->dotProductWith(*vars->lambda);

    componentNorm = rt->infnorm();
    if( print_resids )
    {
       if( 0 == myRank )
       {
//          std::cout << gap << std::endl;
          std::cout << " rt norm = " << componentNorm << std::endl;
       }
    }
    if( componentNorm > norm ) norm = componentNorm;
  }
  if( mcupp > 0 )
  {
    rz->axpy(  1.0, *vars->pi );
#ifdef TIMING
    componentNorm = (vars->pi)->infnorm();
    if( 0 == myRank )
       std::cout << "pi norm " << componentNorm << std::endl;
#endif

    ru->copyFrom( *vars->s );
    ru->axpy( -1.0, prob->supperBound() );
    ru->selectNonZeros( *icupp );
    ru->axpy( 1.0, *vars->u );

    gap += prob->bu->dotProductWith(*vars->pi);

    componentNorm = ru->infnorm();
    if( print_resids )
    {
       if( 0 == myRank )
       {
//          std::cout << gap << std::endl;
          std::cout << " ru norm = " << componentNorm << std::endl;
       }
    }
    if( componentNorm > norm ) norm = componentNorm;
  }
  componentNorm = rz->infnorm();
  if( print_resids )
  {
     if( 0 == myRank )
        std::cout << " rz norm = " << componentNorm << std::endl;

  }
  if( componentNorm > norm ) norm = componentNorm;

  if( nxlow > 0 )
  {
    rv->copyFrom( *vars->x );
    rv->axpy( -1.0, prob->xlowerBound() );
    rv->selectNonZeros( *ixlow );
    rv->axpy( -1.0, *vars->v );

    gap -= prob->blx->dotProductWith(*vars->gamma);

    componentNorm = rv->infnorm();
    if( print_resids )
    {
       if( 0 == myRank )
       {
//          std::cout << gap << std::endl;
          std::cout << " rv norm = " << componentNorm << std::endl;
       }
    }
    if( componentNorm > norm ) norm = componentNorm;
  }
  if( nxupp > 0 )
  {
    rw->copyFrom( *vars->x );
    rw->axpy( -1.0, prob->xupperBound() );
    rw->selectNonZeros( *ixupp );
    rw->axpy(  1.0, *vars->w );

    gap += prob->bux->dotProductWith(*vars->phi);

    componentNorm = rw->infnorm();
    if( print_resids )
    {
       if( 0 == myRank )
       {
//          std::cout << gap << std::endl;
          std::cout << " rw norm = " << componentNorm << std::endl;
       }
    }
    if( componentNorm > norm ) norm = componentNorm;
  }
   
  mDualityGap = gap;
  mResidualNorm = norm;

  if( print_resids && myRank == 0)
  {
     std::cout << "Norm residuals: " << mResidualNorm << "\tduality gap: " << mDualityGap << std::endl;
  }
}
  

void QpGenResiduals::add_r3_xz_alpha(Variables *vars_in, double alpha)
{
  QpGenVars * vars = (QpGenVars *) vars_in;

  if( mclow > 0 ) rlambda->axzpy( 1.0, *vars->t, *vars->lambda );
  if( mcupp > 0 ) rpi    ->axzpy( 1.0, *vars->u, *vars->pi );
  if( nxlow > 0 ) rgamma ->axzpy( 1.0, *vars->v, *vars->gamma );
  if( nxupp > 0 ) rphi   ->axzpy( 1.0, *vars->w, *vars->phi );

  if( alpha != 0.0 ) {
    if( mclow > 0 ) rlambda->addSomeConstants( alpha, *iclow );
    if( mcupp > 0 ) rpi    ->addSomeConstants( alpha, *icupp );
    if( nxlow > 0 ) rgamma ->addSomeConstants( alpha, *ixlow );
    if( nxupp > 0 ) rphi   ->addSomeConstants( alpha, *ixupp );
  }
}

void QpGenResiduals::set_r3_xz_alpha(Variables *vars, double alpha)
{
  this->clear_r3();
  this->add_r3_xz_alpha( vars, alpha );
}
  
void QpGenResiduals::clear_r3()
{
  if( mclow > 0 ) rlambda->setToZero();
  if( mcupp > 0 ) rpi    ->setToZero();
  if( nxlow > 0 ) rgamma ->setToZero();
  if( nxupp > 0 ) rphi   ->setToZero();
}
  
void QpGenResiduals::clear_r1r2()
{
  rQ->setToZero();
  rA->setToZero();
  rC->setToZero();
  rz->setToZero();
  if( nxlow > 0 ) rv->setToZero();
  if( nxupp > 0 ) rw->setToZero();
  if( mclow > 0 ) rt->setToZero();
  if( mcupp > 0 ) ru->setToZero();
}

void QpGenResiduals::project_r3(double rmin, double rmax)
{
  if( mclow > 0 ) {
    rlambda->gondzioProjection( rmin, rmax );
    rlambda->selectNonZeros( *iclow );
  }
  if( mcupp > 0 ) {
    rpi    ->gondzioProjection( rmin, rmax );
    rpi    ->selectNonZeros( *icupp );
  }
  if( nxlow > 0 ) {
    rgamma ->gondzioProjection( rmin, rmax );
    rgamma ->selectNonZeros( *ixlow );
  }
  if( nxupp > 0 ) {
    rphi   ->gondzioProjection( rmin, rmax );
    rphi   ->selectNonZeros( *ixupp );
  }

}
  

int QpGenResiduals::validNonZeroPattern()
{
  if( nxlow > 0 && 
      ( !rv    ->matchesNonZeroPattern( *ixlow ) ||
	!rgamma->matchesNonZeroPattern( *ixlow ) ) ) {
    return 0;
  }

  if( nxupp > 0 &&
      ( !rw  ->matchesNonZeroPattern( *ixupp ) ||
	!rphi->matchesNonZeroPattern( *ixupp ) ) ) {
    return 0;
  }
  if( mclow > 0 &&
      ( !rt     ->matchesNonZeroPattern( *iclow ) ||
	!rlambda->matchesNonZeroPattern( *iclow ) ) ) {
    return 0;
  }

  if( mcupp > 0 &&
      ( !ru ->matchesNonZeroPattern( *icupp ) ||
	!rpi->matchesNonZeroPattern( *icupp ) ) ) {
    return 0;
  }
  
  return 1;
}

QpGenResiduals::~QpGenResiduals()
{
}

void QpGenResiduals::writeToStream(std::ostream& out)
{
  /*
  printf("--------------rQ\n");
  rQ->writeToStream(out);printf("---------------------------\n");
  
  printf("rA\n");
  rA->writeToStream(out);printf("---------------------------\n");
  */
  printf("rC\n");
  rC->writeToStream(out);printf("---------------------------\n");


  printf("rz\n");
  rz->writeToStream(out);printf("---------------------------\n");
  /*  
  if ( mclow > 0 ) {
    printf("rt\n");
    rt->writeToStream(out);printf("---------------------------\n");
    printf("rlambda\n");
    rlambda->writeToStream(out);printf("---------------------------\n");
  }
  if ( mcupp > 0 ) {
    printf("ru\n");
    ru->writeToStream(out);printf("---------------------------\n");
    printf("rpi\n");
    rpi->writeToStream(out);printf("---------------------------\n");
  }

  
  if( nxlow > 0 ) {
    printf("rv\n");
    rv->writeToStream(out);printf("---------------------------\n");
    printf("rgamma\n");
    rgamma->writeToStream(out);printf("---------------------------\n");
  }
  if( nxupp > 0 ) {
    printf("rw\n");
    rw->writeToStream(out);printf("---------------------------\n");
    printf("rphi\n");
    rphi->writeToStream(out);printf("---------------------------\n");
    }
  */
}
