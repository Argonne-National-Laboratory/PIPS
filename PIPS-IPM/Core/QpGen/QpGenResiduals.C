/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifdef PRINT_MAX
#include "StochVector.h"
#include "SimpleVector.h"
#include <cmath>
#endif

#include "QpGenResiduals.h"
#include "QpGenVars.h"
#include "QpGenData.h"

#include "OoqpVector.h"
#include "LinearAlgebraPackage.h"

#include <iostream>
#include <fstream>



using namespace std;

#include "mpi.h"

QpGenResiduals::QpGenResiduals( LinearAlgebraPackage * la,
				long long nx_, long long my_, long long mz_,
				OoqpVector * ixlow_in, OoqpVector * ixupp_in,
				OoqpVector * iclow_in, OoqpVector * icupp_in )
{
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

void QpGenResiduals::calcresids(Data *prob_in, Variables *vars_in)
{
  int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  QpGenVars * vars = (QpGenVars *) vars_in;
  QpGenData * prob = (QpGenData *) prob_in;

  double componentNorm, norm=0.0, gap=0.0;
 
  prob->getg( *rQ );
  prob->Qmult( 1.0, *rQ,  1.0, *vars->x );

  // calculate x^T (g+Qx) - contribution to the duality gap
  gap = rQ->dotProductWith(*vars->x); 

  prob->ATransmult( 1.0, *rQ, -1.0, *vars->y );
  prob->CTransmult( 1.0, *rQ, -1.0, *vars->z );

  vars->gamma->selectNonZeros(*ixlow);
  vars->phi->selectNonZeros( *ixupp );
  if( nxlow > 0 ) rQ->axpy( -1.0, *vars->gamma );
  if( nxupp > 0 ) rQ->axpy(  1.0, *vars->phi );

  componentNorm = rQ->infnorm();
#ifdef TIMING
  double rQtwonorm=rQ->twonorm();
  if(0==myRank)  cout << " rQ infnorm=" << componentNorm 
		      << " | twonorm=" << rQtwonorm << endl;
#endif
  if( componentNorm > norm ) norm = componentNorm;

  prob->getbA( *rA );
  prob->Amult( -1.0, *rA, 1.0, *vars->x );

  //printf("gap2=%20.16f\n", gap);
  
  // contribution -d^T y to duality gap
  gap -= prob->bA->dotProductWith(*vars->y);

  componentNorm = rA->infnorm();
#ifdef TIMING
  if(0==myRank) cout << " rA norm = " << componentNorm << endl;
#endif
  if( componentNorm > norm ) norm = componentNorm;

  rC->copyFrom( *vars->s );
  prob->Cmult( -1.0, *rC, 1.0, *vars->x );

  componentNorm = rC->infnorm();
#ifdef TIMING
  if(0==myRank) cout << " rC norm = " << componentNorm << endl;
#endif
  // cout << " rC norm = " << componentNorm << endl;
  // if( componentNorm > norm ) norm = componentNorm;

  rz->copyFrom( *vars->z );

  if( mclow > 0 ) {
    rz->axpy( -1.0, *vars->lambda );

#ifdef TIMING
    if(0==myRank) cout << "lambda norm " << (vars->lambda)->infnorm() << std::endl;
#endif
    rt->copyFrom( *vars->s );
    rt->axpy( -1.0, prob->slowerBound() );
    rt->selectNonZeros( *iclow );
    rt->axpy( -1.0, *vars->t );
    gap -= prob->bl->dotProductWith(*vars->lambda);

    componentNorm = rt->infnorm();
#ifdef TIMING
    if(0==myRank) cout << " rt norm = " << componentNorm << endl;
#endif
    //cout << " rt norm = " << componentNorm << endl;
    if( componentNorm > norm ) norm = componentNorm;
  }
  if( mcupp > 0 ) { 
    rz->axpy(  1.0, *vars->pi );
#ifdef TIMING
    if(0==myRank) cout  << "pi norm " << (vars->pi)->infnorm() << std::endl;
#endif

#ifdef PRINT_MAX
        double max = 0.0;

        int maxchild = -1;
        int maxindex = -1;
        {

        OoqpVector& lp = *vars->pi;
        StochVector& lam = reinterpret_cast<StochVector&>(lp);
        std::cout << "childrensize " << lam.children.size() << std::endl;
        for( int k = 0; k < lam.children.size(); k++ )
        {
           StochVector& lamc = *lam.children[k];

           SimpleVector& lvec = reinterpret_cast<SimpleVector&>(*lamc.vec);
           for( int i = 0; i < lvec.length(); i++ )
           {
           if( fabs(lvec[i]) > max )
           {
              maxindex = i;
              maxchild = k;
              max = (lvec)[i];
           }
           }
        }
        std::cout << "max at child " << maxchild << " ind " << maxindex << ": " << max << std::endl;

        }


        {

        OoqpVector& lp = *vars->u;
        StochVector& lam = reinterpret_cast<StochVector&>(lp);
        for( int k = 0; k < lam.children.size(); k++ )
        {
           StochVector& lamc = *lam.children[k];

           SimpleVector& lvec = reinterpret_cast<SimpleVector&>(*lamc.vec);
           for( int i = 0; i < lvec.length(); i++ )
           {
              if( i == maxindex && k == maxchild )
              {
                 std::cout << "uval at child " << k  << " ind: " << i << ": " << lvec[i] << std::endl;
              }
           }
        }
        }
#endif

    ru->copyFrom( *vars->s );
    ru->axpy( -1.0, prob->supperBound() );
    ru->selectNonZeros( *icupp );
    ru->axpy( 1.0, *vars->u );

    gap += prob->bu->dotProductWith(*vars->pi);

    componentNorm = ru->infnorm();
#ifdef TIMING
    if(0==myRank) cout << " ru norm = " << componentNorm << endl;
#endif
    //    cout << " ru norm = " << componentNorm << endl;
    if( componentNorm > norm ) norm = componentNorm;
  }
  componentNorm = rz->infnorm();
  //  cout << " rz norm = " << componentNorm << endl;
#ifdef TIMING
  if(0==myRank) cout << " rz norm = " << componentNorm << endl;
#endif
  if( componentNorm > norm ) norm = componentNorm;

  if( nxlow > 0 ) {
    rv->copyFrom( *vars->x );
    rv->axpy( -1.0, prob->xlowerBound() );
    rv->selectNonZeros( *ixlow );
    rv->axpy( -1.0, *vars->v );

    gap -= prob->blx->dotProductWith(*vars->gamma);

    componentNorm = rv->infnorm();
    //    cout << " rv norm = " << componentNorm << endl;
#ifdef TIMING
    if(0==myRank) cout << " rv norm = " << componentNorm << endl;
#endif
    if( componentNorm > norm ) norm = componentNorm;
  }
  if( nxupp > 0 ) {
    rw->copyFrom( *vars->x );
    rw->axpy( -1.0, prob->xupperBound() );
    rw->selectNonZeros( *ixupp );
    rw->axpy(  1.0, *vars->w );

    gap += prob->bux->dotProductWith(*vars->phi);

    componentNorm = rw->infnorm();
#ifdef TIMING
    if(0==myRank) cout << " rw norm = " << componentNorm << endl;
#endif
    //    cout << " rw norm = " << componentNorm << endl;
    if( componentNorm > norm ) norm = componentNorm;
  }
   
  mDualityGap = gap;
  mResidualNorm = norm; 
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


void QpGenResiduals::writeToStream(ostream& out)
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
