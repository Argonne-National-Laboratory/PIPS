#include "RemoteMatTimesVec.h"
#include "SimpleVector.h"
#include "pipsport.h"
#include "mpi.h"

/************************************************************************
 * Implementation of the class RemoteMatTimesVec which perform the 
 * multiplication of a remote matrix with a local vector.
 * 
 */

RemoteMatTimesVec::RemoteMatTimesVec(StochTreePrecond* tree_)
  : stochNode(tree_), tmpVec1(nullptr), delTmpVec1(0)
{};

RemoteMatTimesVec::RemoteMatTimesVec()
  : tmpVec1(nullptr), delTmpVec1(0)
{ }

RemoteMatTimesVec::~RemoteMatTimesVec()
{ 
  if(delTmpVec1==1 && tmpVec1!=nullptr) delete[] tmpVec1;
};

 /** y = beta * y + alpha * A * x */
void RemoteMatTimesVec::doIt(double beta, OoqpVector& y, 
			     double alpha, OoqpVector& x)
{
  int n = y.length(); assert(x.length()==n);

  if(tmpVec1==nullptr) { tmpVec1 = new double[n]; delTmpVec1=1; }

  int mpiComm     = stochNode->commWrkrs;
  int rankPrcnd   = stochNode->rankPrcnd;
  int rankZeroW   = stochNode->rankZeroW;
  int rankMe      = stochNode->rankMe;
  
  MPI_Status status;  

  SimpleVector Ax(tmpVec1, n); Ax.copyFrom(x);
  tmpVec1[n+PREC_EXTRA_DOUBLES-1]=0; // set convergence flag to false

  //printf("RemoteTimesVec rank=[%d] sending %d doubles\n", rankMe, n+PREC_EXTRA_DOUBLES);fflush(stdout);
  MPI_Send(tmpVec1, n+PREC_EXTRA_DOUBLES, MPI_DOUBLE, rankPrcnd, 1, mpiComm);
  //printf("RemoteTimesVec rank=[%d] SENT! receiving %d doubles... \n", rankMe, n);fflush(stdout);
  MPI_Recv(&Ax[0],  n,       MPI_DOUBLE, rankPrcnd, 2, mpiComm, &status);
  //printf("RemoteTimesVec rank=[%d] RECEIVED! finalizing. \n", rankMe);fflush(stdout);
  if(beta==0) {

    if(alpha==1.0) y.copyFrom(Ax);
    else { Ax.scale(alpha); y.copyFrom(Ax); }

  } else {
    y.scale(beta);
    y.axpy(alpha,Ax);  
  }
};

void RemoteMatTimesVec::setTmpVec1(double* buf)
{
  if(tmpVec1!=nullptr && delTmpVec1==1) delete[] tmpVec1;
  delTmpVec1 = 0;
  tmpVec1=buf;
}
