/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#include "sTree.h"
#include "sData.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "SimpleVector.h"
#include "DoubleMatrixTypes.h"
#include "pipsport.h"

#include <cmath>

StochIterateResourcesMonitor sTree::iterMon;

int sTree::rankMe    =-1;
int sTree::rankZeroW = 0;
int sTree::rankPrcnd =-1;
int sTree::numProcs  =-1;

sTree::sTree()
: commP2ZeroW(MPI_COMM_NULL), np(-1), IPMIterExecTIME(-1)
{}


sTree::~sTree() 
{
  for(size_t it=0; it<children.size(); it++)
    delete children[it];
}

int sTree::myl() const
{
   return -1;
}

int sTree::mzl() const
{
   return -1;
}

void sTree::assignProcesses(MPI_Comm comm)
{
  int size;
  MPI_Comm_size(comm, &size);

  vector<int> processes(size);
  for(int p=0; p<size; p++)
    processes[p]=p;

  assignProcesses(comm, processes);
}

#ifndef MIN
#define MIN(a,b) ( (a>b) ? b : a )
#endif

void sTree::assignProcesses(MPI_Comm world, vector<int>& processes)
{
  int ierr;
  commWrkrs = world;
  myProcs = processes;

  int noProcs = processes.size();

  if(children.size()==0) {
    //only one CPU assigned to leafs
    assert(noProcs==1); 
    return;
  }

  if(1==noProcs) {
    for(size_t c=0; c<children.size(); c++)
      children[c]->assignProcesses(world, processes);
    return;
  }
  //here noProcs >=2 so we have to assign children to them
#if 0
  //what is the load for each children
  vector<double> vecChildNodesLoad(children.size()); 
  for(size_t i=0; i<children.size(); i++) {
    vecChildNodesLoad[i] = children[i]->processLoad();
  }
#endif
  //**** solve the asignment problem ****
  //here we'll have the mapping of children to processes
  vector<vector<int> > mapChildNodesToProcs;

  mapChildNodesToProcs.resize(children.size());
  /* old mapping
  assert(children.size() % noProcs == 0);
  for(size_t i=0; i<children.size(); i++) {
    mapChildNodesToProcs[i].resize(1);
    mapChildNodesToProcs[i][0] = i % noProcs;
  }
  */
  // new assignment to agree with BA.cpp. we'll see if this breaks anything
  int nper = children.size()/noProcs;

  /* too many MPI processes? */
  if( nper == 0 )
  {
     std::cout << "too many MPI processes! (max: " << children.size() << ")" << std::endl;
     MPI_Abort(world, 1);
  }

  for(size_t i=0; i<children.size(); i++) {
    mapChildNodesToProcs[i].resize(1);
    mapChildNodesToProcs[i][0] = MIN(i/nper,(size_t)noProcs-1);
  }
  
// #ifdef TIMING    
//   //!log
   // if(0==rankMe) {

   //   int* noduri = new int[noProcs];
   //   for(int i=0; i<noProcs; i++) noduri[i]=0;
    
   //   for(size_t i=0; i<mapChildNodesToProcs.size(); i++)
   //     noduri[mapChildNodesToProcs[i][0]]++;

   //   printf("Nodes: ");
   //   for(int i=0; i<noProcs; i++) {
   //     printf("CPU[%5d]=%5d  ", i, noduri[i]);
   //   }
   //   printf("\n");   
   //   delete[] noduri;
   // }
   //~log
   // #endif 


  MPI_Group mpiWorldGroup; 
  ierr = MPI_Comm_group(commWrkrs, &mpiWorldGroup); assert(ierr == MPI_SUCCESS);
  (void) ierr;
  for( size_t i = 0; i < children.size(); i++ )
  {
     const int noRanks4ThisChild = mapChildNodesToProcs[i].size();
     int * ranksToKeep = new int[noRanks4ThisChild];

     bool isChildInThisProcess = false;

     for( int proc = 0; proc < noRanks4ThisChild; proc++ )
     {
        ranksToKeep[proc] = mapChildNodesToProcs[i][proc];
        if( rankMe == ranksToKeep[proc] )
           isChildInThisProcess = true;
     }

     vector<int> childRanks(noRanks4ThisChild);
     for( int c = 0; c < noRanks4ThisChild; c++ )
        childRanks[c] = ranksToKeep[c];

     if( isChildInThisProcess )
     {
        if( noRanks4ThisChild == 1 )
        {
           children[i]->assignProcesses(MPI_COMM_SELF, childRanks);
        }
        else
        {
           //create the communicator this child should use
           MPI_Comm childComm;
           MPI_Group childGroup;
           ierr = MPI_Group_incl(mpiWorldGroup, noRanks4ThisChild, ranksToKeep,
                 &childGroup);
           assert(ierr==MPI_SUCCESS);

           ierr = MPI_Comm_create(commWrkrs, childGroup, &childComm);
           assert(ierr==MPI_SUCCESS);
           MPI_Group_free(&childGroup); //MPI_Group_free(&mpiWorldGroup);

           //!log printf("----Node [%d] is on proc [%d]\n", i, rankMe);fflush(stdout);
           children[i]->assignProcesses(childComm, childRanks);
        } //END noRanks4ThisChild>1
     }
     else
     { //this Child was not assigned to this CPU
        //!log printf("---Node [%d] not on  proc [%d] \n", i, rankMe);fflush(stdout);
        // continue solving the assignment problem so that
        // each node knows the CPUs the other nodes are on.

        children[i]->assignProcesses(MPI_COMM_NULL, childRanks);
     }

     delete[] ranksToKeep;
  }

  MPI_Group_free(&mpiWorldGroup);
}


double sTree::processLoad() const
{
  //! need a recursive and also a collective call
  if (IPMIterExecTIME<0.0)
    //return (NNZQ+NNZA+NNZB+NNZC+NNZD + N+MY+MZ)/1000.0;
    return (N+MY+MZ)/1000;
  return IPMIterExecTIME;
}

void sTree::GetGlobalSizes(long long& NOut, long long& MYOut, long long& MZOut)
{
  NOut=N; MYOut=MY; MZOut=MZ;
}

/*void sTree::GetLocalSizes(int& nOut, int& myOut, int& mzOut)
{
  nOut=nx(); myOut=my(); mzOut=mz();
}
*/
int sTree::innerSize(int which)
{
  if(which==0) return nx();
  if(which==1) return my();
  assert(which==2);
  return mz();
}

void sTree::syncPrimalVector(StochVector& stVec)
{
  //syncStochVector(stVec,0);
  syncStochVector(stVec);
}

void sTree::syncDualYVector(StochVector& stVec)
{
  //syncStochVector(stVec,1);
  syncStochVector(stVec);
}

void sTree::syncDualZVector(StochVector& stVec)
{
  syncStochVector(stVec);//,2);
}

void sTree::syncStochSymMatrix(StochSymMatrix& mat) 
{
  int syncChildren=0;
  char* marked4Del = new char[children.size()];

  assert(false && "Code needs update to Sync also the Hessian bordering blocks");

  for(size_t it=0; it<children.size(); it++) {
    marked4Del[it]=0;
    
    int syncNeeded; int sending; int partner; 
    children[it]->getSyncInfo(rankMe, syncNeeded, sending, partner);
    if(syncNeeded) {
      syncChildren = 1;
      if(sending) {

	marked4Del[it]=1;

	int dims[3];
	dims[0] = mat.children[it]->diag->size();
	dims[1] = mat.children[it]->diag->getStorageRef().numberOfNonZeros();
	dims[2] = mat.children[it]->n;

	MPI_Send(dims, 3, MPI_INT, partner, children[it]->id(), MPI_COMM_WORLD);
	
	MPI_Send(mat.children[it]->diag->krowM(), dims[0]+1, MPI_INT, partner, 
		 2*children[it]->id(), MPI_COMM_WORLD);
	MPI_Send(mat.children[it]->diag->jcolM(), dims[1], MPI_INT, partner, 
		 3*children[it]->id(), MPI_COMM_WORLD);
	MPI_Send(mat.children[it]->diag->M(), dims[1], MPI_DOUBLE, partner,
		 4*children[it]->id(), MPI_COMM_WORLD);
	
      } else {
	//receiving
	int dims[3]; MPI_Status status;

#ifndef NDEBUG
	const int ierr = MPI_Recv(dims, 3, MPI_INT, partner, children[it]->id(), MPI_COMM_WORLD, &status);
	assert(ierr == MPI_SUCCESS);
#else
	(void) MPI_Recv(dims, 3, MPI_INT, partner, children[it]->id(), MPI_COMM_WORLD, &status);
#endif

	if(mat.children.size() == children.size()) {
	  delete mat.children[it];
	  mat.children[it] = new StochSymMatrix(children[it]->id(), dims[2], dims[0], dims[1],
						children[it]->commWrkrs);
	} else {
	  assert(mat.children.size()==it);
	  mat.AddChild( new StochSymMatrix(children[it]->id(), dims[2], dims[0], dims[1],
					   children[it]->commWrkrs) );
	}

	MPI_Recv(mat.children[it]->diag->krowM(), dims[0]+1, MPI_INT, partner, 
		 2*children[it]->id(), MPI_COMM_WORLD, &status);
	MPI_Recv(mat.children[it]->diag->jcolM(), dims[1], MPI_INT, partner,
		 3*children[it]->id(), MPI_COMM_WORLD, &status);
	MPI_Recv(mat.children[it]->diag->M(), dims[1], MPI_DOUBLE, partner,
		 4*children[it]->id(), MPI_COMM_WORLD, &status);
      }
    }
  }
  if(syncChildren)
    for(size_t it=0; it<children.size(); it++) {
      assert(mat.children.size()==children.size());
      children[it]->syncStochSymMatrix(*mat.children[it]);
    }

  //delete the children marked for deletion
  for(size_t it=0; it<children.size(); it++) {
    if( marked4Del[it]==1) {
      delete mat.children[it];
      mat.children[it] = new StochSymDummyMatrix(children[it]->id());
      //printf(">>>>>>>>>> id %d deleted -> dummy\n", children[it]->id());
    }
  }
  delete[] marked4Del;
}


void sTree::syncStochGenMatrix(StochGenMatrix& mat) 
{
  int syncChildren=0;
  char* marked4Del = new char[children.size()];

  for(size_t it=0; it<children.size(); it++) {
    marked4Del[it]=0;
    
    int syncNeeded; int sending; int partner; 
    children[it]->getSyncInfo(rankMe, syncNeeded, sending, partner);
    if(syncNeeded) {
      syncChildren = 1;
      if(sending) {

	marked4Del[it]=1;

	int dims[8];
	mat.children[it]->Amat->getSize(*dims, *(dims+1));
	dims[2] = mat.children[it]->Amat->getStorageRef().numberOfNonZeros();
	mat.children[it]->Bmat->getSize(*(dims+3), *(dims+4));
	dims[5] = mat.children[it]->Bmat->getStorageRef().numberOfNonZeros();
	dims[6] = mat.children[it]->n; dims[7] = mat.children[it]->m;

	assert(dims[0]==dims[3]);
	MPI_Send(dims, 6, MPI_INT, partner, children[it]->id(), MPI_COMM_WORLD);
	
	MPI_Send(mat.children[it]->Amat->krowM(), dims[0]+1, MPI_INT, partner, 
		 2*children[it]->id(), MPI_COMM_WORLD);
	MPI_Send(mat.children[it]->Amat->jcolM(), dims[2], MPI_INT, partner, 
		 3*children[it]->id(), MPI_COMM_WORLD);
	MPI_Send(mat.children[it]->Amat->M(), dims[2], MPI_DOUBLE, partner,
		 4*children[it]->id(), MPI_COMM_WORLD);

	MPI_Send(mat.children[it]->Bmat->krowM(), dims[3]+1, MPI_INT, partner, 
		 5*children[it]->id(), MPI_COMM_WORLD);
	MPI_Send(mat.children[it]->Bmat->jcolM(), dims[5], MPI_INT, partner,
		 6*children[it]->id(), MPI_COMM_WORLD);
	MPI_Send(mat.children[it]->Bmat->M(), dims[5], MPI_DOUBLE, partner,
		 7*children[it]->id(), MPI_COMM_WORLD);
      } else {
	//receiving
	int dims[8]; MPI_Status status;

#ifndef NDEBUG
	const int ierr = MPI_Recv(dims, 8, MPI_INT, partner, children[it]->id(), MPI_COMM_WORLD, &status);
   assert(ierr == MPI_SUCCESS);
#else
   (void) MPI_Recv(dims, 8, MPI_INT, partner, children[it]->id(), MPI_COMM_WORLD, &status);
#endif

	assert(dims[0]==dims[3]);

	if(mat.children.size() == children.size()) {
	  delete mat.children[it];
	  mat.children[it] = new StochGenMatrix(children[it]->id(), dims[6], dims[7], 
						dims[0], dims[1], dims[2],
						dims[3], dims[4], dims[5],
						children[it]->commWrkrs);
	} else {
	  assert(mat.children.size()==it);
	  mat.AddChild( new StochGenMatrix(children[it]->id(), dims[6], dims[7], 
					   dims[0], dims[1], dims[2],
					   dims[3], dims[4], dims[5],
					   children[it]->commWrkrs) );
	}

	MPI_Recv(mat.children[it]->Amat->krowM(), dims[0]+1, MPI_INT, partner, 
		 2*children[it]->id(), MPI_COMM_WORLD, &status);
	MPI_Recv(mat.children[it]->Amat->jcolM(), dims[2], MPI_INT, partner,
		 3*children[it]->id(), MPI_COMM_WORLD, &status);
	MPI_Recv(mat.children[it]->Amat->M(), dims[2], MPI_DOUBLE, partner,
		 4*children[it]->id(), MPI_COMM_WORLD, &status);

	MPI_Recv(mat.children[it]->Bmat->krowM(), dims[3]+1, MPI_INT, partner, 
		 5*children[it]->id(), MPI_COMM_WORLD, &status);
	MPI_Recv(mat.children[it]->Bmat->jcolM(), dims[5], MPI_INT, partner,
		 6*children[it]->id(), MPI_COMM_WORLD, &status);
	MPI_Recv(mat.children[it]->Bmat->M(), dims[5], MPI_DOUBLE, partner,
		 7*children[it]->id(), MPI_COMM_WORLD, &status);
      }
    }
  }
  if(syncChildren)
    for(size_t it=0; it<children.size(); it++) {
      assert(mat.children.size()==children.size());
      children[it]->syncStochGenMatrix(*mat.children[it]);
    }

  //delete the children marked for deletion
  for(size_t it=0; it<children.size(); it++) {
    if( marked4Del[it]==1) {
      delete mat.children[it];
      mat.children[it] = new StochGenDummyMatrix(children[it]->id());
    }
  }
  delete[] marked4Del;
}

void sTree::syncStochVector(StochVector& stVec)
{
  int syncChildren=0;

  char* marked4Del = new char[children.size()];
  for(size_t it=0; it<children.size(); it++) {
    marked4Del[it] = 0;

    int syncNeeded; int sending; int partner; 
    children[it]->getSyncInfo(rankMe, syncNeeded, sending, partner);
    if(syncNeeded) {
      syncChildren = 1;
      if(sending) {
	// this child needs to be deleted
	marked4Del[it] = 1;
	
	SimpleVector* vec = (SimpleVector*)stVec.children[it]->vec;
	int        dim = vec->length();
	double* buffer = vec->elements();

	//send the dimension
	MPI_Send(&dim, 1, MPI_INT, partner, children[it]->id(), MPI_COMM_WORLD);
	//printf("send %d\n", dim);
	if(dim>0) {
	  //send

#ifndef NDEBUG
	  const int ierr = MPI_Send(buffer, dim, MPI_DOUBLE, partner, children[it]->id(), MPI_COMM_WORLD);
	  assert(ierr == MPI_SUCCESS);
#else
	  (void) MPI_Send(buffer, dim, MPI_DOUBLE, partner, children[it]->id(), MPI_COMM_WORLD);
#endif
	} else { /*printf("zero length vector, NO actual SYNC\n");*/ }
	
	
      } else {
	//receive
	int dim; MPI_Status status;
	MPI_Recv(&dim, 1, MPI_INT, partner, children[it]->id(), MPI_COMM_WORLD, &status);
	//printf("children size: old=%d new=%d\n",
	//stVec.children.size(), children.size());
	//printf("recv %d\n", dim);
	double* buffer;
	if(stVec.children.size() == children.size()) {
	  assert(stVec.children[it]->isKindOf(kStochDummy));
	  delete stVec.children[it];
	  
	  stVec.children[it] = new StochVector(dim, children[it]->commWrkrs);
	  stVec.children[it]->parent = &stVec;
	  buffer = ((SimpleVector*)stVec.children[it]->vec)->elements();
	} else { //children will be added one by one to stVec (it has
		 //no children)
	  assert(stVec.children.size()==it);

	  stVec.AddChild(new StochVector(dim, children[it]->commWrkrs));
	  buffer = ((SimpleVector*)stVec.children[stVec.children.size()-1]->vec)->elements();
	}
	//do nothing for the stoch empty vectors
	if( dim > 0 ) {
#ifndef NDEBUG
	  const int ierr = MPI_Recv(buffer, dim, MPI_DOUBLE, partner, children[it]->id(), MPI_COMM_WORLD, &status);
	  assert(ierr == MPI_SUCCESS);
#else
	  (void) MPI_Recv(buffer, dim, MPI_DOUBLE, partner, children[it]->id(), MPI_COMM_WORLD, &status);
#endif

	} else { /*printf("zero length vector, NO SYNC\n");*/ }
      }
    } else {
      //printf("CPU[%d] Nod[%d] sync NOT NEEDED\n", rankMe, children[it]->id());
    }  
  }


  if(syncChildren)
    for(size_t it=0; it<children.size(); it++) {
      assert(stVec.children.size()==children.size());
      children[it]->syncStochVector(*stVec.children[it]);
    }

  //delete the children marked for deletion
  for(size_t it=0; it<children.size(); it++) {
    if( marked4Del[it]==1) {
      delete stVec.children[it];
      stVec.children[it] = new StochDummyVector();
    }
  }

  delete[] marked4Del;
}

StochVector* sTree::newPrimalVector() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* x = new StochVector(nx(), commWrkrs);
  assert(x!=nullptr);

  for(size_t it=0; it<children.size(); it++) {
    StochVector* child = children[it]->newPrimalVector();
    x->AddChild(child);
  }
  return x;
}

StochVector* sTree::newDualYVector() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  //length of linking part
  int yl = (np == -1) ? myl() : -1;

  StochVector* y = new StochVector(my(), yl, commWrkrs, -1);

  for(size_t it=0; it<children.size(); it++) {
    StochVector* child = children[it]->newDualYVector();
    y->AddChild(child);
  }
  return y;
}

StochVector* sTree::newDualZVector() const
{
  //is this node a dead-end for this process?
  if(commWrkrs == MPI_COMM_NULL)
    return new StochDummyVector();

  //length of linking part

  int zl = (np == -1) ? mzl() : -1;

  StochVector* z = new StochVector(mz(), zl, commWrkrs, -1);

  for(size_t it=0; it<children.size(); it++) {
    StochVector* child = children[it]->newDualZVector();
    z->AddChild(child);
  }
  return z;
}


StochVector* sTree::newPrimalVectorEmpty() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* x = new StochVector(0, commWrkrs);

  for(size_t it=0; it<children.size(); it++) {
    StochVector* child = children[it]->newPrimalVector();
    x->AddChild(child);
  }
  return x;
}

StochVector* sTree::newDualYVectorEmpty() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* y = new StochVector(0, commWrkrs);

  for(size_t it = 0; it < children.size(); it++) {
    StochVector* child = children[it]->newDualYVector();
    y->AddChild(child);
  }
  return y;
}

StochVector* sTree::newDualZVectorEmpty() const
{
  //is this node a dead-end for this process?
  if(commWrkrs == MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* z = new StochVector(0, commWrkrs);

  for(size_t it = 0; it < children.size(); it++) {
    StochVector* child = children[it]->newDualZVector();
    z->AddChild(child);
  }
  return z;
}


StochVector* sTree::newRhs()
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  int locmyl = (np == -1) ? myl() : 0;
  int locmzl = (np == -1) ? mzl() : 0;

  locmyl = max(locmyl, 0);
  locmzl = max(locmzl, 0);

  StochVector* rhs = new StochVector(nx() + my() + mz() + locmyl + locmzl, commWrkrs);

  for(size_t it=0; it<children.size(); it++) {
    StochVector* child = children[it]->newRhs();
    rhs->AddChild(child);
  }

  return rhs;
}

void sTree::startMonitors()
{
  iterMon.recIterateTm_start();
  startNodeMonitors();
}

void sTree::stopMonitors()
{
  iterMon.recIterateTm_stop();
  stopNodeMonitors();
}

void sTree::startNodeMonitors()
{
  resMon.reset();
  for(size_t i=0; i<children.size(); i++)
    children[i]->startNodeMonitors();
}

void sTree::stopNodeMonitors()
{
  for(size_t i=0; i<children.size(); i++)
    children[i]->stopNodeMonitors();

  resMon.computeTotal();
}

void sTree::toMonitorsList(list<NodeExecEntry>& lstExecTm)
{
  lstExecTm.push_back(resMon.eTotal);

  for(size_t i=0; i<children.size(); i++)
    children[i]->toMonitorsList(lstExecTm);
}

void sTree::fromMonitorsList(list<NodeExecEntry>& lstExecTm)
{
  resMon.eTotal = lstExecTm.front();
  lstExecTm.pop_front();

  for(size_t i=0; i<children.size(); i++)
    children[i]->fromMonitorsList(lstExecTm);
}

void sTree::syncMonitoringData(vector<double>& vCPUTotal)
{

  list<NodeExecEntry> lstExecTm;
  this->toMonitorsList(lstExecTm);

  int noNodes = lstExecTm.size(); 
  int nCPUs   = vCPUTotal.size();

  double* recvBuf = new double[noNodes+nCPUs];
  double* sendBuf = new double[noNodes+nCPUs];
  
  list<NodeExecEntry>::iterator iter = lstExecTm.begin();
  for(int it=0; it<noNodes; it++) { sendBuf[it] = iter->tmChildren; iter++; }
  for(int it=noNodes; it<noNodes+nCPUs; it++) sendBuf[it] = vCPUTotal[it-noNodes];

  
  MPI_Allreduce(sendBuf, recvBuf, noNodes+nCPUs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  iter = lstExecTm.begin();
  for(int it=0; it<noNodes; it++) { iter->tmChildren = recvBuf[it]; iter++; }
  for(int it=noNodes; it<noNodes+nCPUs; it++) vCPUTotal[it-noNodes] = recvBuf[it];

  if(children.size()>0) {
    //local time MPI_MAX, but only for nonleafs
    iter = lstExecTm.begin();
    for(int it=0; it<noNodes; it++) { sendBuf[it] = iter->tmLocal; iter++; }

    MPI_Allreduce(sendBuf, recvBuf, noNodes, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    iter = lstExecTm.begin();
    for(int it=0; it<noNodes; it++) { iter->tmLocal = recvBuf[it]; iter++; }
  }
  delete[] recvBuf; delete[] sendBuf;

  //populate the tree with the global data
  this->fromMonitorsList(lstExecTm);

  //compute syncronized total time for each node of the tree, i.e., local+childs+subchilds
  computeNodeTotal(); //updates this->IPMIterExecTIME
}


bool sTree::balanceLoad()
{
  return false; //disabled for now
  //before synchronization, compute the total time recorded on this CPU
  //updates this->IPMIterExecTIME
  computeNodeTotal();

  int nCPUs; MPI_Comm_size(commWrkrs, &nCPUs);
  vector<double> cpuExecTm(nCPUs, 0.0);
  cpuExecTm[rankMe] = this->IPMIterExecTIME;

  //if(sleepFlag && rankMe==1) {cpuExecTm[1] += 4.0;}

  this->syncMonitoringData(cpuExecTm);

  //!log
#if defined(STOCH_TESTING) && !defined(TIMING)
  if(!rankMe) {
    printf("Iteration times per process:\n");
    for(int it=0; it<nCPUs; it++) printf("%8.4f ", cpuExecTm[it]);
    printf("\n\n");
    this->displayExecTimes(0);
  }
#endif
  if(nCPUs==1) return 0;

  double total = 0.0; double maxLoad=0, minLoad=1.e+10;
  for(int it=0; it<nCPUs; it++) {
    total += cpuExecTm[it];
    if(maxLoad<cpuExecTm[it]) maxLoad = cpuExecTm[it];
    if(minLoad>cpuExecTm[it]) minLoad = cpuExecTm[it];
  }

  if(maxLoad<1.0) return 0;

  double balance = max(maxLoad/total*nCPUs, total/nCPUs/minLoad);
  if(balance<1.3) { //it is OK, no balancing

    //decide if balancing is needed due to the 'fat' nodes
    total = 0.0; maxLoad=0;
    for(size_t i=0; i<children.size(); i++) {
      total += children[i]->IPMIterExecTIME;
      
      //if(!rankMe) printf("%7.4f ", children[i]->IPMIterExecTIME);
      
      if(maxLoad<children[i]->IPMIterExecTIME) 
	maxLoad = children[i]->IPMIterExecTIME;
    }
    //if(!rankMe) printf("\n");//aaa
    
    balance = maxLoad/total*children.size();

    //if(!rankMe) printf("CPU node balance=%g\n", balance);

    if(balance<2.00) return 0;
    //else if(!rankMe) printf("!!! Balancing NEEDED: node balance=%g\n", balance);
  } else {
    //if(!rankMe) printf("!!! Balancing NEEDED: CPU balance=%g\n", balance);
  }
  return 0;
  //save the current MPI related information
  saveCurrentCPUState();

  cpuExecTm.clear();

  vector<int> ranks(nCPUs);
  for(int i=0; i<nCPUs; i++) ranks[i]=i; 
  assignProcesses(MPI_COMM_WORLD, ranks);
  return 1;
}

#define maSend 1
#define maRecv 0

void sTree::getSyncInfo(int rank, int& syncNeeded, int& sendOrRecv, int& toFromCPU )
{
  // was this node previously assigned to cpu 'rank'?
  int  wasAssigned=isInVector(rank, myOldProcs);
  // is currently assigned to cpu 'rank'?
  int  isAssigned=isInVector(rank, myProcs);

  //if(0==wasAssigned && 0==isAssigned) return;
  //if(1==wasAssigned && 1==isAssigned) return;
  syncNeeded=0;  sendOrRecv = 0; toFromCPU = -1;
  if(isAssigned!=wasAssigned) {
    syncNeeded=1;
    
    if(wasAssigned) {
      assert(0==isAssigned); assert(myProcs.size()>0);
      //where is this node assigned?
      toFromCPU=myProcs[0];
      sendOrRecv = maSend;
    } else {
      assert(1==isAssigned); assert(myOldProcs.size()>0);
      toFromCPU=myOldProcs[0];
      sendOrRecv = maRecv;
    }
  }
}

int sTree::isInVector(int elem, const vector<int>& vec)
{
  for(size_t i=0; i<vec.size(); i++)
    if(elem==vec[i]) return 1;
  return 0;
}

void sTree::computeNodeTotal()
{
  if(0==children.size())
    this->IPMIterExecTIME = resMon.eTotal.tmChildren;
  else {

    this->IPMIterExecTIME = resMon.eTotal.tmLocal;
    for(size_t i=0; i<children.size(); i++) {
      children[i]->computeNodeTotal();
      
      this->IPMIterExecTIME += children[i]->IPMIterExecTIME;
    }
  }
}

void sTree::saveCurrentCPUState()
{
  myOldMpiComm = commWrkrs;
  myOldProcs = myProcs;

  for(size_t i=0; i<children.size(); i++)
    children[i]->saveCurrentCPUState();
}

#ifdef DEADCODE // STOCH_TESTING
void sTree::displayProcessInfo(int rank)
{
  char szTab[5000];
  for(int i=0; i<5000; i++) szTab[i]='\0';
  if(rankMe==rank)
    displayProcessInfo(szTab);
}

void sTree::displayProcessInfo(char* tab)
{
  printf("%sNode[%d] -> Procs", tab, id());
  for(size_t i=0; i<myProcs.size(); i++) {
    printf("%4d ", myProcs[i]);
  }
  printf("\n");

  int len=strlen(tab);
  tab[len] = tab[len+1] = tab[len+2] = tab[len+3] = ' ';

  for(size_t i=0; i<children.size(); i++)
    children[i]->displayProcessInfo(tab);

  tab[len]= '\0';
}

void sTree::displayVectorVsTreeStructure(StochVector& stVec, int rank)
{
  char szTab[5000];
  for(int i=0; i<5000; i++) szTab[i]='\0';
 
  if(rankMe==rank)
    displayVectorVsTreeStructure(stVec, rankMe, szTab);
}

void sTree::displayVectorVsTreeStructure(StochVector& stVec, int myRank, char* tab)
{
  char szType[50];
  if(stVec.isKindOf(kStochVector)) sprintf(szType, "stoch sz%d", stVec.vec->length());
  else if(stVec.isKindOf(kStochDummy)) strcpy(szType, "dummy");
  else strcpy(szType, "unkn");

  char szShouldBe[50];
  if(isInVector(myRank, myProcs))
    strcpy(szShouldBe, "stoch");
  else 
    strcpy(szShouldBe, "dummy");

  printf("%sNode[%d] on %d should be %s  |  Vec %s    childs %d\n", 
	 tab, id(), myRank, szShouldBe, szType, stVec.children.size());

  if(stVec.isKindOf(kStochDummy)) return;

  int len=strlen(tab);
  tab[len] = tab[len+1] = tab[len+2] = tab[len+3] = ' ';

  for(size_t i=0; i<children.size(); i++)
    children[i]->displayVectorVsTreeStructure(*stVec.children[i], myRank, tab);

  tab[len]= '\0';
}

void sTree::displayMatVsTreeStructure(StochGenMatrix& mat, int rank)
{
  char szTab[5000];
  for(int i=0; i<5000; i++) szTab[i]='\0';
 
  if(rankMe==rank)
    displayMatVsTreeStructure(mat, rankMe, szTab);
}

void sTree::displayMatVsTreeStructure(StochGenMatrix& mat, int myRank, char* tab)
{

  char szType[50];
  if(mat.isKindOf(kStochGenMatrix)) strcpy(szType, "stoch");
  else if(mat.isKindOf(kStochGenDummyMatrix)) strcpy(szType, "dummy");
  else strcpy(szType, "unkn");

  char szShouldBe[50];
  if(isInVector(myRank, myProcs))
    strcpy(szShouldBe, "stoch");
  else 
    strcpy(szShouldBe, "dummy");

  printf("%sNode[%d] on %d should be %s  |  Mat %s    childs %d \n", 
	 tab, id(), myRank, szShouldBe, szType, mat.children.size());

  if(mat.isKindOf(kStochGenDummyMatrix)) return;

  int len=strlen(tab);
  tab[len] = tab[len+1] = tab[len+2] = tab[len+3] = ' ';

  for(size_t i=0; i<children.size(); i++)
    children[i]->displayMatVsTreeStructure(*mat.children[i], myRank, tab);

  tab[len]= '\0';
}

void sTree::displayMatVsTreeStructure(StochSymMatrix& mat, int rank)
{
  char szTab[5000];
  for(int i=0; i<5000; i++) szTab[i]='\0';
 
  if(rankMe==rank)
    displayMatVsTreeStructure(mat, rankMe, szTab);
}

void sTree::displayMatVsTreeStructure(StochSymMatrix& mat, int myRank, char* tab)
{
  char szType[50];
  if(mat.isKindOf(kStochSymMatrix)) strcpy(szType, "stoch");
  else if(mat.isKindOf(kStochSymDummyMatrix)) strcpy(szType, "dummy");
  else strcpy(szType, "unkn");

  char szShouldBe[50];
  if(isInVector(myRank, myProcs))
    strcpy(szShouldBe, "stoch");
  else 
    strcpy(szShouldBe, "dummy");

  printf("%sNode[%d] on %d should be %s  |  id[%d] Mat %s   size %d    childs %d\n", 
	 tab, id(), myRank, szShouldBe, mat.id, szType, mat.mat->size(), mat.children.size());

  if(mat.isKindOf(kStochSymDummyMatrix)) return;

  int len=strlen(tab);
  tab[len] = tab[len+1] = tab[len+2] = tab[len+3] = ' ';

  for(size_t i=0; i<children.size(); i++)
    children[i]->displayMatVsTreeStructure(*mat.children[i], myRank, tab);

  tab[len]= '\0';
}


void sTree::displayExecTimes(int rank)
{
  char szTab[5000];
  for(int i=0; i<5000; i++) szTab[i]='\0';

  if(rankMe==rank) displayExecTimes(szTab);
}

void sTree::displayExecTimes(char* tab)
{
  resMon.computeTotal();
  printf("%sNod[%4d]-server=%8.6fsec  child=%8.6fsec  total=%8.6f  schurMult=%8.6f\n", 
	 tab, id(),
	 resMon.eTotal.tmLocal, resMon.eTotal.tmChildren,
	 IPMIterExecTIME, 
	 np==-1?resMon.eMult.tmLocal:resMon.eMult.tmChildren);

	 //	 resMon.eFact.tmLocal,  resMon.eLsolve.tmLocal);
	 //	 resMon.eDsolve.tmLocal, resMon.eLtsolve.tmLocal,
	 //resMon.eFact.tmChildren,  resMon.eLsolve.tmChildren, 
	 //resMon.eDsolve.tmChildren, resMon.eLtsolve.tmChildren);
	 
  int len=strlen(tab);
  tab[len] = tab[len+1] = tab[len+2] = tab[len+3] = ' ';
  
  for(size_t i=0; i<children.size(); i++)
    children[i]->displayExecTimes(tab);
  
  tab[len]= 0;
}
#endif

