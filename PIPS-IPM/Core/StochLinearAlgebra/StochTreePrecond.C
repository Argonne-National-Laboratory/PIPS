#include "StochTreePrecond.h"
#include "StochResourcePlanner.h"

#include <iostream>
#include <assert.h>

using namespace std;

StochTreePrecond::StochTreePrecond(StochInputTree* root)
  : StochTree(root), commWorkers(MPI_COMM_NULL)
{ }
StochTreePrecond::StochTreePrecond()
  : StochTree(), commWorkers(MPI_COMM_NULL)
{ }

StochTreePrecond::~StochTreePrecond()
{ }



void StochTreePrecond::assignProcesses()
{
  int iErr; int* ranks;
  MPI_Group worldGroup;  MPI_Group P2SWGroup;  MPI_Group WrkGroup;

  int size; MPI_Comm_size(MPI_COMM_WORLD, &size); assert(size>=2);

  vector<int> processes(size); for(int p=0; p<size; p++) processes[p]=p;
  rankPrcnd = size-1; rankZeroW = 0;

  iErr = MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);
  
  ranks = new int[size-1];
  //////////////////////////////////////////////////////////////////
  // create a communicator with the special worker and the precond
  ////////////////////////////////////////////////////////////////// 
  ranks[0]=rankZeroW; ranks[1]=rankPrcnd;
  //create the group
  iErr = MPI_Group_incl(worldGroup, 2, ranks, &P2SWGroup);
  //create the new communicator
  iErr = MPI_Comm_create(MPI_COMM_WORLD, P2SWGroup, &commP2ZeroW); assert(iErr==MPI_SUCCESS);
  //free the group
  MPI_Group_free(&P2SWGroup);

  //////////////////////////////////////////////////////////////////
  // create a communicator with the  workers
  ////////////////////////////////////////////////////////////////// 
  for(int i=0; i<size-1; i++) ranks[i]=i;
  //create the group
  iErr = MPI_Group_incl(worldGroup, size-1, ranks, &WrkGroup);
  //create the new communicator
  iErr = MPI_Comm_create(MPI_COMM_WORLD, WrkGroup, &commWorkers); assert(iErr==MPI_SUCCESS);
  //free the group
  MPI_Group_free(&WrkGroup);


  delete[] ranks;
  assignProcesses(MPI_COMM_WORLD, processes);
}


void StochTreePrecond::assignProcesses(MPI_Comm workersComm, vector<int>& procs)
{
  int ierr;
  
  commWrkrs = workersComm;
  myProcs = procs;
 
  int noProcs = procs.size();

  //only one CPU assigned to leafs
  if(children.size()==0) { assert(noProcs==1); return; }

  //if(1==noProcs) {
  //  for(int c=0; c<children.size(); c++)
  //    children[c]->assignProcessesPrecond(workersComm, procs);
  //  return;
  //}
    
  ///////////////////////////////////////////////////////////
  // here noProcs >=2 so we have to assign children to them
  ///////////////////////////////////////////////////////////

  //what is the load for each children
  vector<double> vChildsLoad(children.size()); 
  for (size_t i=0; i<children.size(); i++) {
    vChildsLoad[i] = children[i]->processLoad();
  }

  //////////////////////////////////////////////////////////////
  // solve the asignment problem 
  //////////////////////////////////////////////////////////////
  vector<vector<int> > mChilds2Procs;
  vector<int> procsNoPr(noProcs-1);
  for(int i=0; i<noProcs; i++)
    if(procs[i]!=rankPrcnd) procsNoPr[i] = procs[i];

  /////////////////////////////////////////////////////////////
  // assign child nodes to 0,1,...,P-1
  /////////////////////////////////////////////////////////////
  StochResourcePlanner planner; double balance;
  planner.assignProcesses(procsNoPr, vChildsLoad,
			  mChilds2Procs, balance); 
  assert(vChildsLoad.size() == mChilds2Procs.size());

  /////////////////////////////////////////////////////////////
  // 'manually" assign child nodes to P
  /////////////////////////////////////////////////////////////
  assignToPreconditioner(vChildsLoad, mChilds2Procs);

  
  //!log --------------------------------------------
  if(0==rankMe) {
    int* noduri = new int[noProcs]; int noCopii=children.size();
    for(int i=0; i<noProcs; i++) noduri[i]=0;
    
    for(int i=0; i<noCopii; i++) noduri[mChilds2Procs[i][0]]++;
    
    printf("Nodes: ");
    for(int i=0; i<noProcs; i++) printf("CPU[%3d]=%3d  ", i, noduri[i]);
    printf("\n");  delete[] noduri;
  }
  //~log --------------------------------------------

  /////////////////////////////////////////////////////////////
  // recursively create communicators for each child
  ////////////////////////////////////////////////////////////

  MPI_Group mpiWorldGroup; 
  ierr = MPI_Comm_group(MPI_COMM_WORLD, &mpiWorldGroup); assert(ierr==MPI_SUCCESS);
  for (size_t i=0; i<children.size(); i++) {

    int isChildInProc=0;
    int noRanksOfChild = mChilds2Procs[i].size();
    int* ranksToKeep = new int[noRanksOfChild];

    for(int proc=0; proc<noRanksOfChild; proc++) {
      ranksToKeep[proc] = mChilds2Procs[i][proc];
      if(rankMe==ranksToKeep[proc])
	isChildInProc=1;
    }   

    vector<int> childRanks(noRanksOfChild);
    for(int c=0; c<noRanksOfChild; c++) childRanks[c]=ranksToKeep[c];

    if(isChildInProc) {
   
      if(noRanksOfChild==1) {

	children[i]->assignProcesses(MPI_COMM_SELF,childRanks); 

      } else {
	//create the communicator this child should use
	MPI_Comm  childComm; MPI_Group childGroup;
	ierr = MPI_Group_incl(mpiWorldGroup, noRanksOfChild, ranksToKeep, &childGroup);
	assert(ierr==MPI_SUCCESS);
	delete[] ranksToKeep;	
	
	ierr = MPI_Comm_create(mpiWorldGroup, childGroup, &childComm); assert(ierr==MPI_SUCCESS);
	MPI_Group_free(&childGroup); //MPI_Group_free(&mpiWorldGroup);
	
	//!log printf("----Node [%d] is on proc [%d]\n", i, rankMe);fflush(stdout); 
	children[i]->assignProcesses(childComm,childRanks); 

      } //END noRanksOfChild>1

    } else { //this Child was not assigned to this CPU
      delete[] ranksToKeep; 
      //!log printf("---Node [%d] not on  proc [%d] \n", i, rankMe);fflush(stdout);
      // continue solving the assignment problem so that 
      // each node knows the CPUs the other nodes are on.
      
      children[i]->assignProcesses(MPI_COMM_NULL,childRanks); 
    }
  }
  MPI_Group_free(&mpiWorldGroup);
}

void 
StochTreePrecond::assignToPreconditioner(vector<double>& vChildsLoad, 
					 vector<vector<int> >& mChilds2Procs)
{
  // do the "manual" assignment for the preconditioner process by taking 
  // from workers and assign to preconditioner

  int takefromall=0;
  if(!takefromall) { 
    //*****************************************
    // take one scenario from the zero-worker
    //*****************************************
    int rankToTakeFrom=rankZeroW; 
    // find the minimum-load scenario from the zero-worker process
    int minPos=-1; double minVal=1.0e9; int totalChilds=0;

    do {
      for(size_t nod=0; nod<mChilds2Procs.size(); nod++) {
	
	if(mChilds2Procs[nod][0]==rankToTakeFrom) {
	  totalChilds++;
	  if(vChildsLoad[nod]<minVal) {
	    minPos = nod; minVal = vChildsLoad[nod];
	  }
	}
      } assert(minPos>=0);

      if(totalChilds>1) break;

      rankToTakeFrom++; minPos=-1; minVal=1.0e9; totalChilds=0;

      if(rankToTakeFrom>=rankPrcnd) {
	cout << "PANIC!!! StochTreePrecond can NOT assign work to PRCND proc.\n" << 
	  "The number of scenarios is too small for the requested number of processes.\n"; 
	break;
      }
    } while(true);

    //for 'minPos' child replace the rank 'rankZeroW' of zero-worker  
    //                   with    the rank 'rankPrcnd' of preconditioner
    mChilds2Procs[minPos][0] = rankPrcnd;
  }
}
/*bool StochTree::balanceLoadPrecond()
{
  //recursively get total the execution time
  computeNodeTotal();

  int nCPUs; MPI_Comm_size(MPI_COMM_WORLD, &nCPUs);
  vector<double> cpuExecTm(nCPUs, 0.0);
  cpuExecTm[rankMe] = this->IPMIterExecTIME;

  //if(sleepFlag && rankMe==1) {cpuExecTm[1] += 4.0;}

  this->syncMonitoringData(cpuExecTm);

  //!log
  if(!rankMe) {
    printf("Iteration times per process:\n");
    for(int it=0; it<nCPUs; it++) printf("%8.4f ", cpuExecTm[it]);
    printf("\n\n");
#ifdef STOCH_TESTING
    this->displayExecTimes(0);
#endif
  }
  if(nCPUs<=10) return 0;
  
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
    for(int i=0; i<children.size(); i++) {
      total += children[i]->IPMIterExecTIME;
      
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

  //save the current MPI related information
  saveCurrentCPUState();

  cpuExecTm.clear();

  vector<int> ranks(nCPUs);
  for(int i=0; i<nCPUs; i++) ranks[i]=i;
  assignProcsPrecond();//(MPI_COMM_WORLD, ranks);
  return 1;
}
*/
