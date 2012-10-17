/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2010 Argonne National Laboratory. See Copyright Notification.  */

#include "StochTree.h"
#include "QpGenStoch.h"
#include "sData.h"

#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "SimpleVector.h"
#include "DoubleMatrixTypes.h"
#include "StochResourcePlanner.h"
#include <cmath>

using namespace std;

#ifdef DEBUG
extern int sleepFlag;
#endif


#ifndef UCTRANS // see note in smlParDriver.C
#define UCTRANS
#endif

StochIterateResourcesMonitor StochTree::iterMon;

int StochTree::rankMe    =-1;
int StochTree::rankZeroW = 0;
int StochTree::rankPrcnd =-1;
int StochTree::numProcs  =-1;

StochTree::StochTree() 
  : commP2ZeroW(MPI_COMM_NULL), data(NULL), tree(NULL), fakedata(NULL), 
     np(-1), IPMIterExecTIME(-1)
{
  if(-1==rankMe) MPI_Comm_rank(MPI_COMM_WORLD, &rankMe);
  if(-1==numProcs) MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
}

StochTree::StochTree(StochInputTree* inputTree)
  : commP2ZeroW(MPI_COMM_NULL), tree(NULL), fakedata(NULL), 
    np(-1), IPMIterExecTIME(-1)
{
  if(-1==rankMe) MPI_Comm_rank(MPI_COMM_WORLD, &rankMe);
  if(-1==numProcs) MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  data = inputTree->nodeInput;
#ifndef POOLSCEN
  for(size_t it=0; it<inputTree->children.size(); it++)
    children.push_back(new StochTree(inputTree->children[it]));
#else
	tree = inputTree;
#endif
}

// np==-1 is used to indicate the root node. these can't be root nodes
StochTree::StochTree(const vector<StochInputTree::StochInputNode*> &localscens)
  : commP2ZeroW(MPI_COMM_NULL), data(NULL), tree(NULL), scens(localscens), 
    np(0), IPMIterExecTIME(-1) 
{
  if(-1==rankMe) MPI_Comm_rank(MPI_COMM_WORLD, &rankMe);
  if(-1==numProcs) MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	fakedata = new StochInputTree::StochInputNode();
	real_children.reserve(scens.size());
	for(size_t i = 0; i < scens.size(); i++) {
		real_children.push_back(new StochTree(scens[i]));

	}
}
StochTree::StochTree(StochInputTree::StochInputNode* data_)
  : commP2ZeroW(MPI_COMM_NULL), data(data_), tree(NULL), fakedata(NULL),
    np(0), IPMIterExecTIME(-1)
{
  if(-1==rankMe) MPI_Comm_rank(MPI_COMM_WORLD, &rankMe);
  if(-1==numProcs) MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
}



StochTree::~StochTree() 
{
  for(size_t it=0; it<children.size(); it++)
    delete children[it];
	if (fakedata) delete fakedata;
	for(size_t i = 0; i < real_children.size(); i++)
		delete real_children[i];
}

// this is usually called before assigning processes
void StochTree::computeGlobalSizes()
{
  if (data) {
    N  = data->n;
    MY = data->my;
    MZ = data->mz;
    
    NNZQ = data->nnzQ;
    NNZA = data->nnzA;
    NNZB = data->nnzB;
    NNZC = data->nnzC;
    NNZD = data->nnzD;
  } else {
    N = MY = MZ = NNZQ = NNZA = NNZB = NNZC = NNZD = 0;
  }
  if (tree && np == -1) {
    for(size_t it=0; it<tree->children.size();it++) {
      N += tree->children[it]->nodeInput->n;
      MY += tree->children[it]->nodeInput->my;
      MZ += tree->children[it]->nodeInput->mz;
      
      NNZQ += tree->children[it]->nodeInput->nnzQ;
      NNZA += tree->children[it]->nodeInput->nnzA;
      NNZB += tree->children[it]->nodeInput->nnzB;
      NNZC += tree->children[it]->nodeInput->nnzC;
      NNZD += tree->children[it]->nodeInput->nnzD;
    }
  } else if (fakedata) {
    fakedata->n = fakedata->my = fakedata->mz = fakedata->nnzQ = fakedata->nnzA = fakedata->nnzB = fakedata->nnzC = fakedata->nnzD = 0;	
    for(size_t it=0; it<scens.size();it++) {
      fakedata->n += scens[it]->n;
      fakedata->my += scens[it]->my;
      fakedata->mz += scens[it]->mz;
      fakedata->nnzQ += scens[it]->nnzQ;
      fakedata->nnzA += scens[it]->nnzA;
      fakedata->nnzB += scens[it]->nnzB;
      fakedata->nnzC += scens[it]->nnzC;
      fakedata->nnzD += scens[it]->nnzD;
      real_children[it]->np = np;
    }
    N += fakedata->n;
    MY += fakedata->my;
    MZ += fakedata->mz;
    NNZQ += fakedata->nnzQ;
    NNZA += fakedata->nnzA;
    NNZB += fakedata->nnzB;
    NNZC += fakedata->nnzC;
    NNZD += fakedata->nnzD;
  }
  for(size_t it=0; it<children.size(); it++) {
    children[it]->np = this->data->n;
    children[it]->computeGlobalSizes();
    N  += children[it]->N;
    MY += children[it]->MY;
    MZ += children[it]->MZ;
    
    //nnz stuff
    NNZQ += children[it]->NNZQ;
    NNZA += children[it]->NNZA;
    NNZB += children[it]->NNZB;
    NNZC += children[it]->NNZC;
    NNZD += children[it]->NNZD;
  }
}


void StochTree::assignProcesses()
{
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  vector<int> processes(size);
  for(int p=0; p<size; p++)
    processes[p]=p;

  assignProcesses(MPI_COMM_WORLD, processes);
}


#ifndef POOLSCEN
void StochTree::assignProcesses(MPI_Comm world, vector<int>& processes)
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

  //what is the load for each children
  vector<double> vecChildNodesLoad(children.size()); 
  for(size_t i=0; i<children.size(); i++) {
    vecChildNodesLoad[i] = children[i]->processLoad();
  }

  //if(sleepFlag) vecChildNodesLoad[0] = 4*vecChildNodesLoad[50];

  //**** solve the asignment problem ****
  //here we'll have the mapping of children to processes
  vector<vector<int> > mapChildNodesToProcs;
#ifndef UCTRANS  
  StochResourcePlanner planner; double balance;
  planner.assignProcesses(processes, vecChildNodesLoad,
			  mapChildNodesToProcs, balance);
  // override the planner in the trivial case of equal loads
  bool eq = true;
  for(size_t i=0; i < children.size(); i++) {
    if (fabs(vecChildNodesLoad[i]-vecChildNodesLoad[0]) > 1E-4) {
      eq = false; break;
    }
  }
  if (eq && (children.size() % noProcs == 0)) {
    for(size_t i=0; i<mapChildNodesToProcs.size();i++) {
      mapChildNodesToProcs[i][0] = i % noProcs;
    }
  }
#else
  assert(children.size() % noProcs == 0);
  mapChildNodesToProcs.resize(children.size());
  for(size_t i=0; i<children.size(); i++) {
    mapChildNodesToProcs[i].resize(1);
    mapChildNodesToProcs[i][0] = i % noProcs;
  }
#endif

//   if(0==rankMe) {
//     //printf("CPU[%d] processes reassigned -> balance=%8.3e\n", rankMe, balance);
//     //printMap(mapChildNodesToProcs, "----------\n", "[%d][%d]=%d  ");
//     //!log
//     //printMap(mapChildNodesToProcs, "----------\n", "[%d][%d]=%d  "); 
//     int* noduri = new int[noProcs];
//     for(int i=0; i<noProcs; i++) noduri[i]=0;
    
//     for(size_t i=0; i<mapChildNodesToProcs.size(); i++)
//       noduri[mapChildNodesToProcs[i][0]]++;
    
//     printf("Nodes: ");
//     for(int i=0; i<noProcs; i++) {
//       printf("CPU[%5d]=%5d  ", i, noduri[i]);
//     }
//     printf("\n");
    
//     delete[] noduri;
//   }
//   //~log

  MPI_Group mpiWorldGroup; 
  ierr = MPI_Comm_group(MPI_COMM_WORLD, &mpiWorldGroup); assert(ierr==MPI_SUCCESS);
  for(size_t i=0; i<children.size(); i++) {

    int isChildInThisProcess=0;
    int noRanks4ThisChild = mapChildNodesToProcs[i].size();
    int * ranksToKeep = new int[noRanks4ThisChild];

    for(int proc=0; proc<noRanks4ThisChild; proc++) {
      ranksToKeep[proc] = mapChildNodesToProcs[i][proc];
      if(rankMe==ranksToKeep[proc])
	isChildInThisProcess=1;
    }   

    vector<int> childRanks(noRanks4ThisChild);
    for(int c=0; c<noRanks4ThisChild; c++) childRanks[c]=ranksToKeep[c];

    if(isChildInThisProcess) {
   
      if(noRanks4ThisChild==1) {
	children[i]->assignProcesses(MPI_COMM_SELF,childRanks); 
      } else {
	//create the communicator this child should use
	MPI_Comm  childComm; MPI_Group childGroup;
	ierr = MPI_Group_incl(mpiWorldGroup, noRanks4ThisChild, ranksToKeep, &childGroup);
	assert(ierr==MPI_SUCCESS);
	delete[] ranksToKeep;
	
	
	ierr = MPI_Comm_create(mpiWorldGroup, childGroup, &childComm); assert(ierr==MPI_SUCCESS);
	MPI_Group_free(&childGroup); //MPI_Group_free(&mpiWorldGroup);
	
	//!log printf("----Node [%d] is on proc [%d]\n", i, rankMe);fflush(stdout); 
	children[i]->assignProcesses(childComm,childRanks); 
      } //END noRanks4ThisChild>1

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

#else
// Combine scenarios
void StochTree::assignProcesses(MPI_Comm world, vector<int>& processes)
{
  int ierr;
  commWrkrs = world;
  myProcs = processes;

  int noProcs = processes.size();

  if (1 == noProcs && tree==NULL) {
    for(size_t i = 0; i < real_children.size(); i++) {
      real_children[i]->assignProcesses(commWrkrs, myProcs);
    }
    return;
  }
  
  //here noProcs >=2 so we have to assign children to them
  
  //what is the load for each children
  vector<double> vecChildNodesLoad(tree->children.size()); 
  for(size_t i=0; i<tree->children.size(); i++) {
    //vecChildNodesLoad[i] = data->children[i]->processLoad();
  	vecChildNodesLoad[i] = 1;
	}

  //if(sleepFlag) vecChildNodesLoad[0] = 4*vecChildNodesLoad[50];

  //**** solve the asignment problem ****
  //here we'll have the mapping of children to processes
  vector<vector<int> > mapChildNodesToProcs;
#ifndef UCTRANS  
  StochResourcePlanner planner; double balance;
  planner.assignProcesses(processes, vecChildNodesLoad,
			  mapChildNodesToProcs, balance);
  // override the planner in the trivial case of equal loads
  bool eq = true;
  for(size_t i=0; i < tree->children.size(); i++) {
    if (fabs(vecChildNodesLoad[i]-vecChildNodesLoad[0]) > 1E-4) {
      eq = false; break;
    }
  }
  if (eq && (tree->children.size() % noProcs == 0)) {
    for(size_t i=0; i<mapChildNodesToProcs.size();i++) {
      mapChildNodesToProcs[i][0] = i % noProcs;
    }
  }
#else
  assert(tree->children.size() % noProcs == 0);
  mapChildNodesToProcs.resize(tree->children.size());
  for(size_t i=0; i<tree->children.size(); i++) {
    mapChildNodesToProcs[i].resize(1);
    mapChildNodesToProcs[i][0] = i % noProcs;
  }
#endif


//   if(0==rankMe) {
//     //printf("CPU[%d] processes reassigned -> balance=%8.3e\n", rankMe, balance);
//     //printMap(mapChildNodesToProcs, "----------\n", "[%d][%d]=%d  ");
//     //!log
//     //printMap(mapChildNodesToProcs, "----------\n", "[%d][%d]=%d  "); 
//     int* noduri = new int[noProcs];
//     for(int i=0; i<noProcs; i++) noduri[i]=0;
    
//     for(size_t i=0; i<mapChildNodesToProcs.size(); i++)
//       noduri[mapChildNodesToProcs[i][0]]++;
    
//     printf("Nodes: ");
//     for(int i=0; i<noProcs; i++) {
//       printf("CPU[%5d]=%5d  ", i, noduri[i]);
//     }
//     printf("\n");
    
//     delete[] noduri;
//   }
//   //~log

	
  vector<vector<int> > mapProcsToChildNodes(noProcs);
  
  for(size_t i=0; i<tree->children.size(); i++) {
    
    int noRanks4ThisChild = mapChildNodesToProcs[i].size();
    assert(noRanks4ThisChild == 1); // don't support splitting a child over procs
    int childRank = mapChildNodesToProcs[i][0];
    mapProcsToChildNodes[childRank].push_back(i);
  }
  assert(children.size() == 0);
  children.reserve(noProcs);
  for(int p=0; p<noProcs;p++) {
    vector<StochInputTree::StochInputNode*> scenarios(mapProcsToChildNodes[p].size());
    for(size_t i=0; i<mapProcsToChildNodes[p].size();i++) {
      scenarios[i] = tree->children[mapProcsToChildNodes[p][i]]->nodeInput;
    }
    children.push_back(new StochTree(scenarios));
    vector<int> s(1,p);
    if (p == rankMe) {
      children[p]->assignProcesses(MPI_COMM_SELF, s);
    } else {
      children[p]->assignProcesses(MPI_COMM_NULL, s);
    }
  }
  tree = NULL;
  computeGlobalSizes(); // recompute so children nodes have values
}
#endif




/*void StochTree::assignProcesses(MPI_Comm myWorld)
{   
  commWrkrs = myWorld;
  assert(false);
  //return if the null communicator was passed, this means this node
  //and its subnodes do not belong to this process.
  if(commWrkrs==MPI_COMM_NULL)
    return;

  //also leaf nodes don't assign processes
  if (children.size()==0)
    return;

  int noProcs;
  int ierr = MPI_Comm_size(myWorld, &noProcs);
  assert(noProcs>=1);

  //if this tree non-leaf node runs in only one process then so do the children
  if(1==noProcs) {
    for(size_t i=0; i<children.size(); i++)
      children[i]->assignProcesses(myWorld);
    return;
  }

  //here noProcs >=2 so we have to assign children to them

  //what is the load for each children
  vector<double> vecChildNodesLoad(children.size()); 
  for(size_t i=0; i<children.size(); i++) {
    vecChildNodesLoad[i] = children[i]->processLoad();
  }
  
  // **** solve the asignment problem ****
  //here we'll have the mapping of children to processes
  vector<vector<int> > mapChildNodesToProcs;
  
  //solve the assignment problem
  StochResourcePlanner planner;
  planner.assignProcesses(noProcs, vecChildNodesLoad,
			  mapChildNodesToProcs); 

  // create a new communicator for each child
  // and delegate further CPUs assignment to children

  //get the MPI_Group of current world
  MPI_Group myGroup;
  ierr = MPI_Comm_group(myWorld, &myGroup);

  //also get the rank of the current process
 
  for(size_t i=0; i<children.size(); i++) {
    //create a sub-group by using MPI_Group_incl and assign it to
    //child IF THE CHILD BELONGS TO THIS PROCESS.

    //unfortunately we need a temporary C-vector here

    int isChildInThisProcess=0;

    int noRanks = mapChildNodesToProcs[i].size();
    int * ranksToKeep = new int[noRanks];

    for(int proc=0; proc<noRanks; proc++) {
      ranksToKeep[proc] = mapChildNodesToProcs[i][proc];
      if(rankMe==ranksToKeep[proc])
	isChildInThisProcess=1;
    }    

    if(isChildInThisProcess) {

      printf("node %d got cpu %d\n", i, rankMe);

      MPI_Comm  childComm;
      MPI_Group childGroup;

      ierr = MPI_Group_incl(myGroup, noRanks, ranksToKeep, &childGroup);
      delete[] ranksToKeep;
    

      //create the communicator and give it to the child
      //cout << id() << " calling MPI_Comm_create for children " << children[i]->id() << endl;
      ierr = MPI_Comm_create(myWorld, childGroup, &childComm);
      if(ierr!=MPI_SUCCESS) {
	cout << "Error " << ierr << " in MPI_Comm_create in process " << rankMe << endl;
      }

      MPI_Group_free(&childGroup);
      
      children[i]->assignProcesses(childComm);

    } else {
      delete[] ranksToKeep;
      children[i]->assignProcesses(MPI_COMM_NULL);
    }
  }
  MPI_Group_free(&myGroup);
}
*/

double StochTree::processLoad() const
{
  //! need a recursive and also a collective call
  if (IPMIterExecTIME<0.0)
    return (NNZQ+NNZA+NNZB+NNZC+NNZD + N+MY+MZ)/1000.0;
  return IPMIterExecTIME;
}

void StochTree::GetGlobalSizes(int& NOut, int& MYOut, int& MZOut)
{
  NOut=N; MYOut=MY; MZOut=MZ;
}

void StochTree::GetLocalSizes(int& nOut, int& myOut, int& mzOut)
{
  nOut=nx(); myOut=my(); mzOut=mz();
}

StochSymMatrix* StochTree::createQ() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochSymDummyMatrix(id());
  
  if (!fakedata) {
    if(data->nnzQ<0)
      data->fnnzQ(data->user_data, data->id, &data->nnzQ);
    
    StochSymMatrix* Q = 
      new StochSymMatrix(data->id, 
			 N, 
			 data->n, 
			 data->nnzQ,
			 commWrkrs);
    
    data->fQ(data->user_data, data->id, 
	     Q->mat->krowM(), 
	     Q->mat->jcolM(),
	     Q->mat->M());  
    
    for(size_t it=0; it<children.size(); it++) {
      StochSymMatrix* child = children[it]->createQ();
      Q->AddChild(child);
    }
    return Q;
  } else {
    assert(real_children.size() > 0);
    vector<StochSymMatrix*> v(real_children.size());
    for(size_t i = 0; i<real_children.size(); i++) {
      v[i] = real_children[i]->createQ();
    }
    StochSymMatrix *out = new StochSymMatrix(v);
    for(size_t i = 0; i<real_children.size(); i++) delete v[i];
    return out;
  }
}




StochGenMatrix* StochTree::createA() const
{

  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL) {
    return new StochGenDummyMatrix(id());
  }

  StochGenMatrix* A = NULL;
  if (!fakedata) {
    if (np==-1) {

      //data->fnnzA(data->user_data, data->id, &nnzA);
      if (data->nnzA<0)
        data->fnnzA(data->user_data, data->id, &data->nnzA);
      data->nnzB=0;

      //this is the root; populate B with A's data
      //B_0 is the A_0 from the theoretical form
      A = new StochGenMatrix(data->id, 
           N, MY, 
           data->my, np, data->nnzB,
           data->my, data->n,  data->nnzA,
           commWrkrs);
      //populate the submatrices A and B
      //data->fA(data->user_data, data->id, A->Amat->krowM(), A->Amat->jcolM(), A->Amat->M());
      data->fA(data->user_data, data->id, A->Bmat->krowM(), A->Bmat->jcolM(), A->Bmat->M());
    } else {

      if (data->nnzA<0)
        data->fnnzA(data->user_data, data->id, &data->nnzA);
      if (data->nnzB<0)
        data->fnnzB(data->user_data, data->id, &data->nnzB);

      A = new StochGenMatrix(data->id, 
           N, MY, 
           data->my, np, data->nnzA, 
           data->my, data->n,  data->nnzB,
           commWrkrs);
      //populate the submatrices A and B
      data->fA(data->user_data, data->id, A->Amat->krowM(), A->Amat->jcolM(), A->Amat->M());
      data->fB(data->user_data, data->id, A->Bmat->krowM(), A->Bmat->jcolM(), A->Bmat->M());
    }

    for(size_t it=0; it<children.size(); it++) {
      StochGenMatrix* child = children[it]->createA();
      A->AddChild(child);
    }
  } else {
    vector<StochGenMatrix*> v(real_children.size());
#ifdef UCTRANS
    v[0] = real_children[0]->createA();
    for(size_t i = 1; i<real_children.size();i++) v[i] = v[0];
#else
    for(size_t i = 0; i<real_children.size(); i++) {
      v[i] = real_children[i]->createA();
    }
#endif
    A = new StochGenMatrix(v);
#ifdef UCTRANS
    delete v[0];
#else
    for(size_t i = 0; i<real_children.size(); i++) delete v[i];
#endif

  }
  return A;
}

StochGenMatrix* StochTree::createC() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochGenDummyMatrix(id());


  StochGenMatrix* C = NULL;
  if (!fakedata) {
    if (np==-1) {

      if(data->nnzC<0)
        data->fnnzC(data->user_data, data->id, &data->nnzC);
      data->nnzD=0;
      
      //this is the root; populate D with C's data
      //D_0 is the C_0 from the theoretical form
      C = new StochGenMatrix(data->id, 
           N, MZ, 
           data->mz, np, data->nnzD, 
           data->mz, data->n,  data->nnzC, commWrkrs);

      data->fC(data->user_data, data->id, C->Bmat->krowM(), C->Bmat->jcolM(), C->Bmat->M());

    } else {
      if(data->nnzC<0)
        data->fnnzC(data->user_data, data->id, &data->nnzC);
      if(data->nnzD<0)
        data->fnnzD(data->user_data, data->id, &data->nnzD);

      C = new StochGenMatrix(data->id, 
           N, MZ, 
           data->mz, np, data->nnzC, 
           data->mz, data->n,  data->nnzD,
           commWrkrs);

      //populate the submatrices C and D
      data->fC(data->user_data, data->id, C->Amat->krowM(), C->Amat->jcolM(), C->Amat->M());
      data->fD(data->user_data, data->id, C->Bmat->krowM(), C->Bmat->jcolM(), C->Bmat->M());
    }
      
    for(size_t it=0; it<children.size(); it++) {
      StochGenMatrix* child = children[it]->createC();
      C->AddChild(child);
    }
  } else {
    vector<StochGenMatrix*> v(real_children.size());
#ifdef UCTRANS
    v[0] = real_children[0]->createC();
    for(size_t i = 1; i<real_children.size();i++) v[i] = v[0];
#else
    for(size_t i = 0; i<real_children.size(); i++) {
      v[i] = real_children[i]->createC();
    }
#endif
    C = new StochGenMatrix(v);
#ifdef UCTRANS
    delete v[0];
#else
    for(size_t i = 0; i<real_children.size(); i++) delete v[i];
#endif
  }
  return C;
}


int StochTree::innerSize(int which)
{
  if(which==0) return nx();
  if(which==1) return my();
  assert(which==2);
  return mz();
}

int StochTree::nx() const {
  if (data) return data->n;
  else return fakedata->n;
}

int StochTree::my() const {
  if (data) return data->my;
  else return fakedata->my;
}

int StochTree::mz() const {
  if (data) return data->mz;
  else return fakedata->mz;
}
// not sure what this is used for
int StochTree::id() const {
  if (data) return data->id;
  else return 0;
}



void StochTree::syncPrimalVector(StochVector& stVec)
{
  //syncStochVector(stVec,0);
  syncStochVector(stVec);
}

void StochTree::syncDualYVector(StochVector& stVec)
{
  //syncStochVector(stVec,1);
  syncStochVector(stVec);
}

void StochTree::syncDualZVector(StochVector& stVec)
{
  syncStochVector(stVec);//,2);
}

void StochTree::syncStochSymMatrix(StochSymMatrix& mat) 
{
  int ierr; int syncChildren=0;
  char* marked4Del = new char[children.size()];

  for(size_t it=0; it<children.size(); it++) {
    marked4Del[it]=0;
    
    int syncNeeded; int sending; int partner; 
    children[it]->getSyncInfo(rankMe, syncNeeded, sending, partner);
    if(syncNeeded) {
      syncChildren = 1;
      if(sending) {

	marked4Del[it]=1;

	int dims[3];
	dims[0] = mat.children[it]->mat->size();
	dims[1] = mat.children[it]->mat->getStorage()->numberOfNonZeros();
	dims[2] = mat.children[it]->n;

	MPI_Send(dims, 3, MPI_INT, partner, children[it]->id(), MPI_COMM_WORLD);
	
	MPI_Send(mat.children[it]->mat->krowM(), dims[0]+1, MPI_INT, partner, 
		 2*children[it]->id(), MPI_COMM_WORLD);
	MPI_Send(mat.children[it]->mat->jcolM(), dims[1], MPI_INT, partner, 
		 3*children[it]->id(), MPI_COMM_WORLD);
	MPI_Send(mat.children[it]->mat->M(), dims[1], MPI_DOUBLE, partner,
		 4*children[it]->id(), MPI_COMM_WORLD);

      } else {
	//receiving
	int dims[3]; MPI_Status status;
	ierr = MPI_Recv(dims, 3, MPI_INT, partner, children[it]->id(), MPI_COMM_WORLD, &status);

	if(mat.children.size() == children.size()) {
	  delete mat.children[it];
	  mat.children[it] = new StochSymMatrix(children[it]->id(), dims[2], dims[0], dims[1],
						children[it]->commWrkrs);
	} else {
	  assert(mat.children.size()==it);
	  mat.AddChild( new StochSymMatrix(children[it]->id(), dims[2], dims[0], dims[1],
					   children[it]->commWrkrs) );
	}

	MPI_Recv(mat.children[it]->mat->krowM(), dims[0]+1, MPI_INT, partner, 
		 2*children[it]->id(), MPI_COMM_WORLD, &status);
	MPI_Recv(mat.children[it]->mat->jcolM(), dims[1], MPI_INT, partner,
		 3*children[it]->id(), MPI_COMM_WORLD, &status);
	MPI_Recv(mat.children[it]->mat->M(), dims[1], MPI_DOUBLE, partner,
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


void StochTree::syncStochGenMatrix(StochGenMatrix& mat) 
{
  int ierr; int syncChildren=0;
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
	dims[2] = mat.children[it]->Amat->getStorage()->numberOfNonZeros();
	mat.children[it]->Bmat->getSize(*(dims+3), *(dims+4));
	dims[5] = mat.children[it]->Bmat->getStorage()->numberOfNonZeros();
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
	ierr = MPI_Recv(dims, 8, MPI_INT, partner, children[it]->id(), MPI_COMM_WORLD, &status);

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

void StochTree::syncStochVector(StochVector& stVec)
{
  int ierr; int syncChildren=0;

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
	  ierr = MPI_Send(buffer, dim, MPI_DOUBLE, 
			  partner, children[it]->id(), MPI_COMM_WORLD);
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

	  ierr = MPI_Recv(buffer, dim, MPI_DOUBLE, 
			  partner, children[it]->id(), MPI_COMM_WORLD, &status);
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


void StochTree::syncStochVector_old(StochVector& stVec, int which)
{
  int ierr; int syncChildren=0;
  //printf("CPU[%d] Nod[%d] syncVector  number children;%d\n",  
  //rankMe, id(), children.size());


  char* marked4Del = new char[children.size()];
  for(size_t it=0; it<children.size(); it++) {
    marked4Del[it] = 0;

    int syncNeeded; int sending; int partner; 
    children[it]->getSyncInfo(rankMe, syncNeeded, sending, partner);
    if(syncNeeded) {
      syncChildren = 1;
      if(sending) {
	double* buffer = ((SimpleVector*)stVec.children[it]->vec)->elements();

	//do nothing for the stoch empty vectors
	if( stVec.vec->length() > 0 ) {
	  //send
	  ierr = MPI_Send(buffer, children[it]->innerSize(which), MPI_DOUBLE, 
			  partner, children[it]->id(), MPI_COMM_WORLD);
	} else { /*printf("zero length vector, NO SYNC\n");*/}
	// this child needs to be deleted
	marked4Del[it] = 1;
      } else {
	//receive
	//printf("children size: old=%d new=%d\n",
	//stVec.children.size(), children.size());

	double* buffer;
	if(stVec.children.size() == children.size()) {
	  assert(stVec.children[it]->isKindOf(kStochDummy));
	  delete stVec.children[it];
	  if( stVec.vec->length() > 0 ) {
	    stVec.children[it] = new StochVector(children[it]->innerSize(which), 
						 children[it]->commWrkrs);
	  } else {
	    stVec.children[it] = new StochVector(0, children[it]->commWrkrs);

	  }
	  buffer = ((SimpleVector*)stVec.children[it]->vec)->elements();
	} else { //children will be added one by one to stVec (it has
		 //no children)
	  assert(stVec.children.size()==it);

	  stVec.AddChild(new StochVector(children[it]->innerSize(which), 
					 children[it]->commWrkrs));
	  buffer = ((SimpleVector*)stVec.children[stVec.children.size()-1]->vec)->elements();
	}
	//do nothing for the stoch empty vectors
	if( stVec.vec->length() > 0 ) {
#ifdef DEBUG
	if( children[it]->innerSize(which)!=stVec.children[it]->vec->length())
	  printf("sizes does not matches:[%d][%d]\n", 
		 children[it]->innerSize(which),
		 stVec.children[it]->vec->length());
#endif
	  MPI_Status status;
	  ierr = MPI_Recv(buffer, children[it]->innerSize(which), MPI_DOUBLE, 
			  partner, children[it]->id(), MPI_COMM_WORLD, &status);
	} else { /*printf("zero length vector, NO SYNC\n");*/ }
      }
    } else {
      //printf("CPU[%d] Nod[%d] sync NOT NEEDED\n", rankMe, children[it]->id());
    }  
  }


  if(syncChildren)
    for(size_t it=0; it<children.size(); it++) {
      assert(stVec.children.size()==children.size());
      children[it]->syncStochVector_old(*stVec.children[it], which);
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


StochVector* StochTree::createc() const
{

  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();  

  StochVector* c = new StochVector(nx(), commWrkrs);
  double* vData = ((SimpleVector*)c->vec)->elements();  
  if (!fakedata) {
    // populate the node's data with data from user.

    data->fc(data->user_data, data->id, 
       vData, data->n);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createc();
      c->AddChild(child);
    }
  } else {
#ifdef UCTRANS
    int n = scens[0]->n;
    scens[0]->fc(scens[0]->user_data,scens[0]->id,
        vData, scens[0]->n);
    for(size_t i = 1; i < scens.size(); i++) {
      memcpy(vData+i*n,vData,n*sizeof(double));
    }
#else
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fc(scens[i]->user_data,scens[i]->id,
        vData+pos, scens[i]->n);
      pos += scens[i]->n;
    }
#endif
  }

  return c;
}

StochVector* StochTree::createb() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* b = new StochVector(my(), commWrkrs);
  double* vData = ((SimpleVector*)b->vec)->elements();  
  if (!fakedata) {
    data->fb(data->user_data, data->id, 
       vData, data->my);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createb();
      b->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fb(scens[i]->user_data,scens[i]->id,
        vData+pos, scens[i]->my);
      pos += scens[i]->my;
    }
  }

  return b;
}

StochVector* StochTree::createxlow() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* xlow = new StochVector(nx(), commWrkrs);
  double* vData = ((SimpleVector*)xlow->vec)->elements();  
  if (!fakedata) {
    data->fxlow(data->user_data, data->id, 
       vData, data->n);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createxlow();
      xlow->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fxlow(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->n);
      pos += scens[i]->n;
    }
  }
  return xlow;
}

StochVector* StochTree::createixlow() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* ixlow = new StochVector(nx(), commWrkrs);
  double* vData = ((SimpleVector*)ixlow->vec)->elements();  
  if (!fakedata) {
    data->fixlow(data->user_data, data->id, 
       vData, data->n);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createixlow();
      ixlow->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fixlow(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->n);
      pos += scens[i]->n;
    }
  }
  return ixlow;
}

StochVector* StochTree::createxupp() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* xupp = new StochVector(nx(), commWrkrs);
  double* vData = ((SimpleVector*)xupp->vec)->elements();  
  if (!fakedata) {
    data->fxupp(data->user_data, data->id, 
       vData, data->n);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createxupp();
      xupp->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fxupp(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->n);
      pos += scens[i]->n;
    }
  }
  return xupp;
}

StochVector* StochTree::createixupp() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* ixupp = new StochVector(nx(), commWrkrs);
  double* vData = ((SimpleVector*)ixupp->vec)->elements();  
  if (!fakedata) {
    data->fixupp(data->user_data, data->id, 
       vData, data->n);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createixupp();
      ixupp->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fixupp(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->n);
      pos += scens[i]->n;
    }

  }

  return ixupp;
}

StochVector* StochTree::createclow() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* clow = new StochVector(mz(), commWrkrs);
  double* vData = ((SimpleVector*)clow->vec)->elements();  
  if (!fakedata) {  
    data->fclow(data->user_data, data->id, 
         vData, data->mz);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createclow();
      clow->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fclow(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->mz);
      pos += scens[i]->mz;
    }
  }
  return clow;
}

StochVector* StochTree::createiclow() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* iclow = new StochVector(mz(), commWrkrs);
  double* vData = ((SimpleVector*)iclow->vec)->elements();  
  if (!fakedata) {
    data->ficlow(data->user_data, data->id, 
       vData, data->mz);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createiclow();
      iclow->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->ficlow(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->mz);
      pos += scens[i]->mz;
    }
  }
  return iclow;
}

StochVector* StochTree::createcupp() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* cupp = new StochVector(mz(), commWrkrs);
  double* vData = ((SimpleVector*)cupp->vec)->elements();  
    if (!fakedata) {
    data->fcupp(data->user_data, data->id, 
       vData, data->mz);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createcupp();
      cupp->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->fcupp(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->mz);
      pos += scens[i]->mz;
    }

  }
  return cupp;
}

StochVector* StochTree::createicupp() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* icupp = new StochVector(mz(), commWrkrs);
  double* vData = ((SimpleVector*)icupp->vec)->elements();  
  if (!fakedata) {
    data->ficupp(data->user_data, data->id, 
       vData, data->mz);

    for(size_t it=0; it<children.size(); it++) {
      StochVector* child = children[it]->createicupp();
      icupp->AddChild(child);
    }
  } else {
    int pos = 0;
    for(size_t i = 0; i < scens.size(); i++) {
      scens[i]->ficupp(scens[i]->user_data, scens[i]->id,
        vData+pos, scens[i]->mz);
      pos += scens[i]->mz;
    }
  }
  return icupp;
}


StochVector* StochTree::newPrimalVector() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* x = new StochVector(nx(), commWrkrs);

  for(size_t it=0; it<children.size(); it++) {
    StochVector* child = children[it]->newPrimalVector();
    x->AddChild(child);
  }
  return x;
}

StochVector* StochTree::newDualYVector() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* y = new StochVector(my(), commWrkrs);

  for(size_t it=0; it<children.size(); it++) {
    StochVector* child = children[it]->newDualYVector();
    y->AddChild(child);
  }
  return y;
}

StochVector* StochTree::newDualZVector() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* z = new StochVector(mz(), commWrkrs);

  for(size_t it=0; it<children.size(); it++) {
    StochVector* child = children[it]->newDualZVector();
    z->AddChild(child);
  }
  return z;
}


StochVector* StochTree::newPrimalVectorEmpty() const
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

StochVector* StochTree::newDualYVectorEmpty() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* y = new StochVector(0, commWrkrs);

  for(size_t it=0; it<children.size(); it++) {
    StochVector* child = children[it]->newDualYVector();
    y->AddChild(child);
  }
  return y;
}

StochVector* StochTree::newDualZVectorEmpty() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* z = new StochVector(0, commWrkrs);

  for(size_t it=0; it<children.size(); it++) {
    StochVector* child = children[it]->newDualZVector();
    z->AddChild(child);
  }
  return z;
}


StochVector* StochTree::newRhs()
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* rhs = new StochVector(nx() + my() + mz(), commWrkrs);

  for(size_t it=0; it<children.size(); it++) {
    StochVector* child = children[it]->newRhs();
    rhs->AddChild(child);
  }

  return rhs;
}

void StochTree::startMonitors()
{
  iterMon.recIterateTm_start();
  startNodeMonitors();
}

void StochTree::stopMonitors()
{
  iterMon.recIterateTm_stop();
  stopNodeMonitors();
}

void StochTree::startNodeMonitors()
{
  resMon.reset();
  for(size_t i=0; i<children.size(); i++)
    children[i]->startNodeMonitors();
}

void StochTree::stopNodeMonitors()
{
  for(size_t i=0; i<children.size(); i++)
    children[i]->stopNodeMonitors();

  resMon.computeTotal();
}

void StochTree::toMonitorsList(list<NodeExecEntry>& lstExecTm)
{
  lstExecTm.push_back(resMon.eTotal);

  for(size_t i=0; i<children.size(); i++)
    children[i]->toMonitorsList(lstExecTm);
}

void StochTree::fromMonitorsList(list<NodeExecEntry>& lstExecTm)
{
  resMon.eTotal = lstExecTm.front();
  lstExecTm.pop_front();

  for(size_t i=0; i<children.size(); i++)
    children[i]->fromMonitorsList(lstExecTm);
}

void StochTree::syncMonitoringData(vector<double>& vCPUTotal)
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

bool StochTree::balanceLoad()
{
  //before synchronization, compute the total time recorded on this CPU
  //updates this->IPMIterExecTIME
  computeNodeTotal();

  int nCPUs; MPI_Comm_size(MPI_COMM_WORLD, &nCPUs);
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

void StochTree::getSyncInfo(int rank, int& syncNeeded, int& sendOrRecv, int& toFromCPU )
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

int StochTree::isInVector(int elem, const vector<int>& vec)
{
  for(size_t i=0; i<vec.size(); i++)
    if(elem==vec[i]) return 1;
  return 0;
}

void StochTree::computeNodeTotal()
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

void StochTree::saveCurrentCPUState()
{
  myOldMpiComm = commWrkrs;
  myOldProcs = myProcs;

  for(size_t i=0; i<children.size(); i++)
    children[i]->saveCurrentCPUState();
}


#ifdef STOCH_TESTING
void StochTree::runTestNGP()
{

  int noProcs = 397;
  int noNodes = 300;//sizeof(nodes)/sizeof(double);
  vector<double> nodesLoad;
  

  //-----------------------

  /*double nodes[] = { 1.21, 1.205, 1.22, 1.05, 1.02, 
		     1.0, 0.99, 0.9, 1.05, 0.9,
		      0.95, 1.1, 0.88, 0.9, 1.05, 
		      1.1, 1.0, 1.09, 1.01, 1.077};
  */
  double nodes[] = { 10.9101, 11.09205, 11.42, 0.405, 11.02, 
		     1.0, 0.099, 0.9, 1.05, 0.09,
		      21.095, 20.9091, 21.28, 21.09, 20.9995, 
		     1.1, 1.0, 1.09, 1.01, 1.077, 
		     1.21, 1.205, 0.92, 1.05, 1.02, 
		      0.95, 1.1, 0.88, 0.9, 1.05, 
		     0.81, .80, .809, .901, .977,
		     1.28, 0.89, .95,1.021, 1.205, 
		     1.1, 1.0, 1.09, 1.01, 1.077,
		     1.21, 1.205, 1.22, 1.05, 1.02, 
		     1.0, 0.99, 0.9, 1.05, 0.9,
		      1.095, 0.91, 1.18, 1.09, 0.9995, 
		     1.0, 0.99, 0.79, 1.05, 0.9,
		      0.95, 1.1, 0.88, 0.9, 1.05, 
		     1.1, 1.0, 1.09, 1.01, 1.077,
		     1.28, 0.89, .95,1.021, 1.205, 
		     1.1, 1.0, 1.09, 1.01, 1.077, 
		     0.820, 0.8199, 0.83, 1.205, 1.19, //mins
		     1.21, 1.205, 1.22, 1.05, 1.02,
		     1.21, 1.205, 1.22, 1.05, 1.02, 
		     1.0, 0.99, 0.9, 1.05, 0.9,
		      1.095, 0.91, 1.18, 1.09, 0.9995, 
		     1.1, 1.0, 1.09, 1.01, 1.077, 
		     1.21, 1.205, 0.92, 1.05, 1.02, 
		     1.0, 0.99, 0.9, 1.05, 0.9,
		      0.95, 1.1, 0.88, 0.9, 1.05, 
		     0.81, .80, .809, .901, .977,
		     1.28, 0.89, .95,1.021, 1.205, 
		     1.1, 1.0, 1.09, 1.01, 1.077,
		     1.21, 1.205, 1.22, 1.05, 1.02, 
		     1.0, 0.99, 0.9, 1.05, 0.9,
		      1.095, 0.91, 1.18, 1.09, 0.9995, 
		     1.1, 1.0, 1.09, 1.01, 1.077,
		     1.28, 0.89, .95,1.021, 1.205, 
		     1.1, 1.0, 1.09, 1.01, 1.077, 
		     1.0, 0.99, 0.9, 1.05, 0.9,
		     1.21, 1.205, 1.22, 1.05, 1.02, 
		     1.1, 1.178, 0.9, 0.80, 1.1, 
		     1.1, 1.40, 1.09, 1.01, 1.077, 
		     1.0, 0.59, 0.79, 1.05, 0.9, 
		     1.0, 0.99, 0.9, 1.05, 0.9, 
		     0.91, 1.18, 1.09, 0.9995, 1.1, 
		     1.0,1.1, 0.88, 0.9, 1.05,
		     1.0, 0.99, 0.9, 1.05, 0.9,
		      0.95, 1.1, 0.88, 0.9, 1.05, 
		     0.81, .80, .809, .901, .977,
		     1.28, 0.89, .95,1.021, 1.205, 
		     1.21, 1.205, 1.22, 1.05, 1.02, 
		     1.0, 0.99, 0.9, 1.05, 0.9,
		      1.095, 0.91, 1.18, 1.09, 0.9995, 
		     1.1, 1.0, 1.09, 1.01, 1.077, 
		     1.21, 1.205, 1.22, 1.05, 1.02, 
		     1.0, 0.99, 0.79, 1.05, 0.9,
		     1.0, 0.99, 0.9, 1.05, 0.9,
		     1.21, 1.205, 1.22, 1.05, 1.02,
		     1.21, 1.205, 1.22, 1.05, 1.02, 
		     1.0, 0.99, 0.9, 1.05, 0.9,
		      1.095, 0.91, 1.18, 1.09, 0.9995, 
		     1.1, 1.0, 1.09, 1.01, 1.077, 
		     1.21, 1.205, 0.92, 1.05, 1.02, 
		     1.0, 0.99, 0.9, 1.05, 0.9,
		      0.95, 1.1, 0.88, 0.9, 1.05, 
		     0.81, .80, .809, .901, .977,
		     1.28, 0.89, .95,1.021, 1.205, 
		     1.1, 1.0, 1.09, 1.01, 1.077,
		     1.21, 1.205, 1.22, 1.05, 1.02, 
		     1.0, 0.99, 0.9, 1.05, 0.9,
		      1.095, 0.91, 1.18, 1.09, 0.9995, 
		     1.0, 0.99, 0.79, 1.05, 0.9,
		      0.95, 1.1, 0.88, 0.9, 1.05, 
		     1.1, 1.0, 1.09, 1.01, 1.077,
		     1.1, 1.0, 1.09, 1.01, 1.077, 
		     1.0, 0.99, 0.9, 1.05, 0.9,
		     1.21, 1.205, 1.22, 1.05, 1.02, 
		     1.1, 1.178, 0.9, 0.80, 1.1, 
		     1.1, 1.40, 1.09, 1.01, 1.077, 
		     1.0, 0.59, 0.79, 1.05, 0.9, 
		     1.150, 0.99, 0.9, 1.05, 0.9, 
		     0.91, 1.18, 1.09, 0.9995, 1.1, 
		     1.0,1.1, 1.88, .9, 1.05};
  
  nodesLoad.resize(noNodes);for(int it=0; it<noNodes; it++) nodesLoad[it] = nodes[it];

  for(int n=2000; n<7350; n+=250) 
    //int n=1008;
  {
    noProcs=n; double balance; vector<vector<int> > mapping;
    StochResourcePlanner planner;
    planner.setPrintLevel(0);
    planner.assignProcesses(noProcs, nodesLoad, mapping, balance);
    if(balance<1.01 || balance>1.65)
      printf("!!!!!!PROCS[%4d]=%8.4f\n", n, balance);
    if(100*(n/100)==n) 
      printf("PROCS[%4d]=%8.4f\n", n, balance);
  }
}

void StochTree::displayProcessInfo(int rank)
{
  char szTab[5000];
  for(int i=0; i<5000; i++) szTab[i]='\0';
  if(rankMe==rank)
    displayProcessInfo(szTab);
}

void StochTree::displayProcessInfo(char* tab)
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

void StochTree::displayVectorVsTreeStructure(StochVector& stVec, int rank)
{
  char szTab[5000];
  for(int i=0; i<5000; i++) szTab[i]='\0';
 
  if(rankMe==rank)
    displayVectorVsTreeStructure(stVec, rankMe, szTab);
}

void StochTree::displayVectorVsTreeStructure(StochVector& stVec, int myRank, char* tab)
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

void StochTree::displayMatVsTreeStructure(StochGenMatrix& mat, int rank)
{
  char szTab[5000];
  for(int i=0; i<5000; i++) szTab[i]='\0';
 
  if(rankMe==rank)
    displayMatVsTreeStructure(mat, rankMe, szTab);
}

void StochTree::displayMatVsTreeStructure(StochGenMatrix& mat, int myRank, char* tab)
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

void StochTree::displayMatVsTreeStructure(StochSymMatrix& mat, int rank)
{
  char szTab[5000];
  for(int i=0; i<5000; i++) szTab[i]='\0';
 
  if(rankMe==rank)
    displayMatVsTreeStructure(mat, rankMe, szTab);
}

void StochTree::displayMatVsTreeStructure(StochSymMatrix& mat, int myRank, char* tab)
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


void StochTree::displayExecTimes(int rank)
{
  char szTab[5000];
  for(int i=0; i<5000; i++) szTab[i]='\0';

  if(rankMe==rank) displayExecTimes(szTab);
}

void StochTree::displayExecTimes(char* tab)
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

