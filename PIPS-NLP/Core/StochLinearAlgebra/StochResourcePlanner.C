/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "StochResourcePlanner.h"

#include "stdio.h"
#include "assert.h"
#include <math.h>
//#include <stdlib.h>
#include <algorithm>
/* 
 * #include "metis.h" follows below
 * 
 */
extern "C" {
#include "defs.h"
#include "struct.h"
};

using namespace std;

extern "C" 
void metis_wpartgraphkway(int *, idxtype *, idxtype *, idxtype *, idxtype *, 
			  int *, int *, int *, float *, int *, int *, idxtype *);
extern "C" 
void metis_wpartgraphrecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *, 
			       int *, int *, int *, float *, int *, int *, idxtype *);


int StochResourcePlanner::noScaProcesses;
 
StochResourcePlanner::StochResourcePlanner()
  :printLevel(0)
{

}

StochResourcePlanner::~StochResourcePlanner() {};

void 
StochResourcePlanner::assignProcesses(int noProcs, 
				      vector<double>& nodesLoad,
				      vector<vector<int> >& mapping)
{
  double balance;
  vector<int> ranks(noProcs);
  for(int p=0;p<noProcs; p++) ranks[p]=p;

  this->assignProcesses(ranks, nodesLoad, mapping, balance);
}

void 
StochResourcePlanner::assignProcesses(int noProcs, 
				      vector<double>& nodesLoad,
				      vector<vector<int> >& mapping,
				      double& balance)
{
  vector<int> ranks(noProcs);
  for(int p=0;p<noProcs; p++) ranks[p]=p;

  this->assignProcesses(ranks, nodesLoad, mapping, balance);
}

void 
StochResourcePlanner::assignProcesses(vector<int>& vRanks, 
				      vector<double>& nodesLoad,
				      vector<vector<int> >& mapping)
{
  double balance;
  this->assignProcesses(vRanks, nodesLoad, mapping, balance);
}

void 
StochResourcePlanner::assignProcesses(vector<int>& vRanks, 
				      vector<double>& nodesLoad,
				      vector<vector<int> >& mapping,
				      double& balance)
{
  LoadBalancing lb(nodesLoad, vRanks);
  lb.setPrintLevel(printLevel);

  int noProcs = vRanks.size();

  int optMaxNoParts = -1; int optMinElemsInSubpart=-1; double minBalance=1.0e10;
  int allAssigned;

  for(int MNoP=45; MNoP>=20; MNoP-=5) {
    lb.setMaximumNoParts(MNoP);

    for(int mES=2; mES<19; mES+=3) {
      mapping.clear();
      allAssigned=0;
      lb.setMinSizeSubpart(mES);
      lb.assign(mapping, balance, allAssigned);

      if(allAssigned) {
	if(minBalance > balance) { 
	  minBalance=balance;
	  optMaxNoParts=MNoP;
	  optMinElemsInSubpart=mES;
	}
      }	else {
	//if(printLevel>0) 
	printf("ResPlanner: warning: allAssigned=0 for MNoP=%d mES=%d\n", MNoP, mES);
      }
    }
  }
  if(optMaxNoParts<0 || optMinElemsInSubpart<0) {
    printf("ResPlanner: error: allAssigned=0 \n");
    assert(false);
  }

  if(printLevel>9) {
    printf("ResPlanner: optimal partition params found: (%d,%d)\n",
	   optMaxNoParts, optMinElemsInSubpart);
  }

  lb.setMinSizeSubpart(optMinElemsInSubpart);lb.setMaximumNoParts(optMaxNoParts);
  mapping.clear();
  lb.assign(mapping, balance, allAssigned);
 

  if(printLevel>0) {
    printf("ResPlanner: assigned %d CPUs for %lu NODEs (%g,%d)\n",
	   noProcs, nodesLoad.size(), balance, allAssigned);
  }
}

  
void StochResourcePlanner::setPrintLevel(int level)
{
  printLevel=level;
}

//===========================================================================================
//===========================================================================================
//======================   METIS based CPU SCHEDULER  =======================================
//===========================================================================================
//===========================================================================================


/* predicates and structs used in sorting */
struct IdxLoad{
  int idx; double load;
};
bool IdxLoadIncrPred(const IdxLoad& il1, const IdxLoad& il2)
{ return il1.load<il2.load; }
bool IdxLoadIncrPredIdx(const IdxLoad& il1, const IdxLoad& il2)
{ return il1.idx<il2.idx; }

/* LoadBalancing's members */

LoadBalancing::LoadBalancing(vector<double>& vNodesLoad, const int noProcs)
  : nodesLoadW(vNodesLoad)
{
  preAssLim1Low = 0.0125;
  preAssLimUpp  = 0.8;
  ratLess       = 0.8;
  ratGreater    = 1.1;

  maxNoParts=30; minElemsInSubpart=15;

  maxMaxsToConnect = 15;
  maxMinsToConnect = 15;

  forceFullAssignment = 1;

  ranksW.resize(noProcs);
  for(int i=0; i<noProcs; i++) ranksW[i]=i;
}

LoadBalancing::LoadBalancing(vector<double>& vNodesLoad, vector<int>& vRanks)
  : nodesLoadW(vNodesLoad)
{
  preAssLim1Low = 0.0125;
  preAssLimUpp  = 0.8;
  ratLess       = 0.8;
  ratGreater    = 1.1;

  maxNoParts=30; minElemsInSubpart=15;

  maxMaxsToConnect = 15;
  maxMinsToConnect = 15;

  forceFullAssignment = 1;

  int noProcs = vRanks.size();
  ranksW.resize(noProcs);
  for(int i=0; i<noProcs; i++) ranksW[i]=vRanks[i];;
}

LoadBalancing::~LoadBalancing() {}

void LoadBalancing::setForceFullAssignment(int flag)
{
  forceFullAssignment = flag;
}

void LoadBalancing::setMaximumNoParts(int maxParts)
{
  maxNoParts = maxParts;
}

void LoadBalancing::setMinSizeSubpart(int minSize)
{
  minElemsInSubpart = minSize;
}

void LoadBalancing::setCompareCriteria( double rLess, double rGreater)
{
  ratLess=rLess; ratGreater=rGreater;
}

void LoadBalancing::setPrintLevel(int level)
{
  printLevel=level;
}

void LoadBalancing::assign(vector<vector<int> >& mapping, double& balancing, int& allAssigned)
{
  int noNodes = nodesLoadW.size();
  int noProcs = ranksW.size();

  double ratio = 1.0*noNodes/noProcs;
  if(ratio<ratLess) assignNLP(mapping, balancing, allAssigned);
  else if(ratio<ratGreater) assignNEP(mapping, balancing, allAssigned);
  else assignNGP(mapping, balancing,  allAssigned);
}
void LoadBalancing::assignNEP(vector<vector<int> >& mapping, double& balancing, int& allAssigned)
{
  int noNodes = nodesLoadW.size();
  int noProcs = ranksW.size();

  if(printLevel>9) {
    printf("============================================\n");
    printf("  NEP Assign %d NODES to %lu CPUs (%d,%d)\n", 
	   noNodes, ranksW.size(), maxNoParts, minElemsInSubpart);
    printf("===========================================\n");
  }
  double balanceNLP,balanceNGP; int allAssignedNLP,allAssignedNGP;
  vector<vector<int> > mappingNLP(noNodes), mappingNGP(noNodes);

  if(noNodes<=noProcs*1.25) {
    assignNLP(mappingNLP, balanceNLP, allAssignedNLP);
  } else {
    allAssignedNLP=0;
    balanceNLP=1.0e10;
  }

  if(noNodes>=noProcs*0.85) {
    assignNGP(mappingNGP, balanceNGP, allAssignedNGP);
  } else {
    allAssignedNGP=0;
    balanceNGP=1.0e10;
  }

  if(printLevel>9)
    printf("NGP(%g,%d) NLP(%g,%d)\n", balanceNGP, allAssignedNGP, balanceNLP, allAssignedNLP);

  mapping     = balanceNLP<balanceNGP ? mappingNLP:mappingNGP;
  balancing   = balanceNLP<balanceNGP ? balanceNLP:balanceNGP;
  allAssigned = balanceNLP<balanceNGP ? allAssignedNLP:allAssignedNGP;
}

void LoadBalancing::assignNLP(vector<vector<int> >& mapping, double& balance, int& allAssigned)
{
  int noNodes = nodesLoadW.size();
  mapping.resize(noNodes);

  if(printLevel>9) {
    printf("===========================================\n");
    printf("  NLP Assign %d NODES to %lu CPUs (%d,%d)\n", 
	   noNodes, ranksW.size(), maxNoParts, minElemsInSubpart);
    printf("===========================================\n");
  }
  vector<int> nodes(noNodes);
  for(int i=0; i<noNodes; i++) nodes[i]=i;

  //percentages of total load
  this->normalizeVec(nodesLoadW);

  //----------------------------------------------------------------
  //------------------pre-assignment--------------------------------
  //----------------------------------------------------------------
  vector<int> preAssRanks, availRanks;
  vector<int> preAssNodes, availNodes;
  vector<vector<int> >  preAssMapping;
  vector<double> availNodesLoad;
  
  preAssignProcsToNodes(nodes, nodesLoadW, ranksW, 
	    availNodes, availNodesLoad, availRanks, 
	    preAssNodes, preAssRanks, mapping);

  if(printLevel>9) 
    printf("NLP: %lu NODES were pre-assigned to %lu CPUS\n", 
	   preAssNodes.size(), preAssRanks.size());

  if(availNodes.size()>availRanks.size()) {
    if(printLevel>9) printf("NLP NOT APLIC. Aborting! %lu NODES %lu CPUs\n",
			    availNodes.size(),availRanks.size());
    balance=1e10;
    allAssigned=0;
    return;
  }

  assignProcsToNodes(availNodes, availNodesLoad, availRanks, mapping);

  computeBalanceW(mapping, nodesLoadW, ranksW.size(), balance, allAssigned);
  if(printLevel>99) printf("before sweep balance=%g allassigned=%d\n", balance,  allAssigned);
  //lastSweepNLP(availNodes, availNodesLoad, availRanks, mapping);
  lastSweepNLP(nodes, nodesLoadW, ranksW, mapping);

  computeBalanceW(mapping, nodesLoadW, ranksW.size(), balance, allAssigned);
  
  if(printLevel>9) printf("GLOBAL-NLP balance=%g allassigned=%d\n", balance,  allAssigned);
}

void LoadBalancing::assignNGP(vector<vector<int> >& mapping, double& balance, int& allAssigned)
{
  int noNodes = nodesLoadW.size();
  mapping.resize(noNodes);
  
  if(printLevel>9) {
    printf("============================================\n");
    printf("  NGP Assign %d NODES to %lu CPUs (%d,%d)\n", 
	   noNodes, ranksW.size(), maxNoParts, minElemsInSubpart);
    printf("===========================================\n");
  }
  vector<int> nodes(noNodes);
  for(int i=0; i<noNodes; i++) nodes[i]=i;

  //percentages of total load
  this->normalizeVec(nodesLoadW);


  //----------------------------------------------------------------
  //------------------pre-assignment--------------------------------
  //----------------------------------------------------------------
  vector<int> preAssRanks, availRanks;
  vector<int> preAssNodes, availNodes;
  vector<vector<int> >  preAssMapping;
  vector<double> availNodesLoad;

  preAssignNodesToProcs(nodes, nodesLoadW, ranksW, 
	    availNodes, availNodesLoad, availRanks, 
	    preAssNodes, preAssRanks, mapping);

  if(printLevel>9) 
    printf("NGP: %lu NODES were pre-assigned to %lu CPUS\n", 
	   preAssNodes.size(), preAssRanks.size());
	
  if(availNodes.size()<availRanks.size()) {
    if(printLevel>9) printf("NGP NOT APLIC. Aborting! %lu NODES %lu CPUs\n",
			    availNodes.size(),availRanks.size());
    balance=1e10;
    allAssigned=0;
    return;
  }

		  
  //----------------------------------------------------------------
  //---------------------assignment---------------------------------
  //----------------------------------------------------------------
  //vector<vector<int> > assMapping(availNodes.size());
  assignNodesToProcs(availNodes, availNodesLoad, availRanks, mapping);

  computeBalanceW(mapping, nodesLoadW, ranksW.size(), balance, allAssigned);
  if(printLevel>9) printf("GLOBAL-NGP balance=%g allassigned=%d\n", balance, allAssigned);

}


void 
LoadBalancing::assignProcsToNodes(vector<int>& nodes,
			      vector<double>& nodesLoad, 
			      vector<int>& ranks, 
			      vector<vector<int> >& globMappingN2P)
{
  int noNodes = nodesLoad.size(); int noProcs = ranks.size();

  int noSubParts = decideNoSubPartitions(noNodes);
  if(noSubParts>1) {
    //split the CPUs in noSubParts partitions of equal weight and
    //calculate the CPU weights corresponding to each partition
    vector<vector<int> > procsPartition, nodesPartition;
    partition(ranks, noSubParts, procsPartition);


    vector<double> procPartWeights;
    computePartitionWeights(procsPartition, noSubParts, procPartWeights);

    // ----------------------------------------------------------------
    // Partition the nodes accordingly to the CPU weighted partitions.
    // ----------------------------------------------------------------
    partitionWeightVert(nodes, nodesLoad, procPartWeights, noSubParts, nodesPartition);

    //!log
    double procPartBalance; int allProcsAssigned;
    computeBalance(procsPartition, noSubParts, procPartBalance, allProcsAssigned);

    double nodePartBalance; int allNodesAssigned;
    computeBalanceWW(nodesPartition, nodesLoad, 
		    procPartWeights, 
		    noSubParts, nodePartBalance, allNodesAssigned);
    if(printLevel>9) {
      printf("----------------------------------------------------------------------------------\n");
      printf("SUBPARTITION NLP: %lu NODES (w%5.3f,%d) %d CPUs(%5.3f,%d) in %d partitions.\n",
	     nodes.size(), nodePartBalance, allNodesAssigned, noProcs,
	     procPartBalance, allProcsAssigned, noSubParts);
      printf("----------------------------------------------------------------------------------\n");
    }

    //check if each NODES partition has fewer elements than CPUS
    //partition; if no, then assign them manually.
    int NLP = 1;
    for(int ipart=0; ipart<noSubParts && NLP; ipart++) {

      vector<int> nodes1 = extractNodesForPart(ipart, nodes, nodesPartition);
      vector<int> ranks1 = extractNodesForPart(ipart, ranks, procsPartition);
      NLP = (nodes1.size()<=ranks1.size() && nodes1.size()>0 && ranks1.size()>0);
    }
    if(NLP) {
      //partition each subpartition
      for(int ipart=0; ipart<noSubParts; ipart++) {
	
	vector<int> nodes1 = extractNodesForPart(ipart, nodes, nodesPartition);
	vector<int> ranks1 = extractNodesForPart(ipart, ranks, procsPartition);
	
	vector<double> nodesLoad1 = extractElems(nodes, nodesLoad, nodes1);
	
	assignProcsToNodes(nodes1, nodesLoad1, ranks1, globMappingN2P);
	
	if(printLevel>999) printf("~~~END NLP subpart %d\n", ipart);
      }
    } else {
      faPartition4(nodes, nodesLoad, ranks, globMappingN2P);
    }
  } else { //noSubparts==1

    //need a temp buffer that containts the weights of the 
    //CPUs, i.e., 1/N
    vector<double> cpuWeights(ranks.size(), 1.0/ranks.size());
    normalizeVec(nodesLoad);

    vector<vector<int> > mapP2NIdx;
    partitionWeightVert(ranks, cpuWeights, nodesLoad, nodesLoad.size(), mapP2NIdx); 

    double balance; int allAssigned;
    computeBalanceWW(mapP2NIdx,
		     cpuWeights, nodesLoad, 
		     nodesLoad.size(),balance, allAssigned);

    if(printLevel>99) printf("Partitioned NLP %lu CPUS into %lu parts (ww%g,%d) \n", 
			     ranks.size(), nodesLoad.size(),balance, allAssigned); 
    
    if(0==allAssigned) {

      if(printLevel>99) printf("Some NODES were NOT assigned -> enforcing full NLP assignment ....\n");   
   
      faPartition3(nodes, nodesLoad, ranks, mapP2NIdx);
      
      //printMap(mapN2PIdx, "parts", "[%d][%d]=%d");
      computeBalanceWW(mapP2NIdx,
		       cpuWeights, nodesLoad, 
		       nodesLoad.size(),balance, allAssigned);
      
     if(printLevel>999) printf("faPartition3 AFTER  %lu NODES to %lu CPUS (%6.3f,%d)\n", 
			     nodes.size(), ranks.size(), balance, allAssigned);
    }
    updateMapping_N2P_NIdx(globMappingN2P, ranks, nodes, mapP2NIdx);


  }
}

void 
LoadBalancing::lastSweepNLP(vector<int>& nodes,
			vector<double>& nodesLoad, 
			vector<int>& ranks, 
			vector<vector<int> >& gMapping)
{
  int nMaxOne=-1; int nMinMany=-1; int howMany=0;
  double maxOne=0.0; double minMany=1.0e10;
  int level=1;

  int maxCPU=0;

  do{
    nMaxOne=-1; nMinMany=-1; howMany=0; maxOne=0.0; minMany=1.0e10;

    for(size_t n=0; n<nodes.size(); n++) {
      if(gMapping[n].size()-level==0) {
	if(maxOne<nodesLoad[n]) {maxOne = nodesLoad[n]; nMaxOne=n;}
      } else { //gMapping[n].size>=1

	if(gMapping[n].size()-level==1)
	  if(minMany>nodesLoad[n]) {
	    minMany=nodesLoad[n]; 
	    nMinMany=n; howMany=gMapping[n].size();
	  }

	if(maxCPU-gMapping[n].size()<0) maxCPU=gMapping[n].size();
      }
    }
    //do we have to switch?
    if(maxOne>minMany && nMinMany>=0) {
      vector<int> aux = gMapping[nMaxOne];
      gMapping[nMaxOne] = gMapping[nMinMany];
      gMapping[nMinMany] = aux;

      if(printLevel>999)
	printf("sweeping level=%2d: n(%3d,%6.4f) %2lu CPU <---> n(%3d,%6.4f) %2d CPUs\n",
	       level, nMaxOne, maxOne, aux.size(), nMinMany, minMany, howMany);

    } else {level++;}
  } while(level<=maxCPU); 
}

void 
LoadBalancing::assignNodesToProcs(vector<int>& nodes,
			      vector<double>& nodesLoad, 
			      vector<int>& ranks, 
			      vector<vector<int> >& globMappingN2P)
{
  int noProcs = ranks.size();

  int noSubParts = decideNoSubPartitions(noProcs);
  if(noSubParts>1) {

    //split the CPUs in noSubParts partitions of equal weight and
    //calculate the CPU weights corresponding to each partition
    vector<vector<int> > procsPartition, nodesPartition;
    partition(ranks, noSubParts, procsPartition);


    vector<double> procPartWeights;
    computePartitionWeights(procsPartition, noSubParts, procPartWeights);

    // ---------------------
    // Partition the nodes
    // ---------------------
    partitionWeightVert(nodes, nodesLoad, procPartWeights, noSubParts, nodesPartition);


    //!log
    double procPartBalance; int allProcsAssigned;
    computeBalance(procsPartition, noSubParts, procPartBalance, allProcsAssigned);

    double nodePartBalance; int allNodesAssigned;
    computeBalanceWW(nodesPartition, nodesLoad, 
		    procPartWeights, 
		    noSubParts, nodePartBalance, allNodesAssigned);
    if(printLevel>9) {
      printf("----------------------------------------------------------------------------------\n");
      printf("SUBPARTITION NGP: %lu NODES (w%5.3f,%d) %d CPUs(%5.3f,%d) in %d partitions.\n",
	     nodes.size(), nodePartBalance, allNodesAssigned, noProcs,
	     procPartBalance, allProcsAssigned, noSubParts);
      printf("----------------------------------------------------------------------------------\n");
    }

    //check if each NODES partition has more elements than CPUS
    //partition; if no, then assign them manually.
    int NGP = 1;
    for(int ipart=0; ipart<noSubParts && NGP; ipart++) {

      vector<int> nodes1 = extractNodesForPart(ipart, nodes, nodesPartition);
      vector<int> ranks1 = extractNodesForPart(ipart, ranks, procsPartition);
      NGP = (nodes1.size()>=ranks1.size() && nodes1.size()>0 && ranks1.size()>0);

    }
    if(NGP) {
      for(int ipart=0; ipart<noSubParts; ipart++) {
	
	vector<int> nodes1 = extractNodesForPart(ipart, nodes, nodesPartition);
	vector<int> ranks1 = extractNodesForPart(ipart, ranks, procsPartition);
	
	vector<double> nodesLoad1 = extractElems(nodes, nodesLoad, nodes1);
	
	assignNodesToProcs(nodes1, nodesLoad1, ranks1, globMappingN2P);
	
	if(printLevel>999) printf("~~~END NGP subpart %d\n", ipart);
      }
    } else {
      //if(nodes.size()>=ranks.size())
	faPartition2(nodes, nodesLoad, ranks, globMappingN2P);
	// else
	//assignProcsToNodes(nodes, nodesLoad, ranks, globMappingN2P);
    }

  } else { //here only one part needed, number of CPUs is small enough
    //partitionV(nodes, nodesLoad, noProcs, nodesPartition);

    if(nodes.size()>ranks.size()) {

      vector<vector<int> > mapN2PIdx;
      partitionVert(nodes, nodesLoad, ranks.size(), mapN2PIdx);
      
      double balance; int allAssigned;
      computeBalance(mapN2PIdx, ranks.size(), balance, allAssigned);
      if(0==allAssigned && forceFullAssignment==1) {
	//!log
	if(printLevel>99) {
	  printf("Some CPUs node NOT assigned -> forcing full assignment ....\n");
	  printf("faPartition BEFORE %lu NODES to %lu CPUS (%6.3f,%d)\n", 
		 nodes.size(), ranks.size(), balance, allAssigned);
	}	
	
	faPartition(nodes, nodesLoad, ranks.size(), mapN2PIdx);
	

	computeBalance(mapN2PIdx, ranks.size(), balance, allAssigned);
	if(printLevel>99) printf("faPartition AFTER  %lu NODES to %lu CPUS (%6.3f,%d)\n", 
				nodes.size(), ranks.size(), balance, allAssigned);
      }
      updateMapping_N2P_PIdx(globMappingN2P, nodes, ranks, mapN2PIdx);
    }else {   
      vector<vector<int> > mapP2NIdx;
      assignProcsToNodes(nodes, nodesLoad, ranks, globMappingN2P);
      //updateMapping_N2P_NIdx(globMappingN2P, ranks, nodes,
      //mapP2NIdx);
      //printMap(globMappingN2P, "in P2N Glob  map\n", "[%d][%d]=%d\n");
    }
    
  }
}
  


void LoadBalancing::preAssignNodesToProcs(const vector<int>& nodes,
			  const vector<double>& nodesLoad, 
			  const vector<int>& ranks,
			  vector<int>& availNodes, 
			  vector<double>& availNodesLoad,
			  vector<int>& availRanks,
			  vector<int>& preAssNodes, 
			  vector<int>& preAssRanks,
			  vector<vector<int> >& gMapping)
{
 
  int noNodes = nodes.size();
  int noProcs = ranks.size();

  preAssRanks.clear(); preAssNodes.clear();
  availNodes.clear();  availRanks.clear();

  int assignedCPU=noProcs-1;

  //normalize the loads and multiply each by the number of avail CPUs 
  //use this swap to identify and process nodes that need 1 or more CPUs
  for(int itn=0; itn<noNodes; itn++) {
    
    double cpus = nodesLoad[itn]*noProcs;

    if(cpus>1-preAssLim1Low) {
      
      if(cpus<1+preAssLimUpp) {
	//node itn gets 1 cpu only for himself

	vector<int> vCPUs(1); vCPUs[0] = ranks[assignedCPU]; 
	gMapping[itn] = vCPUs; //preAssMapping.push_back( vCPUs );
	preAssNodes.push_back(nodes[itn]);
	preAssRanks.push_back(ranks[assignedCPU]);
	assignedCPU--;

	//printf("PreAssign: node %5d  gets  1 cpu  out of  %g. [%d] assigned.\n", itn, cpus, ranks[assignedCPU+1]);

      }else {
	int ncpu = (int)floor(cpus);
	if ( cpus-ncpu > preAssLimUpp) 
	  ncpu++;

	//make sure this node does not get all the cpus so that the
	//subunitary nodes have NO node to run on.
	if(ncpu==assignedCPU+1) ncpu--;

	vector<int> vCPUs(ncpu); 
	for(int itp=0; itp<ncpu; itp++) {
	  vCPUs[itp] = ranks[assignedCPU]; 
	  preAssRanks.push_back(ranks[assignedCPU]);
	  assignedCPU--;
	}
	gMapping[itn] = vCPUs; //preAssMapping.push_back( vCPUs );
	preAssNodes.push_back(nodes[itn]);
	//printf("PreAssign: node %5d  gets  %d cpus out of  %g. [%d][%d] assigned.\n", 
	//, ncpu, cpus, ranks[assignedCPU+1], ranks[assignedCPU+2]);
      }
    } else {
      //subunitary
      availNodes.push_back(nodes[itn]);
    }
  }
  for(int i=0; i<=assignedCPU; i++) availRanks.push_back(ranks[i]);

  assert(availRanks.size()+preAssRanks.size() == ranks.size());
  assert(availNodes.size()+preAssNodes.size() == nodes.size());
 
  availNodesLoad.clear();
  availNodesLoad.resize(availNodes.size());
  for(size_t i=0; i<availNodes.size(); i++) availNodesLoad[i] = nodesLoad[availNodes[i]];
}


void LoadBalancing::preAssignProcsToNodes(const vector<int>& nodes,
			     const vector<double>& nodesLoad, 
			     const vector<int>& ranks,
			     vector<int>& availNodes, 
			     vector<double>& availNodesLoad,
			     vector<int>& availRanks,
			     vector<int>& preAssNodes, 
			     vector<int>& preAssRanks,
			     vector<vector<int> >& gMapping)
{

  //identify those nodes who should get less than 1 CPU and partition
  //them using NGP strategy
  int noProcs = ranks.size();
  int noNodes = nodes.size(); assert(nodesLoad.size()-noNodes==0);
  
  //identify those nodes who should get less than 1 CPU and partition
  //this using NGP strategy
  vector<double> preAssNodesLoad; double preAssCPUTotal=0.0;
  for(size_t i=0; i<nodes.size(); i++) {
    if(noProcs*nodesLoad[i]<1.0) {
      preAssNodes.push_back(nodes[i]);
      preAssNodesLoad.push_back(nodesLoad[i]);
      preAssCPUTotal += (noProcs*nodesLoad[i]);
    } else {
      availNodes.push_back(nodes[i]);
      availNodesLoad.push_back(nodesLoad[i]);
    }
  }

  int preAssNoProcs = (int)ceil(preAssCPUTotal); assert(preAssNoProcs<noProcs);
  for(size_t i=0; i<ranks.size(); i++)
    if(i-preAssNoProcs<0) preAssRanks.push_back(ranks[i]);
    else availRanks.push_back(ranks[i]);

  if(preAssNoProcs>0) {
    if(preAssNodes.size() >= preAssRanks.size()) {
      if(printLevel>9) 
	printf("preAssignNLP: NGP a total of %lu nodes will be NGPed to %d CPUs.\n",
	       preAssNodes.size(), preAssNoProcs);
      
      assignNodesToProcs(preAssNodes, preAssNodesLoad, preAssRanks, gMapping);
    } else {
      
      printf("preAssignNLP: warning:  a total of %lu nodes will be NLPed to %d CPUs.\n",
	     preAssNodes.size(), preAssNoProcs);

      assignProcsToNodes(preAssNodes, preAssNodesLoad, preAssRanks, gMapping);
    }
  }
}

void LoadBalancing::normalizeVec(vector<double>& vec)
{
  double total=0.0;
  for(size_t i=0; i<vec.size(); i++) total += vec[i];
  for(size_t i=0; i<vec.size(); i++) vec[i] /= total;;
}

void LoadBalancing::computeBalance(const vector<vector<int> >& mapPart, int noSubParts, 
			       double& balance, int& allAssigned)
{
  allAssigned=1;   balance = 0.0;

  double total = 1.0*mapPart.size();

  double* parts = new double[noSubParts];
  for(int i=0; i<noSubParts; i++) parts[i] = 0.0;

  for(size_t i=0; i<mapPart.size(); i++) {

    int CPUsPerThisNode = mapPart[i].size();
    
    assert(0!=CPUsPerThisNode);

    if(1==CPUsPerThisNode) parts[mapPart[i][0]] += 1.0;
    else {

      for(int j=0; j<CPUsPerThisNode; j++)
	parts[mapPart[i][j]] += 1.0/CPUsPerThisNode;
      assert(false);

    }
  }

  for(int i=0; i<noSubParts; i++) {
    //parts[i] /= total;
    if(balance<parts[i]/(total/noSubParts)) balance = parts[i]/(total/noSubParts);

    if(parts[i]==0.0) {allAssigned = 0;}
  }
  
  delete[] parts;
}

void LoadBalancing::computeBalanceW(const vector<vector<int> >& mapPart, const vector<double>& nodeWeights,
				int noSubParts, 
				double& balance, int& allAssigned)
{
  allAssigned=1;   balance = 0.0;

  assert(mapPart.size() == nodeWeights.size());

  double total = 0.0;
  for(size_t i=0; i<nodeWeights.size(); i++) total += nodeWeights[i];
  

  double* parts = new double[noSubParts];
  for(int i=0; i<noSubParts; i++) {parts[i] = 0.0;}

  for(size_t i=0; i<mapPart.size(); i++) {

    int CPUsPerThisNode = mapPart[i].size();
    if(0==CPUsPerThisNode) {
      if(printLevel>0) 
	printf("[%lu] has nothing assigned. subparts=%d\n", i, noSubParts);
      assert(0!=CPUsPerThisNode);
    } 

    if(1==CPUsPerThisNode) {
      parts[mapPart[i][0]] += nodeWeights[i];
      //printf("Node[%3d] w=%7.4f has [%2d] cpus:%3d\n", i, nodeWeights[i], CPUsPerThisNode, mapPart[i][0]);
      

    } else {     
      //printf("Node[%3d] w=%7.4f has [%2d] cpus:", i, nodeWeights[i], CPUsPerThisNode);
      for(int j=0; j<CPUsPerThisNode; j++) {
	parts[mapPart[i][j]] += nodeWeights[i]/CPUsPerThisNode;
	//printf("%3d ", mapPart[i][j]);
      }
      //printf("\n");
    }
    
  }

  for(int i=0; i<noSubParts; i++) {

    if(balance < parts[i]/(total/noSubParts)) {
      balance = parts[i]/(total/noSubParts);
      //printf("CPU[%d]-> balance[%g]\n", i, balance);
    }

    if(parts[i]==0.0) {allAssigned = 0;}
  }
  delete[] parts;
}

void 
LoadBalancing::computeBalanceWW(const vector<vector<int> >& mapPart, 
			   const vector<double>& nodeWeights, const vector<double>& partWeights,
			   int noSubParts, 
			   double& balance, int& allAssigned)
{
  allAssigned=1;  balance = 0.0;

  assert(mapPart.size()    == nodeWeights.size());
  assert(partWeights.size() - noSubParts == 0);

  double total = 0.0;
  for(size_t i=0; i<mapPart.size(); i++) total += nodeWeights[i];

  double* parts = new double[noSubParts];
  for(int i=0; i<noSubParts; i++) parts[i] = 0.0;

  for(size_t i=0; i<mapPart.size(); i++) {
    assert(mapPart[i].size()==1);
    parts[mapPart[i][0]] += (nodeWeights[i] / total);

    //printf("Proc[%d] holds node[%d] of weight %g\n", i,  mapPart[i][0], partWeights[mapPart[i][0]]);
  }
  for(int i=0; i<noSubParts; i++) {
    double received = parts[i];
    if(received==0) {balance=1.0e10; allAssigned=0; break;}

    double askedFor = partWeights[i];
    double ratio = askedFor / received;

    if(balance<ratio) {balance = ratio;}//printf("balance became %g because of the node[%d] \n", ratio,i);}
  }

  //balance /= (total/noSubParts);;
  delete[] parts;
}

vector<int> 
LoadBalancing::extractNodesForPart(int ipart, const vector<int>& nodes, const vector<vector<int> >& parts)
{
  vector<int> ret;

  assert(nodes.size() == parts.size());

  for(size_t i=0; i<parts.size(); i++) {
    
    assert(parts[i].size()==1);
    for(size_t p=0; p<parts[i].size(); p++) {
      if(parts[i][p] == ipart)
	ret.push_back(nodes[i]);
    }
  }
  //assert(ret.size()>=1);

  return ret;
}

vector<int> 
LoadBalancing::extractIndxForPart(int ipart, const vector<int>& nodes, const vector<vector<int> >& parts)
{
  size_t i,p;
  vector<int> ret;

  assert(nodes.size() == parts.size());

  for(i=0; i<parts.size(); i++) {
    
    assert(parts[i].size()==1);
    for(p=0; p<parts[i].size(); p++) {
      if(parts[i][p] == ipart)
	ret.push_back(i);
    }
  }
  assert(ret.size()>=1);

  return ret;
}

vector<double> 
LoadBalancing::extractElems(const vector<int>& nodes,
			const vector<double>& nodesLoad, 
			const vector<int>& nodesSubSet)
{
  vector<double> ret(nodesSubSet.size(), 0.0);;
  int noNodes = nodesLoad.size();
  assert(noNodes - nodesSubSet.size() >= 0);
  assert(noNodes - nodes.size() == 0);

  
  for(size_t j=0; j<nodesSubSet.size(); j++)
    for(int i=0; i<noNodes; i++) {
      if(nodes[i]==nodesSubSet[j])
	ret[j] = nodesLoad[i];
  }
#ifdef DEBUG
  for(size_t j=0; j<nodesSubSet.size(); j++)
    assert(ret[j]==0);
#endif
  return ret;
}


int findElemInVec(int el, const vector<int>& vec) 
{
  int i=0; int n=vec.size();
  while(i<n) 
    if(el==vec[i]) return i;
    else i++;
  return -1;
}

// -----------------------------------------------
// parts and subParts are maps from NODES to RANKS
// -----------------------------------------------
void 
LoadBalancing::updateMapping_N2P_N2P(vector<vector<int> >& parts, const vector<int>& nodes,
				 const vector<vector<int> >& subParts, 
				 const vector<int>& subNodes, const vector<int>& subRanks)
{
  size_t isp;
  assert(subParts.size()==subNodes.size());
  assert(parts.size()==nodes.size());

  for(isp=0; isp<subParts.size(); isp++) {
    int ip = findElemInVec(subNodes[isp], nodes); assert(ip==-1);
    
    assert(parts[ip].size()==0); //this node should not have been
				 //assigned before
    parts[ip] = subParts[isp];
  } 
}

// -----------------------------------------------
// parts    - a map from NODES to RANKS
// subParts - a map from NODES to pos in RANKS
// -----------------------------------------------
void 
LoadBalancing::updateMapping_N2P_PIdx(vector<vector<int> >& parts, 
				  const vector<int>& nodes, const vector<int>& ranks,
				  const vector<vector<int> >& subParts)
{
  size_t isp;
  assert(subParts.size()==nodes.size());

  for(isp=0; isp<nodes.size(); isp++) {
    int ip = nodes[isp];

    //!log 
    //if(parts[ip].size()>0) {
    //  printf("size(parts[%d])=%d\n",ip, parts[ip].size());
    //  for(size_t j=0; j<parts[ip].size(); j++) printf("%d\n", parts[ip][j]);
    //}

    assert(parts[ip].size()==0); //this node should not have been
				 //assigned before
    assert(subParts[isp].size()>0); //all the nodes were assigned
				  //processes to 

    parts[nodes[isp]].resize(subParts[isp].size());

    for(size_t r=0; r<subParts[isp].size(); r++) {
      parts[nodes[isp]][r] = ranks[ subParts[isp][r] ];
    }
  } 
}

// -----------------------------------------------
// parts    - a map from NODES to RANKS
// subParts - a map from RANKS to pos in NODES
// -----------------------------------------------
void 
LoadBalancing::updateMapping_N2P_NIdx(vector<vector<int> >& parts, 
				  const vector<int>& ranks, const vector<int>& nodes,
				  const vector<vector<int> >& subParts)
{
  assert(subParts.size()==ranks.size());
  //printMap(subParts, "parts\n", "[%d][%d]=%d", "\n\n");

  for(size_t p=0; p<ranks.size(); p++) {
    
    assert(0!=subParts[p].size());

    for(size_t n=0; n<subParts[p].size(); n++) {
      parts[nodes[subParts[p][n]]].push_back( ranks[p] );
      //parts[nodes[subParts[p][n]]].push_back( subParts[p][n] );
    }
  }
}

int LoadBalancing::decideNoSubPartitions(int noProcs)
{

  //if we have no more than 20(=maxNoParts) CPUs, we will do no SUBpartitioning
  if(noProcs<=maxNoParts) return 1;

  //!
  return 1;

  int noParts = maxNoParts; 

  //if we have 65 nodes, we are not happy because aprox 65/20=3 CPUs will be
  //in each partition and this causes a bad balance. Reduce the the
  //number of partitions to 12 so we will have 65/12=5 cpus per
  //partition, i.e., better balance.

  assert(maxNoParts>=10);          //for debug :)
  if(maxNoParts<10) maxNoParts=10; //for release

  //at least xxx elements per partition
  while(noProcs/noParts<minElemsInSubpart) {
    noParts--;
  }
  return noParts;
}

// nodes have same weight
void LoadBalancing::computePartitionWeights(const vector<vector<int> >& partition, int nparts,
					vector<double>& weights)
{
  int i;

  int noNodes = partition.size();

  weights.clear(); weights.resize(nparts);
  for(i=0; i<nparts; i++) weights[i] = 0.0;

  for(size_t i=0; i<partition.size(); i++) {
    //for(int c=0; c<partition[i].size(); c++)
    assert(partition[i].size() == 1);
    weights[partition[i][0]] += (1 / (1.0*noNodes)); 
  }
}

void LoadBalancing::partition(const vector<int>& nodes, int nparts, vector<vector<int> >& partitions)
{
  int noNodes = nodes.size();

  if(nparts>1) {
  //------------------------------------------------------------------
  // Graph setup 
  //------------------------------------------------------------------
  int noEdges = 2*(noNodes);

  idxtype* xadj = new idxtype[noNodes+1];
  idxtype *adjncy = new idxtype[noEdges]; //we may overallocate by a litle
  

  //build the xadj and adjncy arrays
  xadj[0] = 0;
  
  for(int inode=0; inode<noNodes; inode++) {
    int edges = xadj[inode];
    //add edges for the previous and next one
    adjncy[edges] = inode==0?noNodes-1:inode-1;
    edges++;
    assert(edges<=noEdges);

    adjncy[edges] = inode==noNodes-1?0:inode+1;
    edges++;
    assert(edges<=noEdges);

    xadj[inode+1] = edges;
  }

  assert(xadj[noNodes]==noEdges);

  //nodes weight
  idxtype *vwgt = NULL;
  //edges weights
  idxtype* adjwgt = NULL;

  //flags & options
  int wgtflag = 0; int numflag=0; int options[5];
  options[0] = 0;

  //parts
  int npartits = nparts;

  float* tpwgts = new float[nparts];
  for (int i=0; i<nparts; i++) 
    tpwgts[i] = 1.0/(1.0*(nparts));

  //output
  int edgecut;
  idxtype* part = new idxtype[noNodes];


  //printf("calling metis...."); fflush(stdout);
  metis_wpartgraphkway(&noNodes, 
			     //metis_wpartgraphkway(&noNodes, 
			    xadj, 
			    adjncy, 
			    vwgt,      //info about the weights of the vertices
			    adjwgt,    //info about the weights of the edges
			    &wgtflag,
			    &numflag,
			    &npartits,
			    tpwgts,
			    options, 
			    &edgecut, 
			    part);
  //printf(" metis finished.\n"); fflush(stdout);

  delete[] tpwgts;
  delete[] vwgt;
  delete[] adjncy;
  delete[] xadj;

  partitions.clear(); partitions.resize(noNodes);

  for(int itn=0; itn<noNodes; itn++) {
    vector<int> vCPUs(1); vCPUs[0] = part[itn];    
    partitions[itn] = vCPUs;
  }

  delete[] part;
  } else {
    assert(nparts==1);
    partitions.clear(); partitions.resize(noNodes);
    for(int itn=0; itn<noNodes; itn++) {
      vector<int> vCPUs(1); vCPUs[0] = 0;
      partitions[itn] = vCPUs;
    }
  }
}

void LoadBalancing::partitionVert(const vector<int>& nodes, const vector<double>& nodesLoad, 
			   int nparts, vector<vector<int> >& partitions)
{
  int noNodes = nodesLoad.size();
  //take a look into the load vector and find the largest noProcs and
  //smallest noProcs loads and connect them in graph
  if(nparts>1) {
  int noMaxs = nparts/2 < maxMaxsToConnect ? nparts/2 : maxMaxsToConnect;
  int noMins = nparts/2 < maxMinsToConnect ? nparts/2 : maxMinsToConnect;

  //printf("partitionWV: assigning %d nodes to %d partitions\n", noNodes, nparts);
  //printf("%d Maxs will be connected to %d mins\n", noMaxs, noMins);

  vector<IdxLoad> idM;
  vector<IdxLoad> idm;

  if(noMaxs>0 && noMins>0) {
    idM.resize(noMaxs);
    idm.resize(noMins);

    assert(noMaxs<noNodes); 
    assert(noMins<noNodes);

    for(int itn=0; itn<noMaxs; itn++) { idM[itn].idx=itn; idM[itn].load = nodesLoad[itn]; }
    for(int itn=0; itn<noMins; itn++) { idm[itn].idx=itn; idm[itn].load = nodesLoad[itn]; }

    sort(idM.begin(), idM.end(), IdxLoadIncrPred);
    sort(idm.begin(), idm.end(), IdxLoadIncrPred);
    

    for(int itn=min(noMaxs,noMins); itn<noNodes; itn++) {
      if(itn>=noMaxs)
	if(nodesLoad[itn]>idM[0].load) {
	  idM[0].idx  = itn;
	  idM[0].load = nodesLoad[itn];
	  sort(idM.begin(), idM.end(), IdxLoadIncrPred);
	}
      if(itn>=noMins)
	if(nodesLoad[itn]<idm[noMins-1].load) {
	  idm[noMins-1].idx  = itn;
	  idm[noMins-1].load = nodesLoad[itn];
	  sort(idm.begin(), idm.end(), IdxLoadIncrPred);
	}
      
    }
    assert(idM.size()-noMaxs==0);
    assert(idm.size()-noMins==0);
  }

  //! check to not be bogus
  sort(idM.begin(), idM.end(), IdxLoadIncrPredIdx);
  sort(idm.begin(), idm.end(), IdxLoadIncrPredIdx);

  //------------------------------------------------------------------
  // Graph setup 
  //------------------------------------------------------------------
  int noEdges = 2*(noNodes+noMaxs*noMins);

  int ncon = 1;
  idxtype* xadj = new idxtype[noNodes+1];
  idxtype *adjncy = new idxtype[noEdges]; //we may overallocate by a
					  //litle in case the mins are
					  //neighbours of  the maxs
  
  //build the xadj and adjncy arrays
  xadj[0] = 0;

  int iMin=0; int iMax=0;
  
  for(int inode=0; inode<noNodes; inode++) {
    int edges = xadj[inode];
    //add edges for the previous and next one
    adjncy[edges] = inode==0?noNodes-1:inode-1;
    edges++;
    assert(edges<=noEdges);
    adjncy[edges] = inode==noNodes-1?0:inode+1;
    edges++;
    assert(edges<=noEdges);

    if(iMax-idM.size()<0) {
      if(inode == idM[iMax].idx) {
	//one of the maxs -> add edges to the mins
	for(size_t j=0; j<idm.size(); j++)
	  if( abs(idm[j].idx-inode)!=1 && idm[j].idx!=inode ) {
	    adjncy[edges] = idm[j].idx; edges++;
	  }
	iMax++;
      }
    }
    assert(edges<=noEdges);
    if(iMin-idm.size()<0) {
      if(inode==idm[iMin].idx) {
	//is one of the mins -> add edges to maxs
	for(size_t j=0; j<idM.size(); j++)
	  if( abs(idM[j].idx-inode)!=1 && idM[j].idx!=inode ) {  
	    //the max and min not the same and also they are not
	    //neighbours (i.e. avoid adding the same edge twice)
	    adjncy[edges] = idM[j].idx; edges++;
	  }
	iMin++;
      }
    }
    assert(edges<=noEdges);
    xadj[inode+1] = edges;
  }

  assert(xadj[noNodes]<=noEdges);
  
  //node weights
  idxtype *vwgt = new idxtype[ncon*noNodes];
  assert(ncon==1);
  for(int i=0; i<noNodes; i++) vwgt[i] = floor(10000*nodesLoad[i]);

  //edges weights
  idxtype* adjwgt = NULL;

  //flags & options
  int wgtflag = 2; int numflag=0; int options[5];
  options[0] = 0;

  float* tpwgts = new float[nparts];
  for (int i=0; i<nparts; i++) 
    tpwgts[i] = 1.0/nparts;

  //output
  int edgecut;
  idxtype* part = new idxtype[noNodes];


  //printf("calling metis...."); fflush(stdout);
  metis_wpartgraphkway(&noNodes, 
			     //metis_wpartgraphkway(&noNodes, 
		       xadj, 
		       adjncy, 
		       vwgt,      //info about the weights of the vertices
		       adjwgt,    //info about the weights of the edges
		       &wgtflag,
		       &numflag,
		       &nparts, 
		       tpwgts,
		       options, 
		       &edgecut, 
		       part);
  //printf(" metis finished.\n"); fflush(stdout);

  delete[] tpwgts;
  delete[] vwgt;
  delete[] adjncy;
  delete[] xadj;

  partitions.clear(); partitions.resize(noNodes);

  vector<double> partLoad(nparts, 0.0);
  for(int itn=0; itn<noNodes; itn++) {
    vector<int> vCPUs(1); vCPUs[0] =  part[itn] ;
    
    partitions[itn] = vCPUs;
  }

  delete[] part;
  } else {
    assert(nparts==1);
    partitions.clear(); partitions.resize(noNodes);
    for(int itn=0; itn<noNodes; itn++) {
      vector<int> vCPUs(1); vCPUs[0] = 0;
      partitions[itn] = vCPUs;
    }
  }
}


void LoadBalancing::partitionWeightVert(const vector<int>& nodes, const vector<double>& nodesLoad, 
			   const vector<double>& partWeights,
			   int nparts,
			   vector<vector<int> >& partitions)
{
  int noNodes = nodesLoad.size();
  //take a look into the load vector and find the largest noProcs and
  //smallest noProcs loads and connect them in graph
  if(nparts>1) {
  int noMaxs = nparts/2 < maxMaxsToConnect ? nparts/2 : maxMaxsToConnect;
  int noMins = nparts/2 < maxMinsToConnect ? nparts/2 : maxMinsToConnect;

  //printf("partitionWV: assigning %d nodes to %d partitions\n", noNodes, nparts);
  // printf("%d Maxs will be connected to %d mins\n", noMaxs, noMins);

  assert(nparts-partWeights.size()==0);

  vector<IdxLoad> idM;
  vector<IdxLoad> idm;

  if(noMaxs>0 && noMins>0) {
    idM.resize(noMaxs);
    idm.resize(noMins);

    assert(noMaxs<noNodes); 
    assert(noMins<noNodes);

    for(int itn=0; itn<noMaxs; itn++) { idM[itn].idx=itn; idM[itn].load = nodesLoad[itn]; }
    for(int itn=0; itn<noMins; itn++) { idm[itn].idx=itn; idm[itn].load = nodesLoad[itn]; }

    sort(idM.begin(), idM.end(), IdxLoadIncrPred);
    sort(idm.begin(), idm.end(), IdxLoadIncrPred);
    

    for(int itn=min(noMaxs,noMins); itn<noNodes; itn++) {
      if(itn>=noMaxs)
	if(nodesLoad[itn]>idM[0].load) {
	  idM[0].idx  = itn;
	  idM[0].load = nodesLoad[itn];
	  sort(idM.begin(), idM.end(), IdxLoadIncrPred);
	}
      if(itn>=noMins)
	if(nodesLoad[itn]<idm[noMins-1].load) {
	  idm[noMins-1].idx  = itn;
	  idm[noMins-1].load = nodesLoad[itn];
	  sort(idm.begin(), idm.end(), IdxLoadIncrPred);
	}
      
    }
    assert(idM.size()-noMaxs==0);
    assert(idm.size()-noMins==0);
  }

  //! check to not be bogus
  sort(idM.begin(), idM.end(), IdxLoadIncrPredIdx);
  sort(idm.begin(), idm.end(), IdxLoadIncrPredIdx);

  //printf("MIN INFO:"); for(int i=0;i<idm.size();i++) printf("[%d %g]", idm[i].idx, idm[i].load); printf("\n");
  //printf("MAX INFO:"); for(int i=0;i<idM.size();i++) printf("[%d %g]", idM[i].idx, idM[i].load); printf("\n");

  //------------------------------------------------------------------
  // Graph setup 
  //------------------------------------------------------------------
  int noEdges = 2*(noNodes+noMaxs*noMins);

  int ncon = 1;
  idxtype* xadj = new idxtype[noNodes+1];
  idxtype *adjncy = new idxtype[noEdges]; //we may overallocate by a
					  //litle in case the mins are
					  //neighbours of  the maxs
  //build the xadj and adjncy arrays
  xadj[0] = 0;

  int iMin=0; int iMax=0;
  
  for(int inode=0; inode<noNodes; inode++) {
    int edges = xadj[inode];
    //add edges for the previous and next one
    adjncy[edges] = inode==0?noNodes-1:inode-1;
    edges++;
    assert(edges<=noEdges);
    adjncy[edges] = inode==noNodes-1?0:inode+1;
    edges++;
    assert(edges<=noEdges);

    if(iMax-idM.size()<0) {
      if(inode == idM[iMax].idx) {
	//one of the maxs -> add edges to the mins
	for(size_t j=0; j<idm.size(); j++)
	  if( abs(idm[j].idx-inode)!=1 && idm[j].idx!=inode ) {
	    adjncy[edges] = idm[j].idx; edges++;
	  }
	iMax++;
      }
    }
    assert(edges<=noEdges);
    if(iMin-idm.size()<0) {
      if(inode-idm[iMin].idx==0) {
	//is one of the mins -> add edges to maxs
	for(size_t j=0; j<idM.size(); j++)
	  if( abs(idM[j].idx-inode)!=1 && idM[j].idx!=inode ) {  
	    //the max and min not the same and also they are not
	    //neighbours (i.e. avoid adding the same edge twice)
	    adjncy[edges] = idM[j].idx; edges++;
	  }
	iMin++;
      }
    }
    assert(edges<=noEdges);
    xadj[inode+1] = edges;
  }

  assert(xadj[noNodes]<=noEdges);
  
  //node weights
  idxtype *vwgt = new idxtype[ncon*noNodes];
  assert(ncon==1);
  for(int i=0; i<noNodes; i++) vwgt[i] = floor(10000*nodesLoad[i]);

  //edges weights
  idxtype* adjwgt = NULL;

  //flags & options
  int wgtflag = 2; int numflag=0; int options[5];
  options[0] = 0;

  float* tpwgts = new float[nparts];
  for (int i=0; i<nparts; i++) 
    tpwgts[i] = partWeights[i];

  //output
  int edgecut;
  idxtype* part = new idxtype[noNodes];

  metis_wpartgraphkway(&noNodes, 
		       //metis_wpartgraphrecursive(&noNodes, 
		       xadj, 
		       adjncy, 
		       vwgt,      //info about the weights of the vertices
		       adjwgt,    //info about the weights of the edges
		       &wgtflag,
		       &numflag,
		       &nparts, 
		       tpwgts,
		       options, 
		       &edgecut, 
		       part);

  delete[] tpwgts;
  delete[] vwgt;
  delete[] adjncy;
  delete[] xadj;

  partitions.clear(); partitions.resize(noNodes);

  vector<double> partLoad(nparts, 0.0);

  for(int itn=0; itn<noNodes; itn++) {
    vector<int> vCPUs(1); vCPUs[0] = part[itn];
    
    partitions[itn] = vCPUs;
  }

  delete[] part;
  } else {
    //assert(nparts==1);
    if (nparts == 1) {
      printf("warning: nparts == 1. ignore if not using metis for assigning scenarios");
    }
    partitions.clear(); partitions.resize(noNodes);
    for(int itn=0; itn<noNodes; itn++) {
      vector<int> vCPUs(1); vCPUs[0] = 0;
      partitions[itn] = vCPUs;
    }
  }

}


struct Node{
  int n;
  double load;
};

bool NodeIncrPred(const Node& n1, const Node& n2)
{
  return n1.load<n2.load;
}

bool NodeDecrPred(const Node& n1, const Node& n2)
{
  return n1.load>n2.load;
}


/*
 * This function is called in the NGP case when we try to subpartition the nodes and
 * CPUs in the same number of subpartitions but we get N<P for one of the subpartitions.
 *
 * The calling code sends to this function the larger sets of nodes
 * and ranks, not the subpartition for which N<P.
 *
 * The first P largest nodes are assigned to each CPU, and then the
 * remaining N-P nodes are assigned, starting with the smallest
 * assigned to the rank containing the smallest of the first P largest  nodes.
 */
void LoadBalancing::faPartition2(const vector<int>& nodes, const vector<double>& nodesLoad, 
			     const vector<int>& ranks, vector<vector<int> >& partitions)
{
  assert(nodes.size()>=ranks.size());
  assert(nodes.size()==nodesLoad.size());

  vector<Node> vNodesS(nodes.size());
  for(size_t i=0;i<nodes.size(); i++){vNodesS[i].n=nodes[i]; vNodesS[i].load=nodesLoad[i]*ranks.size();}
  
  sort(vNodesS.begin(), vNodesS.end(), NodeDecrPred);

  //!log for(size_t i=0;i<nodes.size(); i++) printf("v[%d]=%g\n", vNodesS[i].n, vNodesS[i].load);

  size_t n=0;
  for(; n<ranks.size(); n++)
    partitions[ vNodesS[n].n ].push_back(ranks[n]);

  assert(n==ranks.size());
  
  int r=ranks.size()-1;

  while(n < nodes.size() ) {
    partitions[vNodesS[n].n].push_back(ranks[r]);
    n++; r--;
    if(r<0) r=ranks.size()-1;
  }

}

/* 
 * This function is called in the case NEP or NLP when METIS fails to
 * assigned all the nodes to CPUs.
 *
 * partitions = mapping from CPUS to nodes indexes
 */
void LoadBalancing::faPartition3(const vector<int>& nodes, const vector<double>& nodesLoad, 
			    const vector<int>& ranks, vector<vector<int> >& partitions)
{
  assert(nodes.size()<=ranks.size());
  assert(partitions.size()==ranks.size());
  int ncpus = partitions.size();
  int noNodes = nodes.size();

  int* noCPUs = new int[noNodes];

  while(true) {
    for(int n=0; n<noNodes; n++) noCPUs[n]=0;

    //find the number of cpus each node is assigned to.
    for(int r=0; r<ncpus; r++) 
      for(size_t n=0; n<partitions[r].size(); n++)
	noCPUs[ partitions[r][n]  ] ++;

    //Find the node n0 not present on any CPU
    int n0=-1;
    for(int n=0; n<noNodes && n0<0; n++) if(noCPUs[n]==0) n0 = n;
    if(n0==-1) break; //done here

    //printf("found node %d not assigned\n", n0);

    //Find the CPU rMax with max load sharing nodes with another CPU
    int rMax=-1; double maxLoad=0.0;
    for(int r=0; r<ncpus; r++) {
      double load=0.0; int isSharingNodes=0;
      int N = partitions[r].size(); assert(N>0);
      
      for(int n=0; n<N; n++) {
	load += nodesLoad[partitions[r][n]]/noCPUs[partitions[r][n]];
	if(noCPUs[partitions[r][n]]>1) {isSharingNodes=1;}
      }

      if(maxLoad<load && isSharingNodes) { maxLoad=load; rMax=r; }
    }
    assert(rMax>=0); //should have this here

    //replace node(s) from CPU rMin with node n0. 
    partitions[rMax].clear(); partitions[rMax].push_back(n0);
  }
  delete[] noCPUs;
}

/*
 * This function is called in the NLP case when we try to subpartition the nodes and
 * CPUs in the same number of subpartitions but we get N>P for one of the subpartitions.
 *
 * The calling code sends to this function the larger sets of nodes
 * and ranks, not the subpartition for which N>P.
 *
 * 
 *
 */
void LoadBalancing::faPartition4(const vector<int>& nodes, const vector<double>& nodesLoad, 
			    const vector<int>& ranks, vector<vector<int> >& partitions)
{
  int noNodes = nodes.size(); int noProcs = ranks.size();
  assert(noNodes<=noProcs);

  vector<IdxLoad> vNodes(noNodes);
  for(int i=0; i<noNodes; i++) {
    vNodes[i].idx  = nodes[i]; 
    vNodes[i].load = nodesLoad[i];
  }
  sort(vNodes.begin(), vNodes.end(), IdxLoadIncrPred);

  int n;
  for(n=0; n<noNodes; n++) {
    partitions[vNodes[n].idx].push_back(ranks[n]);
  }

  n=noNodes-1;
  for(int c=noNodes; c<noProcs; c++) {
    partitions[vNodes[n].idx].push_back(ranks[c]);
    n--;
    if(n<0) n=noNodes-1;
  }

}
/*
 * This function is called in the case NGP when no subpartitioning can
 * be done (small # of elems) and METIS fails to assign nodes to all CPUs.
 *
 */
void LoadBalancing::faPartition(const vector<int>& nodes, const vector<double>& nodesLoad, 
			    int nparts, vector<vector<int> >& partitions)
{
  assert(nodes.size()-nparts>=0);
  assert(partitions.size()==nodes.size());

  int*    parts = new int[nparts];
  double* loads = new double[nparts];

  int allAssigned=0;
  while(!allAssigned) {
    for(int p=0; p<nparts; p++) {parts[p]=0; loads[p]=0.0;}

    //populate parts with number of nodes assigned to each part
    //populate loads with total load assigned to part  
    
    for(size_t n=0; n<nodes.size(); n++) {
      for(size_t p=0; p<partitions[n].size(); p++) {
	parts[partitions[n][p]]++;
	loads[partitions[n][p]] += nodesLoad[n]/partitions.size();
      } 
    }
  

    //which part has no node assigned and which part has the max load
    //and at least two nodes
    int pMax=-1; int p0 =-1; double maxLoad=0.0;
    for(int p=0; p<nparts; p++) {

      if(parts[p]==0) p0=p;

      if(maxLoad < loads[p] && parts[p] > 1) {
	pMax=p;
	maxLoad = loads[p];
      } 
    }

    //assign the "largest" node in the "largest" partition to the
    //empty partition
    
    if(p0<0) {allAssigned=1;}
    else {
      assert(pMax>=0);
 
      //one node from partition pMax is sent to partition p0
      vector<int> nodesPmax = extractIndxForPart(pMax, nodes, partitions);
      assert(nodesPmax.size()>0);

      //get the maximum load node from this partition
      int nMax=nodesPmax[0];
      for(size_t i=1; i<nodesPmax.size(); i++)
	if(nodesLoad[nodesPmax[i]]>nodesLoad[nMax]) nMax = nodesPmax[i];

      //put nMax in the partition p0
      assert(partitions[nMax].size()==1);
      partitions[nMax][0] = p0;
    }

  } //end while(!allAssigned)
 
  delete[] parts;
  delete[] loads;
  
}


