/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */


#ifndef STOCH_RESOURCE_PLANNER
#define STOCH_RESOURCE_PLANNER

#include <vector>
#include <stdio.h>

using namespace std;

class StochResourcePlanner {
 public:
  StochResourcePlanner();
  virtual ~StochResourcePlanner();

  void assignProcesses(int noProcesses, 
		       vector<double>& nodesLoad,
		       vector<vector<int> >& mapping);

  void assignProcesses(int noProcesses, 
		       vector<double>& nodesLoad,
		       vector<vector<int> >& mapping,
		       double& balance);

  void assignProcesses(vector<int>& ranks,
		       vector<double>& nodesLoad,
		       vector<vector<int> >& mapping,
		       double& balance);

  void assignProcesses(vector<int>& ranks,
		       vector<double>& nodesLoad,
		       vector<vector<int> >& mapping);
  
  void setPrintLevel(int level);
  

  static int noScaProcesses;
  

 protected:
  int printLevel;
};


class LoadBalancing
{
 public:
  LoadBalancing(vector<double>& vNodesLoad, const int noProcs);
  LoadBalancing(vector<double>& vNodesLoad, vector<int>& vRanks);
  virtual ~LoadBalancing();

  /**
   * The nodes 0,1,...,N-1  having the load specified in the vector
   * vNodeLoad (see constructor) of size N will be assigned to P processes specified by
   * vRanks[0], vRanks[1],...,vRanks[P-1]. A node may be spanned to
   * multiple processes. Each CPUs is assigned to at least one node.
   * 
   * Returns a vector v with N elements. Element v[i] contains a
   * vector of size at least one storing the ranks assigned to node i.
   */
  void assign(std::vector<std::vector<int> >& v, double& balancing, int& allAssigned);

  /**
   * If for the node i, we have  
   *    1-l1low < cpuLoad(i) <= 1+lupp, then assign 1 process to node i
   * otherwise
   *    if n-1+lupp < cpuLoad(i) <=n+lupp, then assign n processes to  i
   *    otherwise assign based on number partitioning algorithm.
   *
   * Here: cpuLoad(i) is the exec time for node i.
   */
  void setPredefinedAssignmentLimits(double l1low, double lupp);
  
  /**
   * Different assignment strategies are used.
   *  - S1 N less than P 
   *  - S2 N close to P
   *  - S3 N greater than P. 
   * 
   * If N/P <= rLess, then S1 is used.
   * If rLess < N/P <= rGreater, then S2 is used.
   * If N/p > rGreater, then S3 is used.
   *
   * Default values: 0.9 and 1.1.
   */
  void setCompareCriteria( double rLess, double rGreater);

  /**
   * Computes and returns the balance bal of the assignment.
   *    bal = max{ load of  process i : i=0,1,...,P} / ( totalLoad/P )
   */ 
  void getBalance(double* balance, double* metisBalance=NULL);


  /**
   * Partitioning into a large number of partitions does NOT always
   * work. In this case, METIS may produce empty subparts. 
   * If processes needs to be assigned to a large number of nodes (or
   * viceversa), then both the nodes and processes are partitioned in 
   * s sets, 2<=s<=maxParts and s smaller partitioning problems are
   * solved. The actual value of s depends also on the value of
   * minSizeSubpart, see @setMinSizeSubpart.
   *
   * s is computed as the largest number less than or equal to
   * maxParts for which each of the s subpartitions have (approx.)
   * more elements than minSizeSubpart. This is enforced to avoid bad
   * balancing caused by assigning 3 processes to 2 nodes (in this
   * case an assignment of 9 procs to 6 nodes give better balance).
   *
   * Default value is 20.
   */
  void setMaximumNoParts(int maxParts);

  /**
   * This function sets  the 'minSizeSubpart' parameter use in
   * computing the number s introduced above.
   *
   * Default value 10.
   */
  void setMinSizeSubpart(int minSize);

  /**
   * When N and P have close values METIS may return empty partitions,
   * i.e., some of the CPUs do not get work(=nodes).
   * Call this function with flag=1 if all CPUs have to be assigned. A
   * simple algorithm is used, nodes from CPUs having larger load are
   * assigned to CPUs having no load.
   *
   * Using this option may not give any computational speed-up and increases
   * the inter-CPU communication. The benefit would be a decrease of
   * the memory requirements in the CPUs which the nodes are taken
   * from.
   *
   * Default value: 1
   */
  void setForceFullAssignment(int flag);

  void setPrintLevel(int level);
 
 protected:
  void assignNLP(std::vector<std::vector<int> >& v, double& balancing, int& allAssigned);
  void assignNEP(std::vector<std::vector<int> >& v, double& balancing, int& allAssigned);
  void assignNGP(std::vector<std::vector<int> >& v, double& balancing, int& allAssigned);

  /** 
   * Assigns nodes to CPUs for the case N>>P.
   */
  void assignNodesToProcs(std::vector<int>& nodes,
			  std::vector<double>& nodesLoad, 
			  std::vector<int>& ranks, 
			  std::vector<std::vector<int> >& mapping);

  /** 
   * Assigns nodes to CPUs for the case N<<P
   */
  void assignProcsToNodes(std::vector<int>& nodes,
			  std::vector<double>& nodesLoad, 
			  std::vector<int>& ranks, 
			  std::vector<std::vector<int> >& mapping);

  void lastSweepNLP(std::vector<int>& nodes,
		    std::vector<double>& nodesLoad, 
		    std::vector<int>& ranks, 
		    std::vector<std::vector<int> >& mapping);
  
  void normalizeVec(std::vector<double>& vec);

  void preAssignNodesToProcs(const std::vector<int>& nodes,
			     const std::vector<double>& nodesLoad, 
			     const std::vector<int>& ranks,
			     std::vector<int>& availNodes, 
			     std::vector<double>& availNodesLoad,
			     std::vector<int>& availRanks,
			     std::vector<int>& preAssNodes, 
			     std::vector<int>& preAssRanks,
			     std::vector<std::vector<int> >& gMapping);

  void preAssignProcsToNodes(const std::vector<int>& nodes,
			     const std::vector<double>& nodesLoad, 
			     const std::vector<int>& ranks,
			     std::vector<int>& availNodes, 
			     std::vector<double>& availNodesLoad,
			     std::vector<int>& availRanks,
			     std::vector<int>& preAssNodes, 
			     std::vector<int>& preAssRanks,
			     std::vector<std::vector<int> >& gMapping);

  int decideNoSubPartitions(int noProcs);

  /**
   * Wrapper around metis_GraphPartKway
   *
   * Solves the number partitioning problem for 'nodes' by generating
   * 'noSubParts' partitions.
   * Nodes have equal weights.
   */
  void partition(const std::vector<int>& nodes, int noSubParts, std::vector<std::vector<int> >& partitions);

  void partitionVert(const std::vector<int>& nodes, const std::vector<double>& nodesLoad, 
			   int nparts, std::vector<std::vector<int> >& partitions);

  void partitionWeightVert(const std::vector<int>& nodes, const std::vector<double>& nodesLoad, 
			   const std::vector<double>& partWeights,
			   int nparts, std::vector<std::vector<int> >& partitions);

  void computeBalance(const std::vector<std::vector<int> >& mapPart, 
		      int noSubParts,
		      double& balance, int& allAssigned);
  
  void computeBalanceW(const std::vector<std::vector<int> >& mapPart, 
		       const std::vector<double>& nodeWeights, 
		       int noSubParts, 
		       double& balance, int& allAssigned);

  void computeBalanceWW(const std::vector<std::vector<int> >& mapPart, 
		       const std::vector<double>& nodeWeights, 
		       const std::vector<double>& partWeights,
		       int noSubParts, 
		       double& balance, int& allAssigned);

  void faPartition(const std::vector<int>& nodes, const std::vector<double>& nodesLoad, 
		   int nparts, std::vector<std::vector<int> >& partitions);

  void faPartition2(const std::vector<int>& nodes, const std::vector<double>& nodesLoad, 
		   const std::vector<int>& ranks, std::vector<std::vector<int> >& partitions);

  void faPartition3(const std::vector<int>& nodes, const std::vector<double>& nodesLoad, 
		    const std::vector<int>& ranks, 
		    std::vector<std::vector<int> >& partitions);

  void faPartition4(const std::vector<int>& nodes, const std::vector<double>& nodesLoad, 
		    const std::vector<int>& ranks, 
		    std::vector<std::vector<int> >& partitions);

  void computePartitionWeights(const std::vector<std::vector<int> >& partition, int nparts,
			       std::vector<double>& weights);

  std::vector<int> extractIndxForPart(int ipart, const std::vector<int>& nodes, 
				      const std::vector<std::vector<int> >& parts);
  std::vector<int> extractNodesForPart(int ipart, const std::vector<int>& nodes, 
				       const std::vector<std::vector<int> >& parts);
  std::vector<double> extractElems(const std::vector<int>& nodes,
				   const std::vector<double>& nodesLoad, 
				   const std::vector<int>& nodesSubSet);

  void updateMapping_N2P_N2P (std::vector<std::vector<int> >& parts, const std::vector<int>& nodes,
			      const std::vector<std::vector<int> >& subParts, 
			      const std::vector<int>& subNodes, const std::vector<int>& subRanks);

  void updateMapping_N2P_PIdx(std::vector<std::vector<int> >& parts, 
			      const std::vector<int>& nodes, const std::vector<int>& ranks, 
			      const std::vector<std::vector<int> >& subParts);

  void updateMapping_N2P_NIdx(std::vector<std::vector<int> >& parts, 
			      const std::vector<int>& nodes, const std::vector<int>& ranks, 
			      const std::vector<std::vector<int> >& subParts);
  
 protected:
  double preAssLim1Low, preAssLimUpp;
  double ratLess, ratGreater;
  int maxNoParts;
  int maxMaxsToConnect, maxMinsToConnect;
  int forceFullAssignment;
  int minElemsInSubpart;

  vector<double>& nodesLoadW;
  vector<int> ranksW;

  int printLevel;

};

/* log functions */
//template <class V>
//void printVector(const V& v, const char* msginit, const char* fmt, const char* msgfin="\n");
//template <class M>
//void printMap(const M& m, const char* msginit, const char* fmt, const char* msgfin="\n");

template <class V>
void printVector(const V& v, const char* msginit, const char* fmt, const char* msgfin="\n")
{
  printf("%s", msginit);
  for(int i=0; i<v.size(); i++) printf(fmt, i, v[i]);
  printf("%s", msgfin);

}

template <class M>
void printMap(const M& m, const char* msginit, const char* fmt, const char* msgfin="\n")
{
  printf("%s", msginit);
  for(int i=0; i<m.size(); i++) {
    for(int j=0; j<m[i].size(); j++)
      printf(fmt, i, j, m[i][j]);
    printf("\n");
  }
  printf("%s", msgfin);
}

#endif
