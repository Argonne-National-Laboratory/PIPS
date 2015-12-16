/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */


#include "StochResourcesMonitor.h"
#include "mpi.h"
#include <assert.h>
#include <fstream>
#include <sstream>

using namespace std;

#define MAX(a,b) ((a > b) ? a : b)

StochNodeResourcesMonitor::StochNodeResourcesMonitor() 
  : eFact(0,0), eLsolve(0,0), eDsolve(0,0), eLtsolve(0,0), eTotal(0,0), eMult(0,0), eReduce(0,0), eReduceScatter(0,0), eBcast(0,0)
{


}

StochNodeResourcesMonitor::~StochNodeResourcesMonitor()
{}

void StochNodeResourcesMonitor::reset()
{
  eFact.clear(); eLsolve.clear();
  eDsolve.clear(); eLtsolve.clear();
  eTotal.clear();
  eMult.clear();
  eReduce.clear();
  eReduceScatter.clear();
  eBcast.clear();
  vcSchur.clear(); vcLsolve.clear();
}


void StochNodeResourcesMonitor::computeTotal()
{
  eTotal.tmLocal = eFact.tmLocal + eLsolve.tmLocal + 
    eDsolve.tmLocal + eLtsolve.tmLocal;

  eTotal.tmChildren = eFact.tmChildren + eLsolve.tmChildren + 
    eDsolve.tmChildren + eLtsolve.tmChildren;

}

void StochNodeResourcesMonitor::recFactTmLocal_start()
{
  tmOpStart = MPI_Wtime();
}

void StochNodeResourcesMonitor::recFactTmLocal_stop()
{
  double tmOpEnd = MPI_Wtime();
  eFact.tmLocal += (tmOpEnd-tmOpStart);
}

void StochNodeResourcesMonitor::recFactTmChildren_start()
{  tmOpStart = MPI_Wtime();  }

void StochNodeResourcesMonitor::recFactTmChildren_stop()
{
  double tmOpEnd = MPI_Wtime();
  eFact.tmChildren += (tmOpEnd-tmOpStart);
}


void StochNodeResourcesMonitor::recLsolveTmLocal_start()
{  tmOpStart = MPI_Wtime();  }

void StochNodeResourcesMonitor::recLsolveTmLocal_stop()
{
  double tmOpEnd = MPI_Wtime();
  eLsolve.tmLocal += (tmOpEnd-tmOpStart);
}
void StochNodeResourcesMonitor::recLsolveTmChildren_start()
{  tmOpStart = MPI_Wtime();  }

void StochNodeResourcesMonitor::recLsolveTmChildren_stop()
{
  double tmOpEnd = MPI_Wtime();
  eLsolve.tmChildren += (tmOpEnd-tmOpStart);
}

void StochNodeResourcesMonitor::recDsolveTmLocal_start()
{  tmOpStart = MPI_Wtime();  }

void StochNodeResourcesMonitor::recDsolveTmLocal_stop()
{
  double tmOpEnd = MPI_Wtime();
  eDsolve.tmLocal += (tmOpEnd-tmOpStart);
}

void StochNodeResourcesMonitor::recDsolveTmChildren_start()
{  tmOpStart = MPI_Wtime();  }

void StochNodeResourcesMonitor::recDsolveTmChildren_stop()
{
  double tmOpEnd = MPI_Wtime();
  eDsolve.tmChildren += (tmOpEnd-tmOpStart);
}

void StochNodeResourcesMonitor::recLtsolveTmLocal_start()
{  tmOpStart = MPI_Wtime();  }

void StochNodeResourcesMonitor::recLtsolveTmLocal_stop()
{
  double tmOpEnd = MPI_Wtime();
  eLtsolve.tmLocal += (tmOpEnd-tmOpStart);
}
void StochNodeResourcesMonitor::recLtsolveTmChildren_start()
{  tmOpStart = MPI_Wtime();  }

void StochNodeResourcesMonitor::recLtsolveTmChildren_stop()
{
  double tmOpEnd = MPI_Wtime();
  eLtsolve.tmChildren += (tmOpEnd-tmOpStart);
}


void StochNodeResourcesMonitor::recSchurCom_start(double size, stCommType type)
{
  NodeCommEntry commMon(0.0, size, type);
  vcSchur.push_back(commMon);
  tmOpStart = MPI_Wtime();
}

void StochNodeResourcesMonitor::recSchurCom_stop(double size, stCommType type)
{
  double tmOpEnd = MPI_Wtime();
  vcSchur[vcSchur.size()-1].time = tmOpEnd-tmOpStart;
}

void StochNodeResourcesMonitor::recLsolveCom_start(double size, stCommType type)
{
  NodeCommEntry commMon(0.0, size, type);
  vcLsolve.push_back(commMon);
  tmOpStart = MPI_Wtime();
}

void StochNodeResourcesMonitor::recLsolveCom_stop(double size, stCommType type)
{
  double tmOpEnd = MPI_Wtime();
  vcLsolve[vcSchur.size()-1].time = tmOpEnd-tmOpStart;
}

//record the dense matrix multiplication within Schur complement computation
void StochNodeResourcesMonitor::recSchurMultChildren_start()
{  eMult.recChildren_start(); }

void StochNodeResourcesMonitor::recSchurMultChildren_stop()
{ eMult.recChildren_stop(); }

void StochNodeResourcesMonitor::recSchurMultLocal_start()
{  eMult.recLocal_start(); }

void StochNodeResourcesMonitor::recSchurMultLocal_stop()
{ eMult.recLocal_stop(); }

void StochNodeResourcesMonitor::recReduceTmLocal_start()
{  eReduce.recLocal_start(); }

void StochNodeResourcesMonitor::recReduceTmLocal_stop()
{ eReduce.recLocal_stop(); }


void StochNodeResourcesMonitor::recReduceScatterTmLocal_start()
{  eReduceScatter.recLocal_start(); }

void StochNodeResourcesMonitor::recReduceScatterTmLocal_stop()
{ eReduceScatter.recLocal_stop(); }

void StochNodeResourcesMonitor::recBcastTmLocal_start()
{ eBcast.recLocal_start(); }

void StochNodeResourcesMonitor::recBcastTmLocal_stop()
{ eBcast.recLocal_stop(); }

//**************************************************
//*************  NodeExecEntry   *******************
//**************************************************
NodeExecEntry::NodeExecEntry(double localTime, double childrenTime) 
    : tmLocal(localTime), tmChildren(childrenTime)
{};

void NodeExecEntry::clear() 
{ tmLocal=0.0; tmChildren=0.0;}

void NodeExecEntry::recLocal_start()
{ tmOpStartLocal = MPI_Wtime(); }

void NodeExecEntry::recLocal_stop()
{
  double tmOpEnd = MPI_Wtime();
  tmLocal += (tmOpEnd-tmOpStartLocal);
};

void NodeExecEntry::recChildren_start()
{ tmOpStartChildren = MPI_Wtime(); }

void NodeExecEntry::recChildren_stop()
{
  double tmOpEnd = MPI_Wtime();
  tmChildren += (tmOpEnd-tmOpStartChildren);
};

//**************************************************
//*************  ITERATE monitor *******************
//**************************************************
StochIterateResourcesMonitor::StochIterateResourcesMonitor() {}
StochIterateResourcesMonitor::~StochIterateResourcesMonitor() {}

void
StochIterateResourcesMonitor::recIterateTm_start()
{
  tmIterateStart=MPI_Wtime();
}

void 
StochIterateResourcesMonitor::recIterateTm_stop()
{
  tmIterate = MPI_Wtime()-tmIterateStart;
}


// returns the "high water mark" of resident memory usage
// during the execution of this process
unsigned long StochIterateResourcesMonitor::getMaxMemUsage()
{
  ifstream fd("/proc/self/status");
  string line;
  const string numbers("0123456789");
  size_t pos;
  unsigned long hwm, rss, *toset;
  while (!fd.eof()) {
    toset = 0;
    getline(fd, line);
    string numstr;
    pos = line.find("VmHWM");
    if (pos != string::npos) {
	    toset = &hwm;
    } 
    pos = line.find("VmRSS");
    if (pos != string::npos) {
	    toset = &rss;
    }
    if (toset) {
	    numstr = line.substr(
			    line.find_first_of(numbers),line.find_last_of(numbers)+1);
	    istringstream s(numstr);
	    s >> *toset;
    }
  }
  fd.close();
  
  // according to documentation, rss can be larger than hwm
  return MAX(hwm, rss);
  
}
