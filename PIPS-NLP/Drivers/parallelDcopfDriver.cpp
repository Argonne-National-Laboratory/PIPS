/* PIPS-NLP                                                           	*
 * Author:  Nai-Yuan Chiang                                       	*
 * (C) 2015 Argonne National Laboratory. 				*/

#include "dcopflowInput.hpp"
#include "graphPart.h"
#include "ps.h"

#include "NlpPIPSIpmInterface.h"

#include "sFactoryAug.h"
#include "sFactoryAugAggregationPrecond.h"

#include "FilterIPMStochSolver.h"

#include "sNlpInfoFIX.h"

#include "pipsOptions.h"



using namespace std;
extern int gOuterSolve;
extern int gInnerSCsolve;
extern int gisNLP;

int solve               (const string& dataname, int nscen);



int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);

  if(argc<3) {
    if (mype == 0) printf("Usage: %s [opf data file] [num partitions][outer solve (opt)] \n",argv[0]);
    return 1;
  }
  if (mype == 0) cout << argv[0] << " starting..." << endl;  

  string dataname(argv[1]);
  int nPart 	= atoi(argv[2]);


  pipsOptions *pipsOpt = new pipsOptions();
  pipsOpt->readFile();
  pipsOpt->defGloOpt();  

  gisNLP = 0;
  
  //  gOuterSolve = 0:  
  //  gOuterSolve = 3:  Default solve - Schur complement based decomposition 
  //  gOuterSolve = 4:  BICG: use full mat as preconditioner
  //  gOuterSolve = 5:  BICG: Use Aggregation Preconditioner to solve 1st stage 

  if(argc>=4) {
    gOuterSolve = atoi(argv[3]);
  }
  assert (gOuterSolve>=3 && gOuterSolve<=5);

  //solve_usingCallbacks(dataname, nPart);
  solve(dataname, nPart);


  MPI_Finalize();
  return 0;
}


int solve(const string& dataname, int nPart) {
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  int nprocs; MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if(0==mype) cout << "Using a total of " << nprocs << " MPI processes." << endl;

  unsigned filelen;
  char *filedata;



  graphPart *graphP = new graphPart();
  DCPS *dcpowersys  = new DCPS();

  string fname = dataname;
  
  filelen = fname.length() + 1;
  filedata = new char[filelen];
  memcpy(filedata,dataname.c_str(),filelen);


  dcpowersys->PSReadDCData(filedata);
  dcpowersys->writeGraphFile(filedata,nPart);

  graphP->computeGPfromGraphFile(filedata);
  dcpowersys->setPartitioning(graphP->id_part,graphP->actparts,graphP->ncut_line);

  dcopflowInput* s = new dcopflowInput(dcpowersys,MPI_COMM_WORLD);
  if (mype == 0) cout << "dcopflowInput created .." << endl;

  if(gOuterSolve>=5){
	NlpPIPSIpmInterface<sFactoryAugAggregationPrecond, FilterIPMStochSolver, sNlpInfoFIX> pipsIpm(*s);		
	if (mype == 0) cout << "PIPSIpmInterface created .." << endl;
	pipsIpm.go();
  }else{
	NlpPIPSIpmInterface<sFactoryAug, FilterIPMStochSolver, sNlpInfoFIX> pipsIpm(*s);		
	if (mype == 0) cout << "PIPSIpmInterface created .." << endl;
	pipsIpm.go();
  }

  delete s;
  delete dcpowersys;
  delete[] filedata;
  delete graphP;
  
  return 0;
}


