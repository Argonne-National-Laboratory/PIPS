/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */
#include <stdio.h>
#include <stdlib.h>

#include "rawInput.hpp"
#include "PIPSIpmInterface.h"
#include "OOQPRecourseInterface.hpp"
#include "QpGenSparseMa57.h"
#include "MehrotraSolver.h"
#include "sFactoryAugSchurLeaf.h"
#include "MehrotraStochSolver.h"

#include <string>
#include <sstream>

using namespace std;


int main(int argc, char ** argv) {
  MPI_Init(&argc, &argv);
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);

  if(argc<6) {
    if (mype == 0) printf("Usage: %s [rawdump batches directory] [rawdump root name] [num batches] [num scenarios per batch] [output dir]\n", argv[0]);
    return 1;
  }

  if(mype==0) cout << argv[0] << " starting ..." << endl;
  
  string datadirname(argv[1]);
  string datarootname(argv[2]);
  string datarootnameMean=datarootname;
  int nbatch = atoi(argv[3]);
  int nscen = atoi(argv[4]);
  string outputdir(argv[5]);

  int nprocs; MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if(nprocs<nbatch) {
    if(mype == 0) cout << "Number of procs has to be greater or equal to number of batches" << endl;
    return 1;
  }

  if( nbatch* (nprocs/nbatch) != nprocs ) {
    if(mype == 0) cout << "Number of procs [" << nprocs << "] has to be multiple of the number of batches[" << nprocs/nbatch << "]" << endl;
    return 1;
  }

  if( (nprocs/nbatch) * (nscen/(nprocs/nbatch)) != nscen ) {
    if(mype == 0) cout << "Number of scens [" << nscen << "] has to be a multiple of the number of processes per batch [" << nprocs/nbatch << "]" << endl;
    return 1;
  }
  int nbatchprocs=nprocs/nbatch;
  int color = mype/nbatchprocs; 
  
  MPI_Comm commBatch;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &commBatch);
  assert(commBatch!=MPI_COMM_NULL);
  int myBatchPe; MPI_Comm_rank(commBatch, &myBatchPe);
  printf("mype=%d mynewpe=%d\n", mype, myBatchPe);
 
  {
    stringstream ss; ss << datadirname << (color+1) << "/" << datarootname;
    datarootname=ss.str();
  }
  {
    stringstream ss; ss << datadirname << (color+1) << "/mean-prob/" << datarootnameMean;
    datarootnameMean=ss.str();
  }

  printf("mype=%d mean prob.from [%s]\n", mype, datarootnameMean.c_str());
  rawInput* sMean = new rawInput(datarootnameMean, 1, MPI_COMM_SELF);
  //rawInput* sMean = new rawInput(datarootname, 4, MPI_COMM_SELF);
  std::vector<double> firstStageSol;
  /*
  {
    PIPSIpmInterface<sFactoryAugSchurLeaf, MehrotraStochSolver> pipsIpm(*sMean, MPI_COMM_SELF);
    delete sMean;
    
    pipsIpm.go();   

    firstStageSol = pipsIpm.getFirstStagePrimalColSolution();
  }

  // save the 1st stage solution for each batch
  if(myBatchPe==0) {
    stringstream ss2; ss2<<outputdir << "/batch-" << (1+color) << "-out_primal_1stStage.txt";
    string sFile=ss2.str();
    cout << "saving 1st stage sol to " << sFile << endl;
    ofstream file1stStg(sFile.c_str());
    file1stStg << scientific; file1stStg.precision(16);
    for(size_t i=0; i<firstStageSol.size(); i++)
      file1stStg << firstStageSol[i] << endl;
    file1stStg.close();
  }

  */
  rawInput* s = new rawInput(datarootname, nscen, commBatch);

  //code to load the solution
  {
    stringstream ss2; ss2<<outputdir << "/batch-" << (1+color) << "-out_primal_1stStage.txt";
    string filename = ss2.str();
    ifstream file(filename.c_str());
    assert(file.is_open());

    firstStageSol.resize(s->nFirstStageVars());
    for(int d=0; d<s->nFirstStageVars(); d++)
      file >> firstStageSol[d];
    cout << "first stage sol loaded from file" << endl;
  }
  
  
  for(int scen=0; scen<nscen; scen++) {
    if( scen%nbatchprocs == myBatchPe ) {
      printf("Proc [%d][%d] does scen [%d] in batch [%s]\n", 
	     mype, myBatchPe, scen, datarootname.c_str());
      
      OOQPRecourseInterface<MehrotraSolver,QpGenSparseMa57> ooqpRecourse(*s, scen, firstStageSol);
      ooqpRecourse.go();
    }
  }
  delete s;


  // if(mype==0) cout << "Saving solution" << endl;
  // for(int s=0; s<nscen; s++) {
    
  //   std::vector<double> duals = pipsIpm.getSecondStageDualRowSolution(s);
  //   if(duals.size()) {
  //     stringstream ss1; ss1<<outputdir << "/batch-" << (1+color) << "-out_duals_scen"<<(s+1)<<".txt";
  //     cout << "saving duals to " << ss1.str() << endl;
  //     ofstream fileduals(ss1.str().c_str());
      
  //     for(size_t i=0; i<duals.size(); i++)
  // 	fileduals << duals[i] << endl;
  //     fileduals.close();
  //   }
    
  //   std::vector<double> primals = pipsIpm.getSecondStagePrimalColSolution(s);
  //   if(primals.size()) {
  //     stringstream ss2; ss2<<outputdir << "/batch-" << (1+color) << "-out_primals_scen"<<(s+1)<<".txt";
  //     cout << "saving primals to " << ss2.str() << endl;
  //     ofstream fileprimals(ss2.str().c_str());
  //     std::vector<double> primals = pipsIpm.getSecondStagePrimalColSolution(s);
  //     for(size_t i=0; i<primals.size(); i++)
  // 	fileprimals << primals[i] << endl;
  //     fileprimals.close();
  //   }
  // }
  // if(mynewpe==0) {
  //   std::vector<double> firstStageSol = pipsIpm.getFirstStagePrimalColSolution();
  //   stringstream ss2; ss2<<outputdir << "/batch-" << (1+color) << "-out_primal_1stStage.txt";
  //   string sFile=ss2.str();
  //   cout << "saving 1st stage sol to " << sFile << endl;
  //   ofstream file1stStg(sFile.c_str());
  //   for(size_t i=0; i<firstStageSol.size(); i++)
  //     file1stStg << firstStageSol[i] << endl;
  //   file1stStg.close();
  // }
  // cout << "Solution saved" << endl;
 

  MPI_Finalize();
  return 0;
}

