/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */
#include <stdio.h>
#include <stdlib.h>

#include "rawInput.hpp"
#include "PIPSIpmInterface.h"
#include "OOQPRecourseInterface.hpp"
#include "QpGenSparseMa57.h"
#include "QpGenSparseMa27.h"
#include "MehrotraSolver.h"
#include "GondzioSolver.h"
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

  //if(mype==0) cout << argv[0] << " starting ..." << endl;
  
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
  //printf("mype=%d mynewpe=%d\n", mype, myBatchPe);
 
  {
    stringstream ss; ss << datadirname << (color+1) << "/" << datarootname;
    datarootname=ss.str();
  }
  {
    stringstream ss; ss << datadirname << (color+1) << "/mean-prob/" << datarootnameMean;
    datarootnameMean=ss.str();
  }

  printf("mype=%d mean batch %d  prob.from [%s]\n", mype, color+1, datarootnameMean.c_str());
  rawInput* sMean = new rawInput(datarootnameMean, 1, MPI_COMM_SELF);
  //rawInput* sMean = new rawInput(datarootname, 4, MPI_COMM_SELF);
  std::vector<double> firstStageSol;
  std::vector<double> fstStageDual;
  {
    //PIPSIpmInterface<sFactoryAugSchurLeaf, MehrotraStochSolver> pipsIpm(*sMean, MPI_COMM_SELF);
    PIPSIpmInterface<sFactoryAug, MehrotraStochSolver> pipsIpm(*sMean, MPI_COMM_SELF);
    delete sMean;
    
    pipsIpm.go();   

    firstStageSol = pipsIpm.getFirstStagePrimalColSolution();
    fstStageDual  = pipsIpm.getFirstStageDualRowSolution();
    if(myBatchPe==0)
      printf("mype=%d mean batch %d   1stStageObjective=%20.12f TotalObjective=%20.12f\n", 
	     mype, color+1, 
	     pipsIpm.getFirstStageObjective(),
	     pipsIpm.getObjective());
  }

  // save the 1st stage solution for each batch
  if(myBatchPe==0) {
    {
	stringstream ss2; ss2<<outputdir << "/batch-" << (1+color) << "-out_primal_1stStage.txt";
	string sFile=ss2.str();
	cout << "saving 1st stage sol to " << sFile << endl;
	ofstream file1stStg(sFile.c_str());
	file1stStg << scientific; file1stStg.precision(16);
	for(size_t i=0; i<firstStageSol.size(); i++)
	    file1stStg << firstStageSol[i] << endl;
	file1stStg.close();
    }
    sleep(2);
    {
	stringstream ss2; ss2<<outputdir << "/batch-" << (1+color) << "-out_dual_1stStage.txt";
	string sFile=ss2.str();
	cout << "saving 1st stage dual sol to " << sFile << endl;
	ofstream file1stStg(sFile.c_str());
	file1stStg << scientific; file1stStg.precision(16);
	for(size_t i=0; i<firstStageSol.size(); i++)
	    file1stStg << fstStageDual[i] << endl;
	file1stStg.close();
    }

  }
  MPI_Barrier(commBatch);
  

  rawInput* s = new rawInput(datarootname, nscen, commBatch);

  //load 1st stage solution
  {
    stringstream ss2; ss2<<outputdir << "/batch-" << (1+color) << "-out_primal_1stStage.txt";
    string filename = ss2.str();
    ifstream file(filename.c_str());
    assert(file.is_open());

    firstStageSol.resize(s->nFirstStageVars());
    for(int d=0; d<s->nFirstStageVars(); d++)
      file >> firstStageSol[d];
    //cout << "first stage sol loaded from file" << endl;
  }
  
  
  for(int scen=0; scen<nscen; scen++) {
    if( scen%nbatchprocs == myBatchPe ) {
      printf("Proc [%d][%d] does scen [%d] in batch [%s]\n", 
	     mype, myBatchPe, scen, datarootname.c_str());
      
      pips_options::setIntParameter("OUTER_SOLVE", 1); //iter.refin.
      OOQPRecourseInterface<MehrotraSolver,QpGenSparseMa27> ooqpRecourse(*s, scen, firstStageSol);
      ooqpRecourse.go();

      printf("Proc [%d][%d] scen %d in batch %d    RecourseObjective=%g\n",
	     mype, myBatchPe, scen, color+1, ooqpRecourse.getObjective());

      //////////////////////////////////////////////////////////
      // save primal recourse solution
      //////////////////////////////////////////////////////////
      {
	stringstream ss; ss << outputdir << "/batch-" << (1+color) << "-out_primals_recou" << (scen+1) << ".txt";
	cout << "saving recourse pb primals to " << ss.str() << endl;
	ofstream fileprimals(ss.str().c_str());
	std::vector<double> primals = ooqpRecourse.getPrimalColSolution();
	for(size_t i=0; i<primals.size(); i++)
	  fileprimals << primals[i] << endl;
	fileprimals.close();
      }
      //////////////////////////////////////////////////////////
      // save dual solution of the recourse problem
      //////////////////////////////////////////////////////////
      {
	stringstream ss; ss << outputdir << "/batch-" << (1+color) << "-out_duals_recou" << (scen+1) << ".txt";
	cout << "saving recourse pb duals to " << ss.str() << endl;
	ofstream fileduals(ss.str().c_str());
	std::vector<double> duals = ooqpRecourse.getDualRowSolution();
	for(size_t i=0; i<duals.size(); i++)
	  fileduals << duals[i] << endl;
	fileduals.close();
      }
    }
  }
  delete s;
 

  MPI_Finalize();
  return 0;
}
