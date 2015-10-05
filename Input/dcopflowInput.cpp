/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/
 

#include "dcopflowInput.hpp"

#include "ps.h"
#include "graphPart.h"

#include <fstream>
#include <sstream>


using namespace std;


void _fromTupleToCSC( const int ncol, int nnz, vector<int> colStart, 
							vector<int> rowIdx, vector<int> colIdx, vector<double> elts)
{
  int nnz_ori = nnz;

  vector<int> temp_rowID;
  temp_rowID.resize(nnz_ori);

  vector<int> next_ele_in_col;
  next_ele_in_col.resize(ncol+1);  
  
  for(int i=0; i<nnz_ori; i++){
	colStart[colIdx[i]+1]++;
  }
  for(int i=1; i<ncol+1; i++){
	colStart[i]+=colStart[i-1];
  }
  for(int i=0; i<ncol+1; i++){
	next_ele_in_col[i] = colStart[i];
  }

  for(int k=0;k<nnz_ori;k++){
	temp_rowID[next_ele_in_col[colIdx[k]]] = rowIdx[k];
	next_ele_in_col[colIdx[k]]++;
  }

  for(int i=0; i<ncol; i++){
	assert(next_ele_in_col[i] == colStart[i+1]);
  }  

  for(int i=0; i<nnz_ori; i++){
    rowIdx[i] = temp_rowID[i];
  }

}

void dcopflowInput::scenData::initialize(int nvar, int ncons) {
	collb.resize(nvar);
	colub.resize(nvar);
	objLin.resize(nvar);
	rowlb.resize(ncons);
	rowub.resize(ncons);
	rownames.resize(ncons);
	colnames.resize(nvar);
	didLoad = true;
}

void dcopflowInput::scenData::copyFrom(scenData s2) {
	collb=s2.collb;
	colub=s2.colub;
	objLin=s2.objLin;
	rowlb=s2.rowlb;
	rowub=s2.rowub;
	rownames=s2.rownames;
	colnames=s2.colnames;
	didLoad = s2.didLoad;
}



dcopflowInput::dcopflowInput(const DCPS* powersys, MPI_Comm comm) 
	: dcopf(NULL)
{
  MPI_Comm_rank(comm,&mype_);

  dcopf = new DCOPFLOW(powersys);
  
  nScenarios_= powersys->Nparts;
  
//  if(mype_ == 0)
  { 
    parseZeroData();
  }

}

dcopflowInput::~dcopflowInput() 
{
  if(dcopf) delete dcopf;
}

void dcopflowInput::parseZeroData()
{
  double *tempXl;	
  double *tempXU;

  nFirstStageVars_ = dcopf->Nvar_1st;
  nFirstStageCons_ = dcopf->Ncon_1st;

  nScenarios_ = dcopf->Nparts;
	
  ObjScale = 1.0;

  firstStageData.initialize(nFirstStageVars_, nFirstStageCons_);

  dcopf->ObjGradient_Lin_1st_Partition(&firstStageData.objLin[0]);
  dcopf->VarAndConBounds_1st_Partition(&firstStageData.collb[0],&firstStageData.colub[0],
										&firstStageData.rowlb[0],&firstStageData.rowub[0]);

  CoinBigIndex nnz;  
  vector<CoinBigIndex> colStart;
  vector<int> rowIdx;
  vector<int> colIdx;	
  vector<double> elts;

  // set Amat
  nnz = dcopf->GetJacNNZ_1st_Partition();	
  colStart.clear(); colStart.resize(nFirstStageVars_+1);  
  rowIdx.clear(); rowIdx.resize(nnz);
  colIdx.clear(); colIdx.resize(nnz);
  elts.clear();   elts.resize(nnz);

  dcopf->GetJac_1st_Partition(&rowIdx[0], &colIdx[0], &elts[0]);	
  _fromTupleToCSC(nFirstStageVars_,nnz,colStart,rowIdx,colIdx,elts);
  Amat.copyOf(true, nFirstStageCons_, nFirstStageVars_, nnz, &elts[0], &rowIdx[0], &colStart[0], 0);



  // set QAmat (not required). QTmat is empty, can be computed by getFirstStageHessian in stochInput.cpp
/*
  nnz = 0;		
  rowIdx.resize(nnz,0);
  colIdx.resize(nnz,0);
  elts.resize(nnz);
  QAmat.copyOf(true, dcopf->Nvar_1st,dcopf->Nvar_1st, 0, 0, 0, &colIdx[0], 0);
*/

  // prepare 2nd stage
  localData.resize(nScenarios_);  
  nSecondStageVars_.resize(nScenarios_);
  nSecondStageCons_.resize(nScenarios_);

  Tmat.resize(nScenarios_);
  Wmat.resize(nScenarios_);

//  QTmat.resize(nScenarios_);
  QWmat.resize(nScenarios_);

}




void dcopflowInput::loadLocalScenData(int scen) 
{
  if (localData[scen].didLoad) return;
  localData[scen].didLoad=true;

  nSecondStageVars_[scen] = dcopf->Nvar_2nd[scen];
  nSecondStageCons_[scen] = dcopf->Ncon_2nd[scen];  
  localData[scen].initialize(nSecondStageVars_[scen],nSecondStageCons_[scen]);
	
  dcopf->ObjGradient_Lin_2nd_Partition(scen, &localData[scen].objLin[0]);
  dcopf->VarAndConBounds_2nd_Partition(scen, &localData[scen].collb[0],&localData[scen].colub[0],
											&localData[scen].rowlb[0],&localData[scen].rowub[0]);

  CoinBigIndex nnz, nnz_link;
  vector<CoinBigIndex> colStart, colStart_link;  
  vector<int> rowIdx, rowIdx_link;
  vector<int> colIdx, colIdx_link;
  vector<double> elts, elts_link;

  // set Wmat and Tmat 
  nnz 		= dcopf->GetJacNNZ_2nd_Partition(scen);
  nnz_link 	= dcopf->GetJacNNZ_Link_Partition(scen);

  colStart.clear(); colStart.resize(nSecondStageVars_[scen]+1);	  
  rowIdx.clear(); rowIdx.resize(nnz);	
  colIdx.clear(); colIdx.resize(nnz);
  elts.clear(); elts.resize(nnz);  
  colStart_link.clear(); colStart_link.resize(nFirstStageVars_+1);	  
  rowIdx_link.clear(); rowIdx_link.resize(nnz_link);
  colIdx_link.clear(); colIdx_link.resize(nnz_link);
  elts_link.clear(); elts_link.resize(nnz_link);
  

  dcopf->GetJac_2nd_Link_Partition(	scen, &rowIdx[0], &colIdx[0], &elts[0], 
  									&rowIdx_link[0], &colIdx_link[0], &elts_link[0]);


  _fromTupleToCSC(nSecondStageVars_[scen],nnz,colStart,rowIdx,colIdx,elts);
  Wmat[scen] .copyOf(true, nSecondStageCons_[scen], nSecondStageVars_[scen], nnz, &elts[0], &rowIdx[0], &colStart[0], 0);

  _fromTupleToCSC(nFirstStageVars_,nnz_link,colStart_link,rowIdx_link,colIdx_link,elts_link);
  Tmat[scen] .copyOf(true, nSecondStageCons_[scen], nFirstStageVars_, nnz_link, &elts_link[0], &rowIdx_link[0], &colStart_link[0], 0);

//  doesn't work, cannot set total dim of matrix!
//  Wmat[scen] = CoinPackedMatrix(true, &rowIdx[0], &colIdx[0], &elts[0], nnz);
//  Tmat[scen] = CoinPackedMatrix(true, &rowIdx_link[0], &colIdx_link[0], &elts_link[0], nnz_link);

  // set QWmat.  QTmat is empty, can be computed by getSecondStageCrossHessian in stochInput.cpp
  nnz 		= dcopf->GetHesNNZ_2nd_Partition(scen);
  rowIdx.clear(); rowIdx.resize(nnz);
  colIdx.clear(); colIdx.resize(nnz);
  elts.clear();   elts.resize(nnz);

  dcopf->ObjGradient_Quad_2nd_Partition( scen, &rowIdx[0], &colIdx[0], &elts[0]);

  QWmat[scen] = CoinPackedMatrix(true, &rowIdx[0], &colIdx[0], &elts[0], nnz);
 

}

