/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "amplGenStochInput_AddSlack.hpp"
#include "AmplData_NL.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <limits>
#include <sstream>

#include "asl.h"
#include "asl_pfgh.h"
#include "getstub.h"
#include <climits>

#include "mpi.h"

using namespace std;

extern int gUseReducedSpace;

extern int gNP_Alg;


//  solve problem as 
//  min 1/s f(x_s)
//  st    g(x_s)  = 0
//         T_sx_s = x_0
//
//  this input is differenent from the amplGenStochInput one. 
//  In amplGenStochInput, we do split/reorder the Hessian and Jacobian in this case.
//  Here we need to add dummy variable/constraint


template <typename T>
void FreeAll( T & t ) {
    T tmp;
    t.swap( tmp );
}

amplGenStochInput_AddSlack::amplGenStochInput_AddSlack(const string &datarootname_in, 
					int overrideScenarioNumber, MPI_Comm comm) 
{
  useInputDate = 1;
  
  unsigned filelen;
  char *filedataArgv[2];
  	
  stringstream fname;
  std::string strTemp;
  haveRootNL=false;

  datarootname=datarootname_in;
  
  MPI_Comm_rank(comm,&mype_);	

  // Add the suffix   
  AmplSuffix *amplSuffix = new AmplSuffix();
  amplSuffix->DefineSuffix("pipsNLP_1stStageVar_in", Suffix_Var, Suffix_Int);  

  // no reduced space and graph partitioning  
  if (gUseReducedSpace!=0 || gNP_Alg!=0)
  	assert("adding slack setting is not available in reduced space setting.");
  
  nFirstStageCons_ = 0;
  nFirstStageVars_ = 0;
  
  {
	fname << datarootname << "0.nl";
	ifstream f(fname.str().c_str());

	if (!f.is_open()) {
	  haveRootNL=false;
	}
	else{
	  haveRootNL=true;
	  cout << "  Found " << fname.str() << " from processes " << mype_ << endl;	
	}
	f.close();
	
	assert(haveRootNL==false);
	
	if(mype_==1) cout << "  No first stage NL file! "<< endl;	
 	
	fname.str("");
	fname.clear();
  }

  nScenarios_ = overrideScenarioNumber;
  ObjScale = 1./nScenarios_;
  
  localData.resize(nScenarios_);  
  asl_i.resize(nScenarios_);
  nSecondStageVars_.resize(nScenarios_);
  nSecondStageCons_.resize(nScenarios_);
  LocGloVarMap.resize(nScenarios_);
  LocGloVarIdx.resize(nScenarios_);
  LocLocVarMap.resize(nScenarios_);

  map<int,int>::iterator it;

  // read 1st stage info
  int scen = 0;
//    if(mype_==scen)
  {
	  fname << datarootname << scen+1 << ".nl";
	  ifstream f(fname.str().c_str());
	  if (!f.is_open()) {
	    cout << "  Unable to open " << fname.str() << " from process " << mype_ << endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }else{
	    cout << "  Found " << fname.str() << " from processes " << mype_ << endl;	
	  }
      f.close();
	  
	  filelen = fname.str().length()+1;
	  filedataArgv[0] = new char[9];
      filedataArgv[1] = new char[filelen];
	  memcpy(filedataArgv[0],"pips_nlp",9);
	  memcpy(filedataArgv[1],fname.str().c_str(),filelen);
	  filedataArgv[2]=NULL;
 
	  // Allocate the ampl solver library (ASL) context.
      asl_i[scen] = localData[scen].initASL(filedataArgv, amplSuffix);  	
	  fname.str("");
	  fname.clear();
	  
	  localData[scen].didLoad=true;

	  LocGloVarIdx[scen] = NULL;
	  LocGloVarIdx[scen] = amplSuffix->GetSuffixVal_Int(asl_i[scen], "pipsNLP_1stStageVar_in", Suffix_Var);
	  //we must have 1st stage variable
	  assert(LocGloVarIdx[scen]);

	  //find correct index in C:  in ampl model, index of 1st stage var starts from 1, here we correct it as zero;
	  for (int j=0; j<asl_i[scen]->i.n_var_; j++){
		LocGloVarIdx[scen][j] -= 1;
	  }

	  //count number of 1st stage variables and creat variable map
      int n1stVarTemp=0, n2ndVarTemp=0;
	  for (int j=0; j<asl_i[scen]->i.n_var_; j++){
		if(LocGloVarIdx[scen][j]!=-1){ 
		  LocGloVarMap[scen].insert( pair<int,int>(j,LocGloVarIdx[scen][j]));
		  n1stVarTemp++;
		}else{
		  n2ndVarTemp++;
		}
	  }

	  // set 1st var ub & lb to inf
	  nFirstStageVars_ = n1stVarTemp;
	  assert(nFirstStageCons_==0);
	  firstStageData.initialize(nFirstStageVars_,nFirstStageCons_);
	  for(int j=0; j<nFirstStageVars_; j++){
		firstStageData.collb[j] = -1e25;
		firstStageData.colub[j] = 1e25;
		firstStageData.objGrad[j] = 0;
	  }

	  // allocate space for 2nd stages
	  assert(asl_i[scen]->i.n_var_ == n2ndVarTemp + nFirstStageVars_);
	  nSecondStageVars_[scen] = asl_i[scen]->i.n_var_;	  
      nSecondStageCons_[scen] = asl_i[scen]->i.n_con_ +  nFirstStageVars_;

	  localData[scen].initialize(nSecondStageVars_[scen],nSecondStageCons_[scen]);

	  // get 2nd var ub lb
	  for(int j=0; j<asl_i[scen]->i.n_var_; j++){
		localData[scen].collb[j] = asl_i[scen]->i.LUv_[j];
		localData[scen].colub[j] = asl_i[scen]->i.Uvx_[j];
		localData[scen].objGrad[j] = 1;
	  }

	  // reset 1st stage var bounds from 2nd var ub lb
	  for(it=LocGloVarMap[scen].begin(); it!=LocGloVarMap[scen].end(); it++){
		  firstStageData.collb[it->second] = localData[scen].collb[it->first];
		  firstStageData.colub[it->second] = localData[scen].colub[it->first];
	  }	  

	  // get 2nd con ub lb 
	  for(int j=0; j<asl_i[scen]->i.n_con_; j++){
		localData[scen].rowlb[j] = asl_i[scen]->i.LUrhs_[j];
		localData[scen].rowub[j] = asl_i[scen]->i.Urhsx_[j];
	  }
	  // for dummy constraint
	  for(int j=asl_i[scen]->i.n_con_; j<nSecondStageCons_[scen]; j++){
		localData[scen].rowlb[j] = 0;
		localData[scen].rowub[j] = 0;
	  }
  }

  //********************************************************************//
  LocWmatJacGoffMap.resize(nScenarios_);
  LocTmatJacGoffMap.resize(nScenarios_);
  Tmat.resize(nScenarios_);
  Wmat.resize(nScenarios_);

  LocQAmatHesGoffMap.resize(nScenarios_);
  LocQWmatHesGoffMap.resize(nScenarios_);
  LocQTmatHesGoffMap.resize(nScenarios_);  
  QTmat.resize(nScenarios_);
  QWmat.resize(nScenarios_);

  nnzQ2nd.resize(nScenarios_);
  nnzQCross2nd.resize(nScenarios_);  


  //********************************************************************//
  nnzEqJac1st=0; 
  nnzIneqJac1st=0;
  assert(nFirstStageCons_==0);
  amplRowMap1st = new int[nFirstStageCons_];
  
  nnzEqJac2nd.resize(nScenarios_);
  nnzIneqJac2nd.resize(nScenarios_);
  amplRowMap2nd.resize(nScenarios_);


  //********************************************************************//
  nnzA1st=0;
  nnzC1st=0;
  assert(nFirstStageCons_==0);
  assert(!firstStageData.locASL);

  nnzA2nd.resize(nScenarios_);
  nnzC2nd.resize(nScenarios_);

  nnzALink2nd.resize(nScenarios_,0);
  nnzCLink2nd.resize(nScenarios_,0);
  nnzALoc2nd.resize(nScenarios_,0);
  nnzCLoc2nd.resize(nScenarios_,0);

  JacALinkGoff2nd.resize(nScenarios_);
  JacCLinkGoff2nd.resize(nScenarios_);
  JacALocGoff2nd.resize(nScenarios_);
  JacCLocGoff2nd.resize(nScenarios_);  


  splitMatrices(0);
  getRowMap(0);
  getJacGoffMap(0);

  // set the root nl = NULL
  firstStageData.locASL = asl_i[0];

  dimsEqual = true;
  // initialize the number of vars/cons in 2nd stage
  for (int scen = 1; scen < nScenarios_; scen++) {
	nSecondStageCons_[scen] = 0;
	nSecondStageVars_[scen] = 0;
  }

  delete amplSuffix;
}


// split matrices to Amat and Wmat
void amplGenStochInput_AddSlack::splitMatrices(const int scen)
{

  ASL_pfgh * asl;
  int amplNz=0,	Wnnz  = 0, Tnnz  = 0, Annz  = 0, 
  	  			QWnnz = 0, QTnnz = 0, QAnnzFromScen=0;
  int tempNz = 0;

  vector<int> starts1stStTemp, starts2ndSt, startsLink;
  vector<int> rowIdx2ndSt, rowIdx1stSt, rowIdxLink;
  vector<double> elts2ndSt, elts1stSt, eltsLink;

  vector<int> amplGoffRowIdx;
  vector<double> amplGoffColIdx;
  vector<int> amplGoffColStarts;
  vector<int> amplGoffColLength;
  vector<double> amplGoffElts;

  vector<int> jColStart;
  vector<int> kRowIDX;
  vector<double> kelts;

  cgrad *cg;

  map<int,int>::iterator itVar, itVar_Row; 
  int wrk1stGoff=0, wrk2ndGoff=0, wrkLinkGoff=0;
  int wrkGoff=0;

  int wrkGoffStart=0, addDummyCon=0, nextSearchVar=0;

  if(0==scen){
////////////////////////////	 build Amat   ////////////////////////////////////////////////////////
    // set empty Amat. This is the matrix difined by 1st stage cons ---  this is empty in Curl's setting
    assert(nFirstStageCons_==0);
    jColStart.resize(nFirstStageVars_+1);
    for (int j = 0; j <= nFirstStageVars_; j++) {
	  jColStart[j] = 0;
    }
    Amat.copyOf (true, nFirstStageCons_, nFirstStageVars_, 0, 0, 0, &jColStart[0], 0);

  ////////////////////////////		 build QAmat  ////////////////////////////////////////////////////////
    jColStart.clear();
    jColStart.resize(nFirstStageVars_+1,0);
    QAmat.copyOf(true, nFirstStageVars_, nFirstStageVars_, 0, 0, 0, &jColStart[0], 0);
  }

  // set sparse structure of 2nd stage Jac and linking Jac
  asl = localData[scen].locASL;

  Wnnz  = 0; Tnnz  = 0; 
  QWnnz = 0;  QTnnz = 0; QAnnzFromScen=0;
  wrk1stGoff=0; wrk2ndGoff=0; wrkLinkGoff=0;
  wrkGoff=0;
	
  ////////////////////////////	build Wmat   	////////////////////////////////////////////////////////
  amplNz = nzc;
  Wnnz = amplNz + nFirstStageVars_;
  jColStart.clear();
  jColStart.resize(nSecondStageVars_[scen]+1,0);
  kRowIDX.resize(Wnnz);
  kelts.resize(Wnnz);

  amplGoffRowIdx.clear(); amplGoffRowIdx.resize(amplNz,0);
  amplGoffElts.clear(); amplGoffElts.resize(amplNz,0.);
  amplGoffColStarts.clear(); amplGoffColStarts.resize(n_var+1,0);
  assert(nSecondStageVars_[scen]==n_var);

  // get info from ampl
  for( int j=0; j < n_con; j++){
	for(cg = Cgrad[j]; cg; cg = cg->next){
	  amplGoffRowIdx[cg->goff] = j;
	  amplGoffElts[cg->goff] = cg->coef;
	  amplGoffColStarts[cg->varno+1]++;
	}
  }	
  for(int j=1; j < n_var+1; j++){
	amplGoffColStarts[j] += amplGoffColStarts[j-1];
  }
  assert(amplGoffColStarts[n_var]==amplNz);


  wrkGoffStart=0; wrkGoff=0; addDummyCon=0; nextSearchVar=0;
  for(itVar=LocGloVarMap[scen].begin(); itVar!=LocGloVarMap[scen].end(); itVar++){

	for(int j=wrkGoffStart; j<amplGoffColStarts[itVar->first+1]; j++){
	  kelts[wrkGoff] = amplGoffElts[j];
	  kRowIDX[wrkGoff] = amplGoffRowIdx[j];
	  LocWmatJacGoffMap[scen].insert( pair<int,int>(wrkGoff,j));
	  wrkGoff++;
	}
	  
	// correct  col start
	for(int j=nextSearchVar; j<itVar->first+1; j++){
	  jColStart[j] = amplGoffColStarts[j] + addDummyCon;
	}
	  
	//add elt of the additional constraint;
	kelts[wrkGoff] = 1;
	kRowIDX[wrkGoff] = n_con + addDummyCon;
	wrkGoff++;
	addDummyCon++;

	// correct next col start
	nextSearchVar = itVar->first+1;
	  
	wrkGoffStart = amplGoffColStarts[itVar->first+1];
  }
  assert(addDummyCon == nFirstStageVars_);


  if(wrkGoffStart!=amplNz){
	for(int j=wrkGoffStart; j<amplNz; j++){
	  kelts[wrkGoff] = amplGoffElts[j];
	  kRowIDX[wrkGoff] = amplGoffRowIdx[j];
	  LocWmatJacGoffMap[scen].insert( pair<int,int>(wrkGoff,j));		
	  wrkGoff++;
	}
  }
  for(int j=nextSearchVar; j<n_var+1; j++, nextSearchVar++){
	jColStart[j] = amplGoffColStarts[j] + addDummyCon;
  } 	
  assert(wrkGoff == Wnnz);
  assert(jColStart[n_var] == Wnnz);
  assert(nextSearchVar == n_var+1);


  Wmat[scen].copyOf(true, nSecondStageCons_[scen], nSecondStageVars_[scen], Wnnz, 
				&kelts[0], &kRowIDX[0], &jColStart[0], 0);


  ////////////////////////////	build Tmat	 ////////////////////////////////////////////////////////
  Tnnz = nFirstStageVars_;
  jColStart.clear(); 	jColStart.resize(nFirstStageVars_+1,0);
  kRowIDX.clear(); 	kRowIDX.resize(Tnnz);
  kelts.clear();		kelts.resize(Tnnz);
	
  addDummyCon=0;
	
  for(int j=0; j<nFirstStageVars_+1; j++){
    jColStart[j] = j;
  }
  for(int j=0; j<Tnnz; j++){
	kelts[j] = -1;
  }	
	
  for(itVar=LocGloVarMap[scen].begin(); itVar!=LocGloVarMap[scen].end(); itVar++){
	kRowIDX[itVar->second] =  n_con + addDummyCon++;
  }
  assert(addDummyCon == nFirstStageVars_);
	
  Tmat[scen].copyOf(true, nSecondStageCons_[scen], nFirstStageVars_, Tnnz, 
						&kelts[0], &kRowIDX[0], &jColStart[0], 0);

  ////////////////////////////		 build QTmat  ////////////////////////////////////////////////////////
  jColStart.clear();
  jColStart.resize(nFirstStageVars_+1,0);
  QTmat[scen].copyOf(true, nSecondStageVars_[scen], nFirstStageVars_, 0, 0, 0, &jColStart[0], 0);

  nnzQCross2nd[scen] = 0;
  ////////////////////////////		 build QWmat	////////////////////////////////////////////////////////

  QWnnz = localData[scen].Ampl_nnz_Hessian_Tri();

  jColStart.clear();	jColStart.resize( nSecondStageVars_[scen]+1,0);
  kRowIDX.clear();	kRowIDX.resize(QWnnz);
  kelts.clear();		kelts.resize(QWnnz,0);

  amplGoffRowIdx.clear(); amplGoffRowIdx.resize(QWnnz);
  amplGoffColIdx.clear(); amplGoffColIdx.resize(QWnnz);
  amplGoffElts.clear();
  amplGoffColStarts.clear(); amplGoffColStarts.resize(nSecondStageVars_[scen]+1,0);

  for(int j=0; j<QWnnz;j++){
	amplGoffRowIdx[j] = sputinfo->hrownos[j];
  }
  for(int j=0; j<nSecondStageVars_[scen]+1;j++){
	amplGoffColStarts[j] = sputinfo->hcolstarts[j];
  }

  int findQele=0;
  int *defineEleID = new int[ nSecondStageVars_[scen]];

  // ampl use upper triangle, but transpose to lower triangle for OOQP, both require col-wise form
  for( int j = 0; j < nSecondStageVars_[scen]; j++ ) {
	for( int k = amplGoffColStarts[j]; k < amplGoffColStarts[j+1]; k++ ) {
	 int QrowIDX = amplGoffRowIdx[k];
	  // use only the upper triangle part
	  assert(QrowIDX <= j);
	  jColStart[QrowIDX+1]++;
	  findQele++;
	}
  }
  assert(findQele==QWnnz);

  for(int j=1; j < nSecondStageVars_[scen]+1; j++){
	jColStart[j] += jColStart[j-1];
	defineEleID[j-1] = jColStart[j-1];
  }
  assert(jColStart[nSecondStageVars_[scen]]==QWnnz);
	
  for( int j = 0; j < nSecondStageVars_[scen]; j++ ) {
	for( int k = amplGoffColStarts[j]; k < amplGoffColStarts[j+1]; k++ ) {
	  int GoffInLowerPart = defineEleID[amplGoffRowIdx[k]];
	  kRowIDX[GoffInLowerPart]	=  j;
	  LocQWmatHesGoffMap[scen].insert( pair<int,int>(GoffInLowerPart,k));	
	  defineEleID[amplGoffRowIdx[k]]++;
	}		  
  }
	
  delete [] defineEleID;

  QWmat[scen].copyOf(true, nSecondStageVars_[scen], nSecondStageVars_[scen], QWnnz,
					&kelts[0], &kRowIDX[0], &jColStart[0], 0);

  nnzQ2nd[scen] = QWnnz;

  amplGoffRowIdx.clear();
  amplGoffColStarts.clear(); 

  if(scen == 0){
	nnzQ1st = QAnnzFromScen;
  }

  //////////////////////////////////////////////////////////////////////////////////////////

  FreeAll(starts2ndSt); FreeAll(starts1stStTemp); FreeAll(startsLink);
  FreeAll(rowIdx2ndSt); FreeAll(rowIdx1stSt); FreeAll(rowIdxLink);
  FreeAll(elts2ndSt); FreeAll(elts1stSt); FreeAll(eltsLink);

  FreeAll(amplGoffRowIdx);
  FreeAll(amplGoffColIdx);
  FreeAll(amplGoffColStarts);
  FreeAll(amplGoffColLength);
  FreeAll(amplGoffElts);
	  
}


void amplGenStochInput_AddSlack::getRowMap(const int scen)
{
    int loc_my=0, loc_mz=0;
	amplRowMap2nd[scen] = new int[asl_i[scen]->i.n_con_];
	
	for( int i = 0; i < asl_i[scen]->i.n_con_; i++ ) {
      if ( localData[scen].rowlb[i] == localData[scen].rowub[i] ) {
        amplRowMap2nd[scen][i] = - (loc_my + 1); // Negative values indicate an equality
      	loc_my++;
      } else {
        amplRowMap2nd[scen][i] = loc_mz;
      	loc_mz++;
      }
    }
}



void amplGenStochInput_AddSlack::getJacGoffMap(const int scen)
{
  ASL_pfgh * asl;
  cgrad *cg;
  
  int amplNz;
  map<int,int>::iterator itVar;

  asl = localData[scen].locASL;
  amplNz = nzc;
	
  for(int j=0; j<n_con;j++){
    if( amplRowMap2nd[scen][j] < 0){
	  for(cg = Cgrad[j]; cg; cg = cg->next){
	    nnzA2nd[scen]++;
	  }
    }else{
	  for(cg = Cgrad[j]; cg; cg = cg->next){
	    nnzC2nd[scen]++;	  
	  }
    }
  }

  int Tnnz = Tmat[scen].getNumElements();
  int Wnnz = Wmat[scen].getNumElements();
 
  const int* tRowIDX 	= Tmat[scen].getIndices();   
  const int* wRowIDX 	= Wmat[scen].getIndices();    

  int rowID;
  for( int k = 0; k < Tnnz; k++ ) {
	rowID = tRowIDX[k];
	if(rowID<n_con){
	  if( amplRowMap2nd[scen][rowID] < 0){
	    nnzALink2nd[scen]++;
	  }else{
	    nnzCLink2nd[scen]++;
	  }
	}
  }
  for( int k = 0; k < Wnnz; k++ ) {
    rowID = wRowIDX[k];
	if(rowID<n_con){
	  if( amplRowMap2nd[scen][rowID] < 0){
	    nnzALoc2nd[scen]++;
	  }else{
	    nnzCLoc2nd[scen]++;
	  }
	}
  }
  assert( amplNz  == nnzALoc2nd[scen] + nnzCLoc2nd[scen]);  

  // add eles from dummy constraint
  nnzALoc2nd[scen] += nFirstStageVars_;
  nnzALink2nd[scen] += nFirstStageVars_; 

  int kA = 0, kC = 0, rowIdx=0;
  for(int j=0; j<n_con;j++){
	if( amplRowMap2nd[scen][j] < 0){	
	  for(cg = Cgrad[j]; cg; cg = cg->next){  
	  	JacALocGoff2nd[scen].insert( pair<int,int>(cg->goff,kA++));
	  }
    }else{
	  for(cg = Cgrad[j]; cg; cg = cg->next){
		JacCLocGoff2nd[scen].insert( pair<int,int>(cg->goff,kC++));	
	  }
    }
  }
  assert( kA == nnzA2nd[scen] );
  assert( kC == nnzC2nd[scen] );
  assert( kA+kC==amplNz);
  nnzA2nd[scen]+= nFirstStageVars_;

}

void amplGenStochInput_AddSlack::		loadLocalNLdata(int scen)
{
  
  if (localData[scen].didLoad) return;
  localData[scen].didLoad=true;

  // Add the suffix  
  AmplSuffix *amplSuffix = new AmplSuffix();
  amplSuffix->DefineSuffix("pipsNLP_1stStageVar_in", Suffix_Var, Suffix_Int);  

	  stringstream fname;
	  char *filedataArgv[2];
	  filedataArgv[0] = new char[9];
	  memcpy(filedataArgv[0],"pips_nlp",9);	  
	  unsigned filelen;
	  map<int,int>::iterator it;
	  std::string strTemp;

  	  fname << datarootname << scen+1 << ".nl";
	  ifstream f(fname.str().c_str());
//	  cout << "Try to find " << fname.str() << " from processes " << mype_ << endl;  	  
	  if (!f.is_open()) {
	    cout << "  Unable to open " << fname.str() << " from process " << mype_ << endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	  }else{
	    cout << "  Found " << fname.str() << " from processes " << mype_ << endl;	
	  }
      f.close();

	  filelen = fname.str().length()+1;
      filedataArgv[1] = new char[filelen];
	  memcpy(filedataArgv[1],fname.str().c_str(),filelen);
	  filedataArgv[2]=NULL;
 
	  // Allocate the ampl solver library (ASL) context.
      asl_i[scen] = localData[scen].initASL(filedataArgv, amplSuffix);  	
	  fname.str("");
	  fname.clear();
	  
	  LocGloVarIdx[scen] = NULL;
	  LocGloVarIdx[scen] = amplSuffix->GetSuffixVal_Int(asl_i[scen], "pipsNLP_1stStageVar_in", Suffix_Var);

	  //we must have 1st stage variable
	  assert(LocGloVarIdx[scen]);

	  //find correct index in C:  in ampl model, index of 1st stage var starts from 1, here we correct it as zero;
	  for (int j=0; j<asl_i[scen]->i.n_var_; j++){
		LocGloVarIdx[scen][j] -= 1;
	  }

	  //creat variable map
      int n1stVarTemp=0, n2ndVarTemp=0;
	  for (int j=0; j<asl_i[scen]->i.n_var_; j++){
		if(LocGloVarIdx[scen][j]!=-1){ 
	      //creat local NL to global var map
	      LocGloVarMap[scen].insert( pair<int,int>(j,LocGloVarIdx[scen][j]));
		  n1stVarTemp++;
		}else{
		  n2ndVarTemp++;
		}
	  }

  // allocate space for 2nd stages
  assert(nFirstStageVars_ == n1stVarTemp);
  nSecondStageVars_[scen] = asl_i[scen]->i.n_var_;
  nSecondStageCons_[scen] = asl_i[scen]->i.n_con_ +  nFirstStageVars_; 
  
  if(dimsEqual==true){
  	if(nSecondStageVars_[scen]!=nSecondStageVars_[0] || nSecondStageCons_[scen]!=nSecondStageCons_[0]) dimsEqual=false;
  }



  localData[scen].initialize(nSecondStageVars_[scen],nSecondStageCons_[scen]);

  	  // get 2nd var ub lb
	  for(int j=0; j<asl_i[scen]->i.n_var_; j++){
		localData[scen].collb[j] = asl_i[scen]->i.LUv_[j];
		localData[scen].colub[j] = asl_i[scen]->i.Uvx_[j];
		localData[scen].objGrad[j] = 1;
	  }

	  // get 2nd con ub lb 
	  for(int j=0; j<asl_i[scen]->i.n_con_; j++){
		localData[scen].rowlb[j] = asl_i[scen]->i.LUrhs_[j];
		localData[scen].rowub[j] = asl_i[scen]->i.Urhsx_[j];
	  }
	  // for dummy constraint
	  for(int j=asl_i[scen]->i.n_con_; j<nSecondStageCons_[scen]; j++){
		localData[scen].rowlb[j] = 0;
		localData[scen].rowub[j] = 0;
	  }

  splitMatrices(scen);
  getRowMap(scen);
  getJacGoffMap(scen);

  delete amplSuffix;
}









void amplGenStochInput_AddSlack::doNetworkPart(const int scen, AmplSuffix* amplSuffix )
{
  assert("not done" &&0);
}


