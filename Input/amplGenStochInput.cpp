/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "amplGenStochInput.hpp"
#include "AmplData_NL.hpp"

#include <iostream>
#include <cstring>
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

//  solve problem as general stochastic programming problem
//  min f(x_0) + \sum f(x_i)
//  st    g(x_i)  = 0
//
//  we do split/reorder the Hessian and Jacobian in this case.

template <typename T>
void FreeAll( T & t ) {
    T tmp;
    t.swap( tmp );
}



// here first is the row/col ID, and  second is the col/row ID in row-wise/col-wise setting (Note that first and second are index of col/row, 
// they are not goff!)

void RowColWise_map( int first[], int nnz, int second[], int newGoff[])
{
  int fi, se, j, k, kinc, inc;
  double dtemp;
  int dtemp_goff;
  const int incs[]  = {1, 5, 19, 41, 109, 209, 505,
		       929, 2161, 3905, 8929, 16001, INT_MAX};
  
  for ( k = 0; incs[k] <= nnz/2; k++ ) ;

  kinc = k - 1;

  for( ; kinc >= 0; kinc-- ) {
	// Loop over all increments
	inc = incs[kinc];

	for ( k = inc; k < nnz; k++ ) {
	  dtemp_goff = newGoff[k];
	  fi = first[ k ];
	  se = second[ k];
	  for( j = k; j >= inc; j -= inc ) {
	  	if ( fi < first[j - inc] || 
			 ( fi == first[j - inc] &&
		       se < second[ j - inc ]) ) {
		  first[j]   	= first[j - inc];
		  second[j]  	= second[j - inc];
		  newGoff[j] = newGoff[j - inc];
		} else {
		  break;
		}
	  } 
	  first[j]   	= fi;
	  second[j]  	= se;
	  newGoff[j] = dtemp_goff;
	}
  } // End loop over all increments
}


amplGenStochInput::amplGenStochInput(const string &datarootname_in, 
					int overrideScenarioNumber, MPI_Comm comm) 
{
  useInputDate = 1;
  
  unsigned filelen;
  char *filedataArgv[3];
  	
  stringstream fname;
  std::string strTemp;
  haveRootNL=false;

  datarootname=datarootname_in;
  
  MPI_Comm_rank(comm,&mype_);	

  // Add suffixes   
  AmplSuffix *amplSuffix = new AmplSuffix();
  amplSuffix->DefineSuffix("pipsNLP_1stStageVar_in", Suffix_Var, Suffix_Int);  
  
  if(gUseReducedSpace>0){	
  	amplSuffix->DefineSuffix("pipsNLP_DecisionVar_in", Suffix_Var, Suffix_Int);  
  }

  if(gNP_Alg!=0){
  	amplSuffix->DefineSuffix("pipsNLP_VarPartIdx_in", Suffix_Var, Suffix_Int); 
    amplSuffix->DefineSuffix("pipsNLP_ConPartIdx_in", Suffix_Con, Suffix_Int); 	
  }
  
  nFirstStageCons_=0;
  nFirstStageVars_=0;
    
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

  if(2==gUseReducedSpace){  
	decisionVarDim.resize(nScenarios_,0);
	schurVarConIDinNL.resize(nScenarios_); 
  }

  for (int scen = 0; scen < 1; scen++) {	
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
	  cout << "MPI#" << mype_ << ":  Start read pfgh_read! " << endl;  
      asl_i[scen] = localData[scen].initASL(filedataArgv, amplSuffix); 
	  cout << "MPI#" << mype_ << ":  Scenario " << scen << "   pfgh_read Done!" << endl; 
	  
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
	  
	  //creat variable map
      int n1stVarTemp=0, n2ndVarTemp=0;
	  for (int j=0; j<asl_i[scen]->i.n_var_; j++){
		if(LocGloVarIdx[scen][j]!=-1){ 
	      //creat local NL to global var map
	      LocGloVarMap[scen].insert( pair<int,int>(j,LocGloVarIdx[scen][j]));
		  n1stVarTemp++;
		}else{
		  //creat local NL to local var map
	      LocLocVarMap[scen].insert( pair<int,int>(j,n2ndVarTemp));
		  n2ndVarTemp++;
		}
	  }

	  // define decision var
	  if(gUseReducedSpace>0){  	
	  	decisionVarDim[scen]=0;
		int *decisionVarIDX = amplSuffix->GetSuffixVal_Int(asl_i[scen], "pipsNLP_DecisionVar_in", Suffix_Var);
		//find correct index in C:	in ampl model, index of decision var starts from 1, here we correct it as zero;
		for(int j=0; j<asl_i[scen]->i.n_var_; j++){
		  decisionVarIDX[j] -=1;
		  if(decisionVarIDX[j]>=0) 
			decisionVarDim[scen]++;
		}
		schurVarConIDinNL[scen] = (int*)malloc(decisionVarDim[scen]*sizeof(int));
	  
		for(int findSCVar=0; findSCVar<decisionVarDim[scen];findSCVar++){
		  for(int j=0; j<asl_i[scen]->i.n_var_; j++){
			if(decisionVarIDX[j] == findSCVar){
			  map<int,int>::iterator itVarTemp = LocLocVarMap[scen].find(j);	 
			  assert(itVarTemp != LocLocVarMap[scen].end());
			  schurVarConIDinNL[scen][findSCVar]=itVarTemp->second;
			  break;
			}
		  }
		}
	  }

	  //get 1st stage var since we do not have root nl file, we set up 1st stage infomation from the 1st scenario
      if(scen==0 && haveRootNL==false){
	  	assert(nFirstStageVars_==0 && nFirstStageCons_==0);
		
		nFirstStageVars_ = n1stVarTemp;
		
		firstStageData.initialize(nFirstStageVars_,nFirstStageCons_);

		// get 1st var ub lb
		for(it=LocGloVarMap[scen].begin(); it!=LocGloVarMap[scen].end(); it++){
	      firstStageData.collb[it->second] = asl_i[scen]->i.LUv_[it->first];
	      firstStageData.colub[it->second] = asl_i[scen]->i.Uvx_[it->first];
	      firstStageData.objGrad[it->second] = 0;
		}
	  }  

	  assert(nFirstStageVars_ == n1stVarTemp);
	  assert(nFirstStageCons_ == 0);

	  nSecondStageVars_[scen] = n2ndVarTemp;	  
      nSecondStageCons_[scen] = asl_i[scen]->i.n_con_;

	  localData[scen].initialize(nSecondStageVars_[scen],nSecondStageCons_[scen]);

	  // get 2nd var ub lb
	  for(it=LocLocVarMap[scen].begin(); it!=LocLocVarMap[scen].end(); it++){
		localData[scen].collb[it->second] = asl_i[scen]->i.LUv_[it->first];
		localData[scen].colub[it->second] = asl_i[scen]->i.Uvx_[it->first];
		localData[scen].objGrad[it->second] = 0;
	  }

	  // get 2nd con ub lb 
	  for(int j=0; j<asl_i[scen]->i.n_con_; j++){
		localData[scen].rowlb[j] = asl_i[scen]->i.LUrhs_[j];
		localData[scen].rowub[j] = asl_i[scen]->i.Urhsx_[j];
	  }
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

  splitMatrices(0);

  //********************************************************************//
  nnzEqJac1st=0; 
  nnzIneqJac1st=0;
  amplRowMap1st = new int[nFirstStageCons_];
  assert(nFirstStageCons_==0);
  
  nnzEqJac2nd.resize(nScenarios_);
  nnzIneqJac2nd.resize(nScenarios_);
  amplRowMap2nd.resize(nScenarios_);

  getRowMap(0);

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

  getJacGoffMap(0);

  // set the root nl = first scenario nl
  firstStageData.locASL = asl_i[0];

  if(gNP_Alg!=0){
	doNetworkPart(0,amplSuffix);	
  }

  dimsEqual = true;
  for (int scen = 1; scen < nScenarios_; scen++) {
	nSecondStageCons_[scen] = 0;
	nSecondStageVars_[scen] = 0;
  }

  delete [] amplRowMap1st;
  delete amplSuffix;
  delete [] filedataArgv[0];
  delete [] filedataArgv[1];
}


// split matrices to Amat and Wmat
void amplGenStochInput::splitMatrices(const int scen)
{

  ASL_pfgh * asl;
  int amplNz=0,	Wnnz  = 0, Tnnz  = 0, Annz  = 0, 
  	  			QWnnz = 0, QTnnz = 0, QAnnzFromScen=0;
  int tempNz = 0;

  vector<int> starts2ndSt, starts1stStTemp, startsLink;
  vector<int> rowIdx2ndSt, rowIdx1stSt, rowIdxLink;
  vector<double> elts2ndSt, elts1stSt, eltsLink;

  vector<int> amplGoffRowIdx;
  vector<double> amplGoffColIdx;
  vector<int> amplGoffColStarts;
  vector<int> amplGoffColLength;
  vector<double> amplGoffElts;

  cgrad *cg;
  
  map<int,int>::iterator itVar, itVar_Row; 
  int wrk1stGoff=0, wrk2ndGoff=0, wrkLinkGoff=0;
  int wrkGoff=0;

////////////////////////////	 build Amat   ////////////////////////////////////////////////////////
  // set empty Amat. This is the matrix difined by 1st stage cons --- assume this is empty
  assert(nFirstStageCons_==0);
  starts1stStTemp.resize(nFirstStageVars_+1);
  for (int j = 0; j <= nFirstStageVars_; j++) {
	starts1stStTemp[j] = 0;
  }
  Amat.copyOf(true, nFirstStageCons_, nFirstStageVars_, 0, 0, 0, &starts1stStTemp[0], 0);

  asl = localData[scen].locASL;

  // set sparse structure of 2nd stage Jac and linking Jac
	Wnnz  = 0; Tnnz  = 0; 
  	QWnnz = 0;  QTnnz = 0; QAnnzFromScen=0;
	wrk1stGoff=0; wrk2ndGoff=0; wrkLinkGoff=0;
	wrkGoff=0;
	
////////////////////////////	build Wmat   and Tmat 	////////////////////////////////////////////////////////
	amplNz = nzc;

	amplGoffRowIdx.clear(); amplGoffRowIdx.resize(amplNz,0);
	amplGoffColIdx.clear(); amplGoffColIdx.resize(amplNz,0);
	amplGoffElts.clear(); amplGoffElts.resize(amplNz,0.);
	amplGoffColStarts.clear(); amplGoffColStarts.resize(n_var+1,0);
	amplGoffColLength.clear(); amplGoffColLength.resize(n_var,0);

	// get info from ampl
	for( int j=0; j < n_con; j++){
	  for(cg = Cgrad[j]; cg; cg = cg->next){
	    amplGoffRowIdx[cg->goff] = j;
		amplGoffColIdx[cg->goff] = cg->varno;
		amplGoffElts[cg->goff] = cg->coef;
		amplGoffColLength[cg->varno]++;

		itVar = LocGloVarMap[scen].find(cg->varno);
	    if(itVar != LocGloVarMap[scen].end()){
		  Tnnz++;
		}else{
		  Wnnz++;
		  itVar = LocLocVarMap[scen].find(cg->varno);
	      if(itVar == LocLocVarMap[scen].end())
		  	assert("impossible!" && 0);
		}
	  }
	}	
	for(int j=1; j < n_var+1; j++){
	  amplGoffColStarts[j] += amplGoffColStarts[j-1] + amplGoffColLength[j-1];
	}
	assert(amplGoffColStarts[n_var]==amplNz);
	assert(Wnnz+Tnnz==amplNz);
	assert(nSecondStageVars_[scen] + nFirstStageVars_==n_var);

	starts2ndSt.clear(); starts2ndSt.resize(nSecondStageVars_[scen]+1,0);
	rowIdx2ndSt.clear(); rowIdx2ndSt.resize(Wnnz);
	elts2ndSt.clear();   elts2ndSt.resize(Wnnz);

	startsLink.clear(); startsLink.resize(nFirstStageVars_+1,0);
	rowIdxLink.clear(); rowIdxLink.resize(Tnnz);
	eltsLink.clear();  eltsLink.resize(Tnnz);

  	for(int jcol=0; jcol < n_var; jcol++){
	  itVar = LocGloVarMap[scen].find(jcol);
	  if(itVar != LocGloVarMap[scen].end()){
	  	// this goff belongs to Tmat	
		startsLink[itVar->second+1] = amplGoffColLength[jcol];
	  }else{
		// this goff belongs to Wmat
		itVar = LocLocVarMap[scen].find(jcol);
		starts2ndSt[itVar->second+1] = amplGoffColLength[jcol];
	  }
	}
	for(int j=1; j < nFirstStageVars_+1; j++){
	  startsLink[j] += startsLink[j-1];
	}
	for(int j=1; j < nSecondStageVars_[scen]+1; j++){
	  starts2ndSt[j] += starts2ndSt[j-1];
	}

	assert(startsLink[nFirstStageVars_]  == Tnnz); 
	assert(starts2ndSt[nSecondStageVars_[scen]] == Wnnz);

	std::vector<int> nextJacDiagInCol; nextJacDiagInCol.resize(nSecondStageVars_[scen],0);
	for(int jj=0;jj<nSecondStageVars_[scen];jj++) nextJacDiagInCol[jj]=starts2ndSt[jj];	
	std::vector<int> nextJacLinkInCol; nextJacLinkInCol.resize(nFirstStageVars_,0);
	for(int jj=0;jj<nFirstStageVars_;jj++) nextJacLinkInCol[jj]=startsLink[jj];


	for(int jcol=0; jcol < n_var; jcol++){
	  itVar = LocGloVarMap[scen].find(jcol);
	  if(itVar != LocGloVarMap[scen].end()){
	  	// this goff belongs to Tmat
  	    for(int k=amplGoffColStarts[jcol]; k<amplGoffColStarts[jcol+1]; k++){
		  wrkGoff = nextJacLinkInCol[itVar->second]++;	
	  	  eltsLink[wrkGoff] = amplGoffElts[k];
		  rowIdxLink[wrkGoff] = amplGoffRowIdx[k];
		  LocTmatJacGoffMap[scen].insert( pair<int,int>(wrkGoff,k));	  
		  wrkLinkGoff++;
	    }	
	  }else{
		// this goff belongs to Wmat
		itVar = LocLocVarMap[scen].find(jcol); 
	    assert(itVar != LocLocVarMap[scen].end());		

	    for(int k=amplGoffColStarts[jcol]; k<amplGoffColStarts[jcol+1]; k++){
		  wrkGoff = nextJacDiagInCol[itVar->second]++;	
	  	  elts2ndSt[wrkGoff] = amplGoffElts[k];
		  rowIdx2ndSt[wrkGoff] = amplGoffRowIdx[k];
		  LocWmatJacGoffMap[scen].insert( pair<int,int>(wrkGoff,k)); 
		  wrk2ndGoff++;
	    }
	  }
	}
	assert(wrk2ndGoff  == Wnnz); 
	assert(wrkLinkGoff == Tnnz);
	for(int jj=0;jj<nSecondStageVars_[scen];jj++) 
	  assert(nextJacDiagInCol[jj]==starts2ndSt[jj+1]);	  
	for(int jj=0;jj<nFirstStageVars_;jj++) 
	  assert(nextJacLinkInCol[jj]==startsLink[jj+1]);

	assert(Wnnz!=0);

    Wmat[scen].copyOf(true, nSecondStageCons_[scen], nSecondStageVars_[scen], Wnnz, 
				&elts2ndSt[0], &rowIdx2ndSt[0], &starts2ndSt[0], 0);
	Tmat[scen].copyOf(true, nSecondStageCons_[scen], nFirstStageVars_, Tnnz, 
				&eltsLink[0], &rowIdxLink[0], &startsLink[0], 0);	




////////////////////////////	   build QAmat	 & QTmat & QWmat	////////////////////////////////////////

	// assume all the scenarios return same structure
	amplNz = localData[scen].Ampl_nnz_Hessian_Tri();

	amplGoffRowIdx.clear(); amplGoffRowIdx.resize(amplNz,0);
	amplGoffColIdx.clear(); amplGoffColIdx.resize(amplNz,0);
	amplGoffElts.clear(); amplGoffElts.resize(amplNz,0.);
	amplGoffColStarts.clear(); amplGoffColStarts.resize(n_var+1,0);
	amplGoffColLength.clear(); amplGoffColLength.resize(n_var,0);

	// realocate space
    starts2ndSt.clear(); starts2ndSt.resize(nSecondStageVars_[scen]+1,0);	
	startsLink.clear();  startsLink.resize(nFirstStageVars_+1,0);
    if(scen == 0){
      starts1stSt.clear(); starts1stSt.resize(nFirstStageVars_+1,0);
    }

	for(int QcolIdx=0; QcolIdx<n_var;QcolIdx++){
  	  //count 1st, 2nd and link size
	  for( int k=sputinfo->hcolstarts[QcolIdx];k<sputinfo->hcolstarts[QcolIdx+1];k++){	  
		int QrowIDX =  sputinfo->hrownos[k];//amplGoffRowIdx[k];
	    itVar = LocGloVarMap[scen].find(QcolIdx);		
		itVar_Row = LocGloVarMap[scen].find(QrowIDX); 
		
	  	if(itVar != LocGloVarMap[scen].end() && itVar_Row != LocGloVarMap[scen].end()){
		  //belong to QA
		  QAnnzFromScen++;
		  if(scen == 0){
		  	if(itVar_Row->second >= itVar->second)
		      starts1stSt[itVar->second+1]++;
			else
			  starts1stSt[itVar_Row->second+1]++;
		  }
		}else if(itVar == LocGloVarMap[scen].end() && itVar_Row == LocGloVarMap[scen].end()){
	      //belong to QW
		  QWnnz++;
		  itVar = LocLocVarMap[scen].find(QcolIdx);
		  itVar_Row = LocLocVarMap[scen].find(QrowIDX);
		  if(itVar_Row->second >= itVar->second)
		    starts2ndSt[itVar->second+1]++;
		  else
			starts2ndSt[itVar_Row->second+1]++;	  
		}else if(itVar != LocGloVarMap[scen].end() && itVar_Row == LocGloVarMap[scen].end()){
		  //belong to QT
		  QTnnz++;	  
		  startsLink[itVar->second+1]++;	
	    }else if(itVar == LocGloVarMap[scen].end() && itVar_Row != LocGloVarMap[scen].end()){
		  //belong to QT
		  QTnnz++;
		  startsLink[itVar_Row->second+1]++;	
		}else{assert("impossible"&&0);}
	  }
	}
	assert(QWnnz+QTnnz+QAnnzFromScen==amplNz);

	if(scen == 0){
	  for(int j=1; j < nFirstStageVars_+1; j++){
	    starts1stSt[j] += starts1stSt[j-1];
      }
	  assert(starts1stSt[nFirstStageVars_]==QAnnzFromScen);
	}

	for(int j=1; j < nSecondStageVars_[scen]+1; j++){
  		starts2ndSt[j] += starts2ndSt[j-1];
	}

    for(int j=1; j < nFirstStageVars_+1; j++){
	  startsLink[j] += startsLink[j-1];
    }
	
	assert(starts2ndSt[nSecondStageVars_[scen]]==QWnnz);
	assert(startsLink[nFirstStageVars_]==QTnnz);


	// realocate space
    rowIdx2ndSt.clear(); rowIdx2ndSt.resize(QWnnz);
	vector<int> colIdx2ndSt; colIdx2ndSt.resize(QWnnz);
    elts2ndSt.clear();   elts2ndSt.resize(QWnnz);

    rowIdxLink.clear();  rowIdxLink.resize(QTnnz);
    eltsLink.clear();    eltsLink.resize(QTnnz);

    if(scen == 0){
	  nnzQ1st = QAnnzFromScen;
    }
	assert(nnzQ1st == QAnnzFromScen);
	rowIdx1stSt.clear(); rowIdx1stSt.resize(QAnnzFromScen);
	elts1stSt.clear();	 elts1stSt.resize(QAnnzFromScen);

////////////////////////////		 build QWmat	////////////////////////////////////////////////////////

	wrk1stGoff=0; wrk2ndGoff=0; wrkLinkGoff=0; 
    wrkGoff=0;
	
	std::vector<int> nextQ1stInCol; nextQ1stInCol.resize(nFirstStageVars_,0);
	for(int jj=0;jj<nFirstStageVars_;jj++) nextQ1stInCol[jj]=starts1stSt[jj];
	std::vector<int> nextQ2ndInCol; nextQ2ndInCol.resize(nSecondStageVars_[scen],0);
	for(int jj=0;jj<nSecondStageVars_[scen];jj++) nextQ2ndInCol[jj]=starts2ndSt[jj];	
	std::vector<int> nextQCroInCol; nextQCroInCol.resize(nFirstStageVars_,0);
	for(int jj=0;jj<nFirstStageVars_;jj++) nextQCroInCol[jj]=startsLink[jj];

	int *LocQWmatHesGoffMap_ColWise = (int*)malloc(QWnnz*sizeof(int));
	int *LocQTmatHesGoffMap_ColWise = (int*)malloc(QTnnz*sizeof(int));
	int *LocQAmatHesGoffMap_ColWise = (int*)malloc(QAnnzFromScen*sizeof(int));
	
	for(int QcolIdx=0; QcolIdx < n_var;QcolIdx++){
  	  //count 1st, 2nd and link size	  
	  for( int k=sputinfo->hcolstarts[QcolIdx];k<sputinfo->hcolstarts[QcolIdx+1];k++){	  
		int QrowIDX = sputinfo->hrownos[k];// amplGoffRowIdx[k];
		itVar = LocGloVarMap[scen].find(QcolIdx);  
		itVar_Row = LocGloVarMap[scen].find(QrowIDX); 
		
	  	if(itVar != LocGloVarMap[scen].end() && itVar_Row != LocGloVarMap[scen].end()){
		  //belong to QA
//		  if(scen == 0){
			if(itVar_Row->second >= itVar->second){
			  wrkGoff = nextQ1stInCol[itVar->second]++;
			  rowIdx1stSt[wrkGoff] = itVar_Row->second;
			}
			else{
			  wrkGoff = nextQ1stInCol[itVar_Row->second]++;
			  rowIdx1stSt[wrkGoff] = itVar->second;
			}
			LocQAmatHesGoffMap_ColWise[wrkGoff] = k;
//		  }
		  wrk1stGoff++;
		}else if(itVar == LocGloVarMap[scen].end() && itVar_Row == LocGloVarMap[scen].end()){
	      //belong to QW
		  itVar = LocLocVarMap[scen].find(QcolIdx); 
		  itVar_Row = LocLocVarMap[scen].find(QrowIDX);
		  if(itVar_Row->second >= itVar->second){
			wrkGoff = nextQ2ndInCol[itVar->second]++;
			rowIdx2ndSt[wrkGoff] = itVar_Row->second;
			colIdx2ndSt[wrkGoff] = itVar->second;
		  }
		  else{
		  	wrkGoff = nextQ2ndInCol[itVar_Row->second]++;
			rowIdx2ndSt[wrkGoff] = itVar->second;
			colIdx2ndSt[wrkGoff] = itVar_Row->second;
		  }
		  LocQWmatHesGoffMap_ColWise[wrkGoff] = k;		  
		  wrk2ndGoff++;		  
		}else if(itVar != LocGloVarMap[scen].end() && itVar_Row == LocGloVarMap[scen].end()){
		  //belong to QT
		  itVar_Row = LocLocVarMap[scen].find(QrowIDX);
		  wrkGoff = nextQCroInCol[itVar->second]++;
		  rowIdxLink[wrkGoff] = itVar_Row->second;	
		  LocQTmatHesGoffMap_ColWise[wrkGoff] = k;
		  wrkLinkGoff++;		  
	    }else if(itVar == LocGloVarMap[scen].end() && itVar_Row != LocGloVarMap[scen].end()){
		  //belong to QT	 
		  itVar = LocLocVarMap[scen].find(QcolIdx);	
		  wrkGoff = nextQCroInCol[itVar_Row->second]++;
		  rowIdxLink[wrkGoff] = itVar->second;			
		  LocQTmatHesGoffMap_ColWise[wrkGoff] = k;
		  wrkLinkGoff++;	
		}else{assert("impossible"&&0);}
	  }
	}
	for(int jj=0;jj<nFirstStageVars_;jj++) 
	  assert(nextQ1stInCol[jj]==starts1stSt[jj+1]);
	for(int jj=0;jj<nSecondStageVars_[scen];jj++) 
	  assert(nextQ2ndInCol[jj]==starts2ndSt[jj+1]);	  
	for(int jj=0;jj<nFirstStageVars_;jj++) {
	  assert(nextQCroInCol[jj]==startsLink[jj+1]);
	}
	assert(wrk1stGoff == QAnnzFromScen &&  wrk2ndGoff == QWnnz && wrkLinkGoff ==QTnnz);

	nnzQ2nd[scen] = QWnnz;
	nnzQCross2nd[scen] = QTnnz;
	
	if(scen == 0){
	  QAmat.copyOf(true, nFirstStageVars_, nFirstStageVars_, QAnnzFromScen,
					  &elts1stSt[0], &rowIdx1stSt[0], &starts1stSt[0], 0);
	}
	QWmat[scen].copyOf(true, nSecondStageVars_[scen], nSecondStageVars_[scen], QWnnz,
					  &elts2ndSt[0], &rowIdx2ndSt[0], &starts2ndSt[0], 0);
	QTmat[scen].copyOf(true, nSecondStageVars_[scen], nFirstStageVars_, QTnnz,
					  &eltsLink[0], &rowIdxLink[0], &startsLink[0], 0);

//////////////////////////////////////////////////////////////////////////////////////////
	// compute GoffMap in row_wise

	std::vector<int> starts_Rowbeg, nextIdxInRowWise; 

	int wrkGoff_row=0;
	wrk1stGoff=0; wrk2ndGoff=0; wrkLinkGoff=0; 

	//////////////////////////////         for QTmat     ///////////////////////////////////////
	// get rowwise for QTmat
	starts_Rowbeg.clear(); starts_Rowbeg.resize(nSecondStageVars_[scen]+1,0);
	for(int k=0;k<QTnnz;k++){
	  starts_Rowbeg[rowIdxLink[k]+1]++;
	}
	for(int jj=1;jj<nSecondStageVars_[scen]+1;jj++) starts_Rowbeg[jj]+=starts_Rowbeg[jj-1];	

	nextIdxInRowWise.clear(); nextIdxInRowWise.resize(nSecondStageVars_[scen]+1,0);
	for(int jj=0;jj<nSecondStageVars_[scen]+1;jj++) nextIdxInRowWise[jj]=starts_Rowbeg[jj];

	// build LocQTmatHesGoffMap    
	for(int locColID=0;locColID<nFirstStageVars_;locColID++) {
	  for(int k=startsLink[locColID];k<startsLink[locColID+1];k++){
		int QrowIDX = rowIdxLink[k];
		wrkGoff_row = nextIdxInRowWise[QrowIDX]++;
		LocQTmatHesGoffMap[scen].insert( pair<int,int>(LocQTmatHesGoffMap_ColWise[k],wrkGoff_row));
		wrkLinkGoff++;
	  }
	}
	assert(wrkLinkGoff ==QTnnz);
	
	for(int jj=0;jj<nSecondStageVars_[scen];jj++) {
	  assert(nextIdxInRowWise[jj]==starts_Rowbeg[jj+1]);
	}
	
	//////////////////////////////         for QWmat     ///////////////////////////////////////
  	// get rowwise for QWmat
	starts_Rowbeg.clear();starts_Rowbeg.resize(nSecondStageVars_[scen]+1,0);

	for(int k=0;k<QWnnz;k++){
	  int rowID = rowIdx2ndSt[k];
	  int colID = colIdx2ndSt[k];
	  assert(rowID>=colID);
	  starts_Rowbeg[rowID+1]++;
	}

	for(int jj=1;jj<nSecondStageVars_[scen]+1;jj++) starts_Rowbeg[jj]+=starts_Rowbeg[jj-1];	

	nextIdxInRowWise.clear(); nextIdxInRowWise.resize(nSecondStageVars_[scen]+1,0);
	for(int jj=0;jj<nSecondStageVars_[scen]+1;jj++) nextIdxInRowWise[jj]=starts_Rowbeg[jj];

	//build LocQWmatHesGoffMap
	for(int locColID=0; locColID < nSecondStageVars_[scen];locColID++) {
	  for(int k=starts2ndSt[locColID];k<starts2ndSt[locColID+1];k++){
		int QrowIDX = rowIdx2ndSt[k];
		wrkGoff_row = nextIdxInRowWise[QrowIDX]++;
		LocQWmatHesGoffMap[scen].insert( pair<int,int>(LocQWmatHesGoffMap_ColWise[k],wrkGoff_row));
		wrk2ndGoff++;
	  }
	}
	assert(wrk2ndGoff ==QWnnz);

	for(int jj=0;jj<nSecondStageVars_[scen];jj++) {
	  assert(nextIdxInRowWise[jj]==starts_Rowbeg[jj+1]);
	}

	//////////////////////////////         for QAmat     ///////////////////////////////////////
  	// get rowwise for QAmat
	starts_Rowbeg.clear();starts_Rowbeg.resize(nFirstStageVars_+1,0);
	for(int k=0;k<QAnnzFromScen;k++){
	  int rowID = rowIdx1stSt[k];
	  int colID=-1; 
	  for(int jc=0;jc<nFirstStageVars_+1;jc++){
	  	if(starts1stSt[jc]<=k && k<starts1stSt[jc+1]){colID=jc; break;}
	  }
	  assert(colID!=-1);

	  if(rowID>=colID)
	  	starts_Rowbeg[rowID+1]++;
	  else
	  	starts_Rowbeg[colID+1]++;
	}
	for(int jj=1;jj<nFirstStageVars_+1;jj++) starts_Rowbeg[jj]+=starts_Rowbeg[jj-1];	


	nextIdxInRowWise.clear(); nextIdxInRowWise.resize(nFirstStageVars_+1,0);
	for(int jj=0;jj<nFirstStageVars_+1;jj++) nextIdxInRowWise[jj]=starts_Rowbeg[jj];

	//build LocQAmatHesGoffMap
	for(int locColID=0; locColID < nFirstStageVars_;locColID++) {
	  for(int k=starts1stSt[locColID];k<starts1stSt[locColID+1];k++){
		int QrowIDX = rowIdx1stSt[k];
		wrkGoff_row = nextIdxInRowWise[QrowIDX]++;
		LocQAmatHesGoffMap[scen].insert( pair<int,int>(LocQAmatHesGoffMap_ColWise[k],wrkGoff_row));
		wrk1stGoff++;
	  }
	}
	assert(wrk1stGoff == QAnnzFromScen);
	
	for(int jj=0;jj<nFirstStageVars_;jj++) {
	  assert(nextIdxInRowWise[jj]==starts_Rowbeg[jj+1]);
	}

	free(LocQTmatHesGoffMap_ColWise);
    free(LocQWmatHesGoffMap_ColWise);
	free(LocQAmatHesGoffMap_ColWise);

	FreeAll(starts2ndSt); FreeAll(starts1stStTemp); FreeAll(startsLink);
	FreeAll(rowIdx2ndSt); FreeAll(rowIdx1stSt); FreeAll(rowIdxLink);
	FreeAll(elts2ndSt); FreeAll(elts1stSt); FreeAll(eltsLink);

	FreeAll(amplGoffRowIdx);
	FreeAll(amplGoffColIdx);
	FreeAll(amplGoffColStarts);
	FreeAll(amplGoffColLength);
	FreeAll(amplGoffElts);
	  
}


void amplGenStochInput::getRowMap(const int scen)
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



void amplGenStochInput::getJacGoffMap(const int scen)
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

  // this is row-wise, correct!
  int T_nvar = Tmat[scen].getMajorDim();
  int W_nvar = Wmat[scen].getMajorDim();
  assert(nFirstStageVars_==T_nvar);
  assert(nSecondStageVars_[scen]==W_nvar);

  int Tnnz = Tmat[scen].getNumElements();
  int Wnnz = Wmat[scen].getNumElements();

  const int* tColstart 	= Tmat[scen].getVectorStarts();  
  const int* tRowIDX 	= Tmat[scen].getIndices();   

  const int* wColstart 	= Wmat[scen].getVectorStarts();  
  const int* wRowIDX 	= Wmat[scen].getIndices();    

  int rowID;
  for( int k = 0; k < Tnnz; k++ ) {
	rowID = tRowIDX[k];
	if( amplRowMap2nd[scen][rowID] < 0){
	  nnzALink2nd[scen]++;
	}else{
	  nnzCLink2nd[scen]++;
	}
  }
  for( int k = 0; k < Wnnz; k++ ) {
    rowID = wRowIDX[k];
	if( amplRowMap2nd[scen][rowID] < 0){
	  nnzALoc2nd[scen]++;
	}else{
	  nnzCLoc2nd[scen]++;
	}
  }
  assert( amplNz == nnzALink2nd[scen] + nnzALoc2nd[scen]+nnzCLink2nd[scen] + nnzCLoc2nd[scen]);


  // for T 
  int *irowT_A 	= (int*)malloc(nnzALink2nd[scen]*sizeof(int));
  int *jcolT_A 	= (int*)malloc(nnzALink2nd[scen]*sizeof(int));
  int *newTGoff_A = (int*)malloc(nnzALink2nd[scen]*sizeof(int));

  int *irowT_C 	= (int*)malloc(nnzCLink2nd[scen]*sizeof(int));
  int *jcolT_C 	= (int*)malloc(nnzCLink2nd[scen]*sizeof(int));
  int *newTGoff_C = (int*)malloc(nnzCLink2nd[scen]*sizeof(int));

  int kA=0, kC=0;
  for( int j = 0; j < T_nvar; j++ ) {
    for( int k = tColstart[j]; k < tColstart[j+1]; k++ ) {
	  rowID = tRowIDX[k];
	  if( amplRowMap2nd[scen][rowID] < 0){
	    irowT_A[kA]=rowID; jcolT_A[kA] = j; newTGoff_A[kA]=k;
		kA++;
	  }else{
	    irowT_C[kC]=rowID; jcolT_C[kC] = j; newTGoff_C[kC]=k;
		kC++;
	  }
    }
  }
  RowColWise_map(irowT_A,kA,jcolT_A,newTGoff_A);
  RowColWise_map(irowT_C,kC,jcolT_C,newTGoff_C);  

  for(int k = 0; k < kA; k++){
	JacALinkGoff2nd[scen].insert( pair<int,int>(LocTmatJacGoffMap[scen][newTGoff_A[k]],k));	
  }
  for(int k = 0; k < kC; k++){
	JacCLinkGoff2nd[scen].insert( pair<int,int>(LocTmatJacGoffMap[scen][newTGoff_C[k]],k));	
  }

  free(irowT_A); free(jcolT_A); free(newTGoff_A);
  free(irowT_C); free(jcolT_C); free(newTGoff_C); 




  // for W
  int *irowW_A 	= (int*)malloc(nnzALoc2nd[scen]*sizeof(int));
  int *jcolW_A 	= (int*)malloc(nnzALoc2nd[scen]*sizeof(int));
  int *newWGoff_A = (int*)malloc(nnzALoc2nd[scen]*sizeof(int));

  int *irowW_C 	= (int*)malloc(nnzCLoc2nd[scen]*sizeof(int));
  int *jcolW_C 	= (int*)malloc(nnzCLoc2nd[scen]*sizeof(int));
  int *newWGoff_C = (int*)malloc(nnzCLoc2nd[scen]*sizeof(int));
  
  kA=0; kC=0;
  for( int j = 0; j < W_nvar; j++ ) {
    for( int k = wColstart[j]; k < wColstart[j+1]; k++ ) {
	  rowID = wRowIDX[k];
	  if( amplRowMap2nd[scen][rowID] < 0){
	    irowW_A[kA]=rowID; jcolW_A[kA] = j; newWGoff_A[kA]=k;
		kA++;
	  }else{
	    irowW_C[kC]=rowID; jcolW_C[kC] = j; newWGoff_C[kC]=k;
		kC++;
	  }
    }
  }  
  RowColWise_map(irowW_A,kA,jcolW_A,newWGoff_A);
  RowColWise_map(irowW_C,kC,jcolW_C,newWGoff_C);    


  for(int k = 0; k < kA; k++){
	JacALocGoff2nd[scen].insert( pair<int,int>(LocWmatJacGoffMap[scen][newWGoff_A[k]],k));	
  }
  for(int k = 0; k < kC; k++){
	JacCLocGoff2nd[scen].insert( pair<int,int>(LocWmatJacGoffMap[scen][newWGoff_C[k]],k));	
  }
  free(irowW_A); free(jcolW_A); free(newWGoff_A);
  free(irowW_C); free(jcolW_C); free(newWGoff_C);

  assert( amplNz == nnzALink2nd[scen] + nnzALoc2nd[scen]+nnzCLink2nd[scen] + nnzCLoc2nd[scen]);

}

void amplGenStochInput::		loadLocalNLdata(int scen)
{
  
  if (localData[scen].didLoad) return;
  localData[scen].didLoad=true;

  // Add the suffix  
  AmplSuffix *amplSuffix = new AmplSuffix();
  amplSuffix->DefineSuffix("pipsNLP_1stStageVar_in", Suffix_Var, Suffix_Int);  
  if(gUseReducedSpace>0){  
	  amplSuffix->DefineSuffix("pipsNLP_DecisionVar_in", Suffix_Var, Suffix_Int);   
  }
  if(gNP_Alg!=0){
  	amplSuffix->DefineSuffix("pipsNLP_VarPartIdx_in", Suffix_Var, Suffix_Int); 
  	amplSuffix->DefineSuffix("pipsNLP_ConPartIdx_in", Suffix_Con, Suffix_Int); 
  }

  stringstream fname;
  char *filedataArgv[3];
  filedataArgv[0] = new char[9];
  memcpy(filedataArgv[0],"pips_nlp",9);	  
  unsigned filelen;
  map<int,int>::iterator it;
  std::string strTemp;

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
  filedataArgv[1] = new char[filelen];
  memcpy(filedataArgv[1],fname.str().c_str(),filelen);
  filedataArgv[2]=NULL;
 
  // Allocate the ampl solver library (ASL) context.
  cout << "MPI#" << mype_ << ":  Start read pfgh_read! " << endl;  
  asl_i[scen] = localData[scen].initASL(filedataArgv, amplSuffix); 
  cout << "MPI#" << mype_ << ":  Scenario " << scen << "   pfgh_read Done!" << endl; 

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
	  //creat local NL to local var map
	  LocLocVarMap[scen].insert( pair<int,int>(j,n2ndVarTemp));
	  n2ndVarTemp++;
	}
  }

  assert(nFirstStageVars_ == n1stVarTemp);
  nSecondStageVars_[scen] = n2ndVarTemp;
  nSecondStageCons_[scen] = asl_i[scen]->i.n_con_;  
  if(dimsEqual==true){
  	if(nSecondStageVars_[scen]!=nSecondStageVars_[0] || nSecondStageCons_[scen]!=nSecondStageCons_[0]) dimsEqual=false;
  }
  
  // define decision var
  if(gUseReducedSpace>0){	
	decisionVarDim[scen]=0;
	int *decisionVarIDX = amplSuffix->GetSuffixVal_Int(asl_i[scen], "pipsNLP_DecisionVar_in", Suffix_Var);
	//find correct index in C:	in ampl model, index of decision var starts from 1, here we correct it as zero;
	for(int j=0; j<asl_i[scen]->i.n_var_; j++){
	  decisionVarIDX[j] -=1;
	  if(decisionVarIDX[j]>=0) 
		decisionVarDim[scen]++;
	}
	schurVarConIDinNL[scen] = (int*)malloc(decisionVarDim[scen]*sizeof(int));
	  
  	for(int findSCVar=0; findSCVar<decisionVarDim[scen];findSCVar++){
	  for(int j=0; j<asl_i[scen]->i.n_var_; j++){
		if(decisionVarIDX[j] == findSCVar){
		  map<int,int>::iterator itVarTemp = LocLocVarMap[scen].find(j);	 
		  assert(itVarTemp != LocLocVarMap[scen].end());
		  schurVarConIDinNL[scen][findSCVar]=itVarTemp->second;
		  break;
		}
	  }
	}
  }

  localData[scen].initialize(nSecondStageVars_[scen],nSecondStageCons_[scen]);
	  
  // get 2nd var ub lb
  for(it=LocLocVarMap[scen].begin(); it!=LocLocVarMap[scen].end(); it++){
	localData[scen].collb[it->second] = asl_i[scen]->i.LUv_[it->first];
	localData[scen].colub[it->second] = asl_i[scen]->i.Uvx_[it->first];
	localData[scen].objGrad[it->second] = 0;
  }

  // get 2nd con ub lb 
  for(int j=0; j<asl_i[scen]->i.n_con_; j++){
	localData[scen].rowlb[j] = asl_i[scen]->i.LUrhs_[j];
	localData[scen].rowub[j] = asl_i[scen]->i.Urhsx_[j];
  }

  splitMatrices(scen);
  getRowMap(scen);
  getJacGoffMap(scen);

  if(gNP_Alg!=0){ 	
	doNetworkPart(scen,amplSuffix);	
  }

  delete amplSuffix;
  delete [] filedataArgv[1];
  
}


void amplGenStochInput::doNetworkPart(const int scen, AmplSuffix* amplSuffix )
{
  ASL_pfgh * asl;
  asl = localData[scen].locASL;
  int nvar = n_var, ncon = n_con;

  NR_partIDX_var.resize(nScenarios_,NULL);
  NR_partIDX_con.resize(nScenarios_,NULL);

  int nb_primal=0, nb_dual=0;

  int *temp_varID = amplSuffix->GetSuffixVal_Int(asl_i[scen], "pipsNLP_VarPartIdx_in", Suffix_Var);
  int *temp_conID = amplSuffix->GetSuffixVal_Int(asl_i[scen], "pipsNLP_ConPartIdx_in", Suffix_Con);	

  int loc_nx = nSecondStageVars_[scen];
  int loc_my=0, loc_mz=0;
  
  for(int j=0; j<ncon;j++){
	if( amplRowMap2nd[scen][j] < 0){
	  loc_my++;
	}else{
	  loc_mz++;
	}
  }

  nb_primal = loc_nx+loc_mz;
  nb_dual = nSecondStageCons_[scen];

  assert(nb_dual==loc_mz+loc_my);


  NR_partIDX_var[scen] = (int*)malloc(nb_primal*sizeof(int));
  NR_partIDX_con[scen] = (int*)malloc(nb_dual*sizeof(int));

  int findID=0, findConID=0;
  for(int i=0;i<nvar;i++){
	if(temp_varID[i]>=0){
	  NR_partIDX_var[scen][findID++]=temp_varID[i];
	  if(findID==loc_nx) break;
	}
  }

  int conID=0;
  for(int i=0;i<ncon;i++){
  	if( amplRowMap2nd[scen][i] < 0){
	  conID = -amplRowMap2nd[scen][i] - 1;
	  NR_partIDX_con[scen][conID]=temp_conID[i];
	  findConID++;
	}else{
	  conID = amplRowMap2nd[scen][i];
	  NR_partIDX_var[scen][conID + loc_nx]=temp_conID[i];
	  findID++;
	  NR_partIDX_con[scen][conID + loc_my]=temp_conID[i];
	  findConID++;
	}
  }
  assert(findConID==nb_dual && findID==nb_primal);
  
}


amplGenStochInput::~amplGenStochInput(){
	for( int i = 0; i < nScenarios_; i++ ) {
		if(amplRowMap2nd[i])
			delete [] amplRowMap2nd[i];

	}	
}
