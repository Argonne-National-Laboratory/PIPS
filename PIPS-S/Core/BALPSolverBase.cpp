#include "BALPSolverBase.hpp"

#include "PIPSLogging.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include "CoinTime.hpp"

using namespace std;

BALPSolverBase::BALPSolverBase(const BAData &data) : nIter(0), startTime(0.0), data(data), status(Uninitialized),
	dualTol(1.0e-7), primalTol(1.0e-7), zeroTol(1.0e-12), primalRelTol(1.0e-9), pivotTol(5.0e-7),
	phase1(false), doreport(false), reportFrequency(10000000), dumpEvery(0),
	replaceFirst(0), firstReplaceSecond(0), secondReplaceFirst(0), replaceSecondSelf(0), 
	replaceSecondOther(0), wasOptimal(0), logIters(0),
	ftranTime(0.0), btranTime(0.0), updateIteratesTime(0.0), ftranDSETime(0.0), invertTime(0.0),
		selectEnteringTime(0.0), selectLeavingTime(0.0), priceTime(0.0), updateColumnTime(0.0)
	{

	x.allocate(data.dims,data.ctx, PrimalVector);
	d.allocate(data.dims,data.ctx, PrimalVector);

	basicIdx.allocate(data.dims,data.ctx, BasicVector);
	infeasList.allocate(data.dims,data.ctx, PrimalVector);
	
	rowVec.allocate(data.dims,data.ctx, PrimalVector);
	ftranVec.allocate(data.dims,data.ctx, BasicVector);
	ftranVec2.allocate(data.dims,data.ctx, BasicVector);
	btranVec.allocate(data.dims,data.ctx, BasicVector);

	states.allocate(data.dims,data.ctx, PrimalVector);
	
	// based on 6.26 in Koberstein thesis
	// TODO: this needs to be revised
	//reinvertFrequency = min(80+data.dims.totalCons()/300,1750);
	reinvertFrequency = 150;
	//reinvertFrequency = 1;

	

	la = new BAPFIParWrapper<BALinearAlgebra,BAPFIPar>(data);
	reinvertFrequency_backup = reinvertFrequency;
	lastBadReinvert = lastGoodReinvert = 0;

}
/*
void BALPSolverBase::initializeBounds() {
	
	// implicitly loops over first stage also
	const vector<int> &localscen = data.l.localScenarios();
	for (unsigned i = 0; i < localscen.size(); i++) {
		int scen = localscen[i];
		int nvar = data.dims.numVars(scen);
		denseVector &l2 = data.l.getVec(scen);
		denseVector &u2 = data.u.getVec(scen);
		denseVector &lTol2 = lTol.getVec(scen);
		denseVector &uTol2 = uTol.getVec(scen);
		for (int idx = 0; idx < nvar; idx++) {
			lTol2[idx] = l2[idx] - primalTol;
			uTol2[idx] = u2[idx] + primalTol;
			//lTol2[idx] = l2[idx]-abs(l2[idx])*primalRelTol-primalTol;
			//uTol2[idx] = u2[idx]+abs(u2[idx])*primalRelTol+primalTol;
		}
	}	

}*/

BALPSolverBase::~BALPSolverBase() {
	if (la) delete la;
	if (logIters) delete logIters;
}

void BALPSolverBase::setupIndices() {

	int nbasic = 0;
	basicIdx.clear();
	const vector<int> &localScen = basicIdx.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		vector<int> &basicIdx2 = basicIdx.getVec(scen);
		const denseFlagVector<variableState> &states2 = states.getVec(scen);
		int nvar = data.dims.numVars(scen);
		for (int j = 0; j < nvar; j++) {
			if (states2[j] == Basic) {
				//printf("%d %d %d %d\n", scen, j, nbasic,basicIdx2.size());
				basicIdx2.push_back(j);
				nbasic++;
			}
		}
	}

	int nbasic1 = basicIdx.getFirstStageVec().size();
	nbasic = data.ctx.reduce(nbasic-nbasic1) + nbasic1;
	//printf("nbasic: %d total cons: %d\n", nbasic, data.dims.totalCons());
	assert(nbasic == data.dims.totalCons());
	//if (data.ctx.mype() == 0) cout << "nbasic first stage: " << nbasic1 << endl;

}


// assume states is already set
// only call with !reinvert if flipping bounds
void BALPSolverBase::initialize(bool reinvert) {

	double t;	
	if (reinvert) { 
		rebuildIndices();
		t = MPI_Wtime();
		la->reinvert(states);
		invertTime += MPI_Wtime()-t;
		lastReinvert = nIter;
		if (logIters) {
			*logIters << CoinCpuTime() << " " << nIter << endl;
		}
	}

	const vector<int> &localScen = basicIdx.localScenarios();
	
	//FILE *f = fopen("../balpvars","w");
	
	// calculate xb = B^{-1}(b-An*xn)
	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		const denseFlagVector<constraintType> &vartype1 = data.vartype.getVec(scen);
		const denseFlagVector<variableState> &states1 = states.getVec(scen);
		denseVector &x1 = x.getVec(scen);
		const denseVector &l1 = myLB().getVec(scen);
		const denseVector &u1 = myUB().getVec(scen);
		int nvar = data.dims.numVars(scen);
		for (int i = 0; i < nvar; i++) {
			if (states1[i] == Basic) {
				x1[i] = 0; // needs to be zero for multiply
			} else if (vartype1[i] == Free) {
				x1[i] = 0;
			} else if (states1[i] == AtLower) {
				x1[i] = l1[i];
			} else {
				x1[i] = u1[i];
			}
			//if (states1[i] != Basic)fprintf(f,"%s %.12e\n", data.names.getVec(scen)[i].c_str(),x1[i]);
		}
	}
	sparseBAVector &xb = ftranVec;
	xb.clear();
	t = MPI_Wtime();
	data.multiply(x,xb,states);
	priceTime += MPI_Wtime() -t;

	/*for (int j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		CoinIndexedVector &xb1 = xb.getVec(scen).v;
		double *xbElts = xb1.denseVector();
		int ncons = data.dims.numCons(scen);
		int nvarReal = data.dims.inner.numVars(scen);
		for (int i = 0; i < ncons; i++) {
			fprintf(f,"%s %.12e\n", data.names.getVec(scen)[nvarReal+i].c_str(),xbElts[i]);
			
		}
	}*/

	xb.negate();
	t = MPI_Wtime();
	la->ftran(xb);
	ftranTime += MPI_Wtime() - t;

	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		denseVector &x1 = x.getVec(scen);
		CoinIndexedVector &xb1 = xb.getVec(scen).v;
		double *xb1Elts = xb1.denseVector();
		const vector<int> &basicIdx1 = basicIdx.getVec(scen);
		for (unsigned i = 0; i < basicIdx1.size(); i++) {
			x1[basicIdx1[i]] = xb1Elts[i];
		//fprintf(f,"%s %.12e\n", data.names.getVec(scen)[basicIdx1[i]].c_str(),xb1Elts[i]);
			
		}
	}
	//fclose(f);
	//assert(nIter == -1);

	// calculate y = B^{-T}c_B
	sparseBAVector &y = btranVec;
	y.clear();

	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		const denseVector &c1 = myObjective().getVec(scen);
		CoinIndexedVector &y1 = y.getVec(scen).v;
		const vector<int> &basicIdx1 = basicIdx.getVec(scen);
		for (unsigned i = 0; i < basicIdx1.size(); i++) {
			double coef = c1[basicIdx1[i]];
			if (coef) {
				y1.insert(i,coef);
			}
		}
	}
	t = MPI_Wtime();
	la->btran(y);
	btranTime += MPI_Wtime() - t;
	
	// calculate dn = cn - An^Ty
	sparseBAVector &dn = rowVec;
	dn.clear();
	t = MPI_Wtime();
	data.multiplyT(y,dn);
	priceTime += MPI_Wtime() - t;
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		denseVector &d1 = d.getVec(scen);
		CoinIndexedVector &dn1 = dn.getVec(scen).v;
		const denseVector &c1 = myObjective().getVec(scen);
		denseFlagVector<variableState> &states1 = states.getVec(scen);
		int nvar = data.dims.numVars(scen);
		for (int j = 0; j < nvar; j++) {
			if (states1[j] != Basic) {
				d1[j] = c1[j] - dn1[j];
			} else {
				d1[j] = 0.0;
			}
		}
	}

	double objvalNew = x.dotWith(myObjective());

	largestPrimalError = calculateLargestPrimalError();
	/*if (data.ctx.mype() == 0 && nIter >= 1 && objvalNew - objval < -1e-7) {
		printf("OBJECTIVE DECREASE: was %.9g now %.9g (Iteration %d)\n",objval,objvalNew,nIter);
	}*/
	
	objval = objvalNew;

	status = Initialized;

	initializeIndexedInfeasibilities();

	if (checkDualFeasible(infeasList)) {
		if (checkPrimalFeasible()) {
			status = Optimal;
		} else {
			status = DualFeasible;
		}
	} else if (checkPrimalFeasible()) {
		status = PrimalFeasible;
	}


}

void BALPSolverBase::setStates(const BAFlagVector<variableState> &state) {

	assert(states.getFirstStageVec().length() == data.dims.numFirstStageVars());
	const vector<int> &localScen = data.l.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		const denseFlagVector<variableState> &s = state.getVec(scen);
		denseFlagVector<variableState> &mys = states.getVec(scen);
		int n = data.dims.numVars(scen);
		assert(s.length() == n);
		for (int j = 0; j < n; j++) {
			mys[j] = s[j];
		}
	}
	status = LoadedFromFile;
	setupIndices();
}


bool BALPSolverBase::checkDualFeasible(BAContainer<vector<int> > &infeas) {
	assert(status != Uninitialized && status != LoadedFromFile);

	infeas.clear();

	int ninfeas = 0;

	const vector<int> &localScen = infeas.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		denseFlagVector<variableState> &states1 = states.getVec(scen);
		const denseFlagVector<constraintType> &vartype1 = data.vartype.getVec(scen);
		vector<int> &infeas1 = infeas.getVec(scen);
		const denseVector &d1 = d.getVec(scen);
		int nvar = data.dims.numVars(scen);
		for (int idx = 0; idx < nvar; idx++) {
			if (states1[idx] == Basic) continue;
			if (vartype1[idx] == Free && (d1[idx] < -dualTol || d1[idx] > dualTol)) {
				infeas1.push_back(idx);
				ninfeas++;
			} else if (states1[idx] == AtLower && d1[idx] < -dualTol) {
				if (vartype1[idx] == Fixed) states1[idx] = AtUpper;
				else { 
					infeas1.push_back(idx); ninfeas++; 
					//printf("lower infeas, %s %e\n",data.names.getVec(scen)[idx].c_str(),d1[idx]);
				}
			} else if (states1[idx] == AtUpper && d1[idx] > dualTol) {
				if (vartype1[idx] == Fixed) states1[idx] = AtLower;
				else { 
					infeas1.push_back(idx); ninfeas++; 
					//printf("upper infeas, %s %e\n",data.names.getVec(scen)[idx].c_str(),d1[idx]);
				}
			}
		}
	}

	int ninfeas1 = infeas.getFirstStageVec().size();
	ninfeas = data.ctx.reduce(ninfeas-ninfeas1) + ninfeas1;
	//if (data.ctx.mype() == 0 && ninfeas) cout << "number dual infeasible: " << ninfeas << endl;
	return (ninfeas == 0);
}



bool BALPSolverBase::checkPrimalFeasible() const {
	
	const vector<int> &localScen = x.localScenarios();
	int found = 0;

	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		const denseVector &x1 = x.getVec(scen);
		const denseVector &l1 = myLB().getVec(scen);
		const denseVector &u1 = myUB().getVec(scen);
		int nvar = data.dims.numVars(scen);
		for (int i = 0; i < nvar; i++) {
			if (x1[i] > u1[i] + primalTol) {
				//printf("(%d,%d) at %.7g above ub %.7g by %.7g\n",scen,i,x1[i],u1[i],x1[i]-u1[i]);
				found = 1; break;
			}
			if (x1[i] < l1[i] - primalTol) {
				//printf("(%d,%d) at %.7g below lb %.7g by %.7g\n",scen,i,x1[i],l1[i],l1[i]-x1[i]);
				found = 1; break;
			}
		}
		if (found) break;
	}
	found = data.ctx.reduce(found);


	return !found;
}


#include <sstream>
#include <iomanip>

void BALPSolverBase::writeBasis(const string &filebase) {
	assert(0);
	// append iterate number
	stringstream ss;
	ss << filebase << nIter;
	string s = ss.str();

	ofstream f(s.c_str());

	// not MPI-safe, will print 1st stage on each process	
	const vector<int> &localScen = basicIdx.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		const denseFlagVector<variableState> &states2 = states.getVec(scen);
		int nvar = data.dims.numVars(scen);
		for (int j = 0; j < nvar; j++) {
			if (states2[j] == Basic) {
				f << data.names.getVec(scen)[j] << endl;
			}
		}
	}


}

void BALPSolverBase::recordPivot(int enterScen, int leaveScen) {
	if (enterScen == -1) {
		if (leaveScen == -1) {
			replaceFirst++;
		} else {
			firstReplaceSecond++;
		}
	} else {
		if (leaveScen == enterScen) {
			replaceSecondSelf++;
		} else if (leaveScen == -1) {
			secondReplaceFirst++;
		} else {
			replaceSecondOther++;
		}
	}
}

double BALPSolverBase::calculateLargestPrimalError() const {
	const vector<int> &localScen = basicIdx.localScenarios();
	CoinBigIndex j;

	double largest = 0.0;
	int nvar1 = data.dims.inner.numFirstStageVars();
	int nrow1 = data.dims.numFirstStageCons();
	const double *Aelts = data.Arow->getElements();
	const int *Aidx = data.Arow->getIndices();
	const CoinBigIndex *Astart = data.Arow->getVectorStarts();
	const denseVector &x1 = x.getFirstStageVec();
	for (int i = 0; i < nrow1; i++) {
		double row = 0.0;
		for (j = Astart[i]; j < Astart[i+1]; j++) {
			row += Aelts[j]*x1[Aidx[j]];
		}
		row -= x1[nvar1+i];
		row = abs(row);
		if (row > largest) largest = row;
	}
		

	for (unsigned k = 1; k < localScen.size(); k++) {
		int scen = localScen[k];
		const double *Welts = data.Wrow[scen]->getElements();
		const int *Widx = data.Wrow[scen]->getIndices();
		const CoinBigIndex *Wstart = data.Wrow[scen]->getVectorStarts();
		const double *Telts = data.Trow[scen]->getElements();
		const int *Tidx = data.Trow[scen]->getIndices();
		const CoinBigIndex *Tstart = data.Trow[scen]->getVectorStarts();
		const denseVector &x2 = x.getSecondStageVec(scen);

		int nrow2 = data.dims.numSecondStageCons(scen);
		int nvar2 = data.dims.inner.numSecondStageVars(scen);
		for (int i = 0; i < nrow2; i++) {
			double row = 0.0;
			for (j = Wstart[i]; j < Wstart[i+1]; j++) {
				row += Welts[j]*x2[Widx[j]];
			}
			for (j = Tstart[i]; j < Tstart[i+1]; j++) {
				row += Telts[j]*x1[Tidx[j]];
			}
			row -= x2[nvar2+i];
			row = abs(row);
			if (row > largest) { 
				largest = row;
			}
		}
	
			
	}
	double l;
	MPI_Allreduce(&largest,&l,1,MPI_DOUBLE,MPI_MAX,data.ctx.comm());

	return l;
}


// load status from file, instead of taking starting basis from data
// file base is base of filename including iteration number, if applicable
void BALPSolverBase::loadStatus(const string &filebase) {

	unsigned filelen; char *filedata;
	if (data.ctx.mype() == 0) {
		string f0name = filebase + "0";
		ifstream f0(f0name.c_str());
		assert(f0.is_open());
		
		string dat((istreambuf_iterator<char>(f0)),istreambuf_iterator<char>()); // read into string
		f0.close();
		filelen = dat.length() + 1;
		filedata = new char[filelen];
		memcpy(filedata,dat.c_str(),filelen);
		MPI_Bcast(&filelen,1,MPI_UNSIGNED,0,data.ctx.comm());
	} else {
		MPI_Bcast(&filelen,1,MPI_UNSIGNED,0,data.ctx.comm());
		filedata = new char[filelen];
	}
	MPI_Bcast(filedata,filelen,MPI_CHAR,0,data.ctx.comm());
	
	string data(filedata);
	delete [] filedata;
	istringstream f1(data);
	loadStatus2(-1,f1);


	const vector<int> &localScen = basicIdx.localScenarios();
	for (unsigned k = 1; k < localScen.size(); k++) {
		int scen = localScen[k];
		stringstream ss;
		ss << filebase << scen+1;
		ifstream f(ss.str().c_str());
		assert(f.is_open());
		loadStatus2(scen,f);
	}
	status = LoadedFromFile;

}
