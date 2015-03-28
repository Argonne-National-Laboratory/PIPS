#include "BALPSolverPrimal.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include "CoinTime.hpp"

#include "PIPSLogging.hpp"

using namespace std;

// structs for MPI
struct doubledouble { double d1,d2; };
struct doubleint { double d; int i;};

BALPSolverPrimal::BALPSolverPrimal(const BAData &data) : BALPSolverBase(data),
	didperturb(false), devexPricing(true)
{

	uPerturb.allocate(data.dims,data.ctx, PrimalVector);
	lPerturb.allocate(data.dims,data.ctx, PrimalVector);	
	edgeWeights.allocate(data.dims,data.ctx, PrimalVector);
	reference.allocate(data.dims,data.ctx, PrimalVector);

	uPerturb.copyFrom(data.u);
	lPerturb.copyFrom(data.l);
	packed.allocate(data.dims,data.ctx);

	dualInfeas.allocate(data.dims,data.ctx, PrimalVector);

	//reinvertFrequency = 10;

}



void BALPSolverPrimal::initializeIndexedInfeasibilities() {
	assert(status != Uninitialized && status != LoadedFromFile);
	dualInfeas.clear();

	int ninfeas = 0;
	double suminf = 0.0;
	double inf0 = 0.0;

//	FILE *f = fopen("../balplu","w");
	const vector<int> &localScen = dualInfeas.localScenarios();
	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		CoinIndexedVector &dualInfeas1 = dualInfeas.getVec(scen);
		const denseVector &d1 = d.getVec(scen);
		const denseFlagVector<variableState> &states1 = states.getVec(scen);
		BAIndex baidx;
		baidx.scen = scen;
		int nvar = data.dims.numVars(scen);
		for (int i = 0; i < nvar; i++) {
			baidx.idx = i;
			if (states1[i] == Basic) continue;
			infeasType t = checkInfeas(baidx); // there is a better way to organize checkInfeas
			if (t != NotInfeasible) {
				double infeas = fabs(d1[i]);
				if (scen != -1) suminf += infeas;
				else inf0 += infeas;
				infeas *= infeas;
				assert(infeas);
				dualInfeas1.insert(i,infeas);
				ninfeas++;
			//	fprintf(f,"%s %e %e %e\n",data.names.getVec(scen)[idx].c_str(),x1[idx],l1[idx],u1[idx]);
			}
		}
	}

	int ninfeas1 = dualInfeas.getFirstStageVec().getNumElements();
	double allsuminf = data.ctx.reduce(suminf) + inf0;
	ninfeas = data.ctx.reduce(ninfeas - ninfeas1) + ninfeas1;
	if (data.ctx.mype() == 0) 
	  PIPS_ALG_LOG_SEV(info)<<boost::format("PIPS-S %d  Obj: %.9g Dual inf %.8g (%d) Elapsed %.1f") % nIter % objval % allsuminf % ninfeas % (MPI_Wtime()-startTime);

//	fclose(f);

}


BAIndex BALPSolverPrimal::price() const {

	const vector<int> &localScen = dualInfeas.localScenarios();

	//long double rmax = 0.0;
	double rmax = 0.0;
	BAIndex maxidx = {-1,-1};
	// don't need to worry about variables that
	// are actually feasible, because they can't be the max
	for (unsigned k = 0; k < localScen.size(); k++) {
		int scen = localScen[k];
		const CoinIndexedVector &dualInfeas1 = dualInfeas.getVec(scen);
		const int * COIN_RESTRICT idx = dualInfeas1.getIndices();
		const double * COIN_RESTRICT vec = dualInfeas1.denseVector();
		int ninfeas = dualInfeas1.getNumElements();
		const denseVector &weights1 = edgeWeights.getVec(scen);
		if (devexPricing) {
			for (int i = 0; i < ninfeas; i++) {
				int j = idx[i];
				double r = vec[j];
				if (r > rmax*weights1[j]) {
					rmax = r/weights1[j];
					maxidx.idx = j;
					maxidx.scen = scen;
				}
			}
		} else {
			for (int i = 0; i < ninfeas; i++) {
				int j = idx[i];
				double r = vec[j];
				if (r > rmax) {
					rmax = r;
					maxidx.idx = j;
					maxidx.scen = scen;
				}
			}
		}
	}
	doubleint my = { rmax, data.ctx.mype() }, best;
	MPI_Allreduce(&my,&best,1,MPI_DOUBLE_INT,MPI_MAXLOC,data.ctx.comm());
	MPI_Bcast(&maxidx,1,MPI_2INT,best.i,data.ctx.comm());

	return maxidx;
}

// Harris's two-pass ratio test
// enterType is template parameter to avoid unnecessary conditionals
template <infeasType enterType> BAIndex BALPSolverPrimal::ratioHarris(const sparseBAVector& pivotCol) {
	assert(enterType != NotInfeasible); // should be static_assert
	const vector<int> &localScen = dualInfeas.localScenarios();
	double thetaMax = 1e20;
	// pass 1
	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		const CoinIndexedVector& pivotCol1 = pivotCol.getVec(scen).v;
		const double * COIN_RESTRICT pivotElts = pivotCol1.denseVector();
		const int * COIN_RESTRICT pivotIdx = pivotCol1.getIndices();
		int nnz = pivotCol1.getNumElements();
		//const denseFlagVector<constraintType> &vartype1 = data.vartype.getVec(scen);
		const double * COIN_RESTRICT x1 = &x.getVec(scen)[0];
		const double * COIN_RESTRICT l1 = &lPerturb.getVec(scen)[0];
		const double * COIN_RESTRICT u1 = &uPerturb.getVec(scen)[0];
		const vector<int> &basicIdx1 = basicIdx.getVec(scen);
		packedVector &packed1 = packed.getPackedVec(j);
		packed1.resizeAndDestroy(nnz);
		double * COIN_RESTRICT packedElts = packed1.denseVector();
		int * COIN_RESTRICT packedIdx = packed1.getIndices();
		int nnzpacked = 0;

		for (int r = 0; r < nnz; r++) {
			int i = pivotIdx[r];
			int idx = basicIdx1[i];
			if (idx == -1) continue;
			double elt = pivotElts[i];
			if (enterType == Above) elt = -elt; // "static if"
			
			double ratio, realRatio;
			// compute both while we have data in cache
			if (elt > 0.0) {
				ratio = (x1[idx] - l1[idx] + primalTol)/elt;
				realRatio = (x1[idx] - l1[idx])/elt;
			} else {
				ratio = (x1[idx] - u1[idx] - primalTol)/elt;
				realRatio = (x1[idx] - u1[idx])/elt;
			}
		
			/*if (ratio < thetaMax) {
				thetaMax = ratio;
			}*/
			if (realRatio < thetaMax) {
				packedElts[nnzpacked] = realRatio;
				packedIdx[nnzpacked++] = i;
				if (ratio < thetaMax) {
					thetaMax = ratio;
				}
			}
		}
		packed1.setNumElements(nnzpacked);
	}

	// reduce (minimum) thetaMax here
	double tmax = thetaMax;
	MPI_Allreduce(&tmax,&thetaMax,1,MPI_DOUBLE,MPI_MIN,data.ctx.comm());

	BAIndex enter = { -1, -1 };
	if (thetaMax == 1e20) {
	  if (data.ctx.mype() == 0) PIPS_ALG_LOG_SEV(info) << "unbounded?";
		return enter; // unbounded
	}

	
	// pass 2
	double maxAlpha = 0.0;
	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		const CoinIndexedVector& pivotCol1 = pivotCol.getVec(scen).v;
		const double * COIN_RESTRICT pivotElts = pivotCol1.denseVector();
		const packedVector &packed1 = packed.getPackedVec(j);
		const double * COIN_RESTRICT packedElts = packed1.denseVector();
		const int * COIN_RESTRICT packedIdx = packed1.getIndices();
		int nnzpacked = packed1.getNumElements();


		for (int r = 0; r < nnzpacked; r++) {
			int i = packedIdx[r];
			double elt = pivotElts[i];
			elt = (elt > 0.0) ? elt : -elt;
			double ratio = packedElts[r];
			
			if (ratio <= thetaMax && elt > maxAlpha) {
				maxAlpha = elt;
				enter.scen = scen;
				enter.idx = i;
			}
		}
	}
	doubleint my = {maxAlpha, data.ctx.mype()}, best;
	MPI_Allreduce(&my,&best,1,MPI_DOUBLE_INT,MPI_MAXLOC,data.ctx.comm());
	MPI_Bcast(&enter,1,MPI_2INT,best.i,data.ctx.comm());
	

	return enter;

}

bool BALPSolverPrimal::quickCheckDualFeasible() const {
	
	const vector<int> &localScen = dualInfeas.localScenarios();
	int found = 0;
	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		const CoinIndexedVector &dualInfeas1 = dualInfeas.getVec(scen);
		const int *idx = dualInfeas1.getIndices();
		const double *vec = dualInfeas1.denseVector();
		int ninfeas = dualInfeas1.getNumElements();
		for (int i = 0; i < ninfeas; i++) {
			int k = idx[i];
			assert(vec[k]);
			if (vec[k] != FEASIBLE_TINY) {
				found = 1;
				break;
			}
		}
		if (found) break;
	}
	found = data.ctx.reduce(found);
	return !found;
}

// iteration of phase II primal simplex
void BALPSolverPrimal::iterate() {
	assert(status == PrimalFeasible || status == Optimal);
	if (quickCheckDualFeasible()) {
		status = Optimal;
		return;
	}

	double t = MPI_Wtime();
	BAIndex enterIdx = price();
	selectEnteringTime += MPI_Wtime() - t;
	infeasType enterType;
	double dualinf;
	if (data.ctx.assignedScenario(enterIdx.scen)) {
		enterType = checkInfeas(enterIdx);
		if (enterType == NotInfeasible) {
			assert(0);
		}
		dualinf = d[enterIdx];
		
			
	}
	if (enterIdx.scen != -1) MPI_Bcast(&dualinf,1,MPI_DOUBLE,data.ctx.owner(enterIdx.scen),data.ctx.comm());
	
	if (dualinf > 0.0) enterType = Above; // Nonbasic leaving upper bound
	else enterType = Below; // Nonbasic leaving lower bound

	// Pivot Column
	sparseBAVector &aq = ftranVec;
	aq.clear();
	data.getCol(aq,enterIdx);

	// FTRAN
	t = MPI_Wtime();
	la->ftran(aq);
	ftranTime += MPI_Wtime() -t;	

	// Ratio test

	t = MPI_Wtime();
	BAIndex leave;
	if (enterType == Above) leave = ratioHarris<Above>(aq);
	else leave = ratioHarris<Below>(aq);

	selectLeavingTime += MPI_Wtime() - t;
	if (leave.idx == -1) {
		status = ProvenUnbounded;
		return;
	}
	// global index of leaving variable
	BAIndex leaveGlobal;
	leaveGlobal.scen = leave.scen;
	leaveGlobal.idx = -1; // dummy value

	// BTRAN
	sparseBAVector &rho = btranVec;
	rho.clear();
	if (data.ctx.assignedScenario(leaveGlobal.scen)) {
		rho.getVec(leaveGlobal.scen).v.insert(leave.idx,1.0);
		leaveGlobal.idx = basicIdx.getVec(leaveGlobal.scen)[leave.idx];
	}
	t = MPI_Wtime();
	la->btranSimplex(rho,leave);
	btranTime += MPI_Wtime() - t;
	//btranDensity += rho.density();
	
	// Pivot row
	sparseBAVector &alpha = rowVec;
	alpha.clear();
	t = MPI_Wtime();
	data.multiplyT(rho,alpha);
	priceTime += MPI_Wtime() - t;
//	printf("PRICE from scen -1 has %d nonzeros\n",alpha.getVec(-1).v.getNumElements());

	

	//aq.checkClean();
	/*FILE *f = fopen("../balpvars","w");
	const vector<int> &localScen = aq.localScenarios();
	for (int j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		CoinIndexedVector &aq1 = aq.getVec(scen).v;
		aq1.scan();
		int *idx = aq1.getIndices();
		double *elts = aq1.denseVector();
		int n = aq1.getNumElements();
		for (int i = 0; i < n; i++) {
			fprintf(f,"%d: %e\n",idx[i],elts[idx[i]]);
		}
	}
	fclose(f);*/

	//ftranDensity += aq.density();


	// EXPAND
	double thetap; // primal step length: \hat x_B = x_B - thetap*(pivot column)
	double ftranPivot;

	infeasType leaveType;
	if (data.ctx.assignedScenario(leave.scen)) {
		ftranPivot = aq[leave];
		bool perturb = false;
		double &xval = x[leaveGlobal], &lval = lPerturb[leaveGlobal], &uval = uPerturb[leaveGlobal]; 
		if (enterType == Below) {
			if (ftranPivot > 0.0) {
				thetap = (xval - lval)/ftranPivot;
				leaveType = Below;
			} else {
				thetap = (xval - uval)/ftranPivot;
				leaveType = Above;
			}
			if (thetap < 0.0) { thetap = 1e-12; perturb = true; }
		} else {
			if (ftranPivot > 0.0) {
				thetap = (xval - uval)/ftranPivot;
				leaveType = Above;
			} else {
				thetap = (xval - lval)/ftranPivot;
				leaveType = Below;
			}
			if (thetap > 0.0) { thetap = -1e-12; perturb = true; }
		}
		if (doreport && data.ctx.mype() == data.ctx.owner(leave.scen)) {
		  std::stringstream ss; 
		  ss << boost::format("Selected %d,%d,%d (%s) to leave (x: %e, pivot: %e) ")
		    % leave.scen % leaveGlobal.idx % leave.idx % data.names[leaveGlobal].c_str() % xval % ftranPivot;
		  if (leaveType == Below) 
		    ss << boost::format("going to LB (%g)") % lval;
		  else 
		    ss << boost::format("going to UB (%g)") % uval;
		  PIPS_ALG_LOG_SEV(info) << ss.str();
		}
		xval -= thetap*ftranPivot;
		if (perturb) {
			didperturb = true;
			if (leaveType == Below) {
				lval = xval;
			} else {
				uval = xval;
			}
		}
	}
	if (leave.scen != -1) {
		doubledouble dd;
		if (data.ctx.mype() == data.ctx.owner(leave.scen)) {
			dd.d1 = ftranPivot;
			dd.d2 = thetap;
		}
		MPI_Bcast(&dd,sizeof(doubledouble),MPI_BYTE,data.ctx.owner(leave.scen),data.ctx.comm());
		ftranPivot = dd.d1;
		thetap = dd.d2;
		
	} 


	double thetad, btranPivot;
	if (data.ctx.assignedScenario(enterIdx.scen)) {
		thetad = d[enterIdx]/alpha[enterIdx];
		
	}
	if (enterIdx.scen != -1) {
		doubledouble dd;
		if (data.ctx.mype() == data.ctx.owner(enterIdx.scen)) {
			dd.d1 = alpha[enterIdx];
			dd.d2 = thetad;
		}
		MPI_Bcast(&dd,sizeof(doubledouble),MPI_BYTE,data.ctx.owner(enterIdx.scen),data.ctx.comm());
		btranPivot = dd.d1;
		thetad = dd.d2;
		
	} else {
		btranPivot = alpha[enterIdx];
	}
	if (data.ctx.owner(enterIdx.scen) == data.ctx.mype() && doreport) {
	  std::stringstream ss;
	  ss << boost::format("Selected %d,%d (%s) to enter (d: %g, thetap: %g) ") 
	    % enterIdx.scen % enterIdx.idx % data.names[enterIdx].c_str()  % dualinf  %  thetap;
	  if (enterType == Below) ss << "was at LB";
	  else ss << "was at UB";
	  PIPS_ALG_LOG_SEV(info) << ss.str();
	}


	

	if (fabs(ftranPivot-btranPivot) > 1e-9*(1+abs(ftranPivot))) {
	  if (data.ctx.mype() == 0) 
	    PIPS_ALG_LOG_SEV(warning)<<boost::format("Warning, pivot values from FTRAN and BTRAN differ significantly: %.7g %.7g") 
	      % ftranPivot % btranPivot;
		// trigger reinversion?
	}
	if (fabs(ftranPivot) < 1e-5 && la->nUpdates() > 0) {
		// drop iteration and reinvert
		//if (data.ctx.mype() == 0) printf("Bad pivot value (%e), dropping and reinverting\n", ftranPivot);
		initialize();
		return;
	}

	// basis change and update
	t = MPI_Wtime();
	updateDuals(alpha,leaveGlobal,enterIdx,thetad);
	updateIteratesTime += MPI_Wtime() - t;

	BAIndex enter = { enterIdx.scen, -1 };
	if (data.ctx.assignedScenario(enter.scen)) {
		enter.idx = static_cast<int>(basicIdx.getVec(enterIdx.scen).size());
	}

	if (data.ctx.assignedScenario(leave.scen)) {
		assert(basicIdx.getVec(leave.scen)[leave.idx] != -1); 
		basicIdx.getVec(leave.scen)[leave.idx] = -1; // indicate invalid index
	}

	t  = MPI_Wtime();
	if (devexPricing) {
		updatePrimals<true>(aq,enterIdx,thetap);
		updateDevexWeights(alpha,ftranPivot,enterIdx,leaveGlobal);
	} else {
		updatePrimals<false>(aq,enterIdx,thetap);
	}
	updateIteratesTime += MPI_Wtime() - t;

	if (data.ctx.assignedScenario(leave.scen)) {
		if (leaveType == Below) {
			x[leaveGlobal] = lPerturb[leaveGlobal];
			states[leaveGlobal] = AtLower;
		} else {
			x[leaveGlobal] = uPerturb[leaveGlobal];
			states[leaveGlobal] = AtUpper;
		}
	}
	if (data.ctx.assignedScenario(enter.scen)) {
		states[enterIdx] = Basic;
		basicIdx.getVec(enterIdx.scen).push_back(enterIdx.idx);
	}

	recordPivot(enterIdx.scen,leaveGlobal.scen);

	t = MPI_Wtime();
	int status = la->replaceColumnNoF(enter,leave,aq);
	updateColumnTime += MPI_Wtime() - t;
	/*if (status != 0 && data.ctx.mype() == 0) {
		printf("bad pivot, going to reinvert\n");
	}*/
	
	if (fabs(ftranPivot) < 1e-4) {
		// keep iteration but reinvert
		//if (data.ctx.mype() == 0) printf("Dangerous pivot value (%e), reinverting\n", ftranPivot);
		nIter += 1;
		initialize();
		nIter -= 1;
		return;
	}
	if (status != 0) {
		initialize();
	}





}


void BALPSolverPrimal::updateDuals(const sparseBAVector &alpha, BAIndex leaveIdx, BAIndex enterIdx, const double thetad) {

	const vector<int> localScen = x.localScenarios();
	for (unsigned k = 0; k < localScen.size(); k++) {
		int scen = localScen[k];
		const CoinIndexedVector& alpha1 = alpha.getVec(scen).v;
		const denseFlagVector<variableState> &states1 = states.getVec(scen);
		denseVector &d1 = d.getVec(scen);

		int alphaNelts = alpha1.getNumElements();
		const int *alphaIdx = alpha1.getIndices();
		const double *alphaElts = alpha1.denseVector();
		for (int j = 0; j < alphaNelts; j++) {
			int idx = alphaIdx[j];
			if (states1[idx] == Basic) continue;
			d1[idx] -= thetad*alphaElts[idx];
			updateInfeas(scen, idx);
		}
	}
	if (data.ctx.assignedScenario(enterIdx.scen)) {
		dualInfeas.getVec(enterIdx.scen).denseVector()[enterIdx.idx] = FEASIBLE_TINY;
		d[enterIdx] = 0.0;
	}
	if (data.ctx.assignedScenario(leaveIdx.scen)) {
		d[leaveIdx] = -thetad;
		updateInfeas(leaveIdx.scen, leaveIdx.idx);
	}

}


// enterIdx: global index of entering variable
// leave: old basic index of leaving variable
template<bool doDevex> void BALPSolverPrimal::updatePrimals(const sparseBAVector &aq, BAIndex enterIdx, double thetap) {
	const vector<int> &localScen = x.localScenarios();
	pivotDevex = 0.0;
	for (unsigned k = 0; k < localScen.size(); k++) {
		int scen = localScen[k];
		const CoinIndexedVector &aq1 = aq.getVec(scen).v;
		double * COIN_RESTRICT x1 = &x.getVec(scen)[0];
		double * COIN_RESTRICT l1 = &lPerturb.getVec(scen)[0];
		double * COIN_RESTRICT u1 = &uPerturb.getVec(scen)[0];
		const vector<int> &basicIdx1 = basicIdx.getVec(scen);
		const int * COIN_RESTRICT aqIdx = aq1.getIndices();
		const double * COIN_RESTRICT aqElts = aq1.denseVector();
		const denseFlagVector<bool> &reference1 = reference.getVec(scen);
		int aqNelts = aq1.getNumElements();

		for (int j = 0; j < aqNelts; j++) {
			int i = aqIdx[j];
			int idx = basicIdx1[i];
			if (idx == -1) { /*assert(aqElts[i] <= COIN_INDEXED_TINY_ELEMENT);*/ continue; }
			//if (scen == leave.scen && i == leave.idx) continue;
			x1[idx] -= thetap*aqElts[i];
			if (x1[idx] < l1[idx] - primalTol) {
				l1[idx] = x1[idx]; didperturb = true;
			} else if (x1[idx] > u1[idx] + primalTol) {
				u1[idx] = x1[idx]; didperturb = true;
			}
			if (doDevex) {
				if (reference1[idx]) pivotDevex += aqElts[i]*aqElts[i];
			}
				
		}
	}
	if (data.ctx.assignedScenario(enterIdx.scen)) {
		x[enterIdx] += thetap;
	}
	if (doDevex) { 
		if (data.ctx.mype() == data.ctx.owner(enterIdx.scen)) // only add once
			if (reference[enterIdx]) pivotDevex += 1.0;
	}
	
}


void BALPSolverPrimal::updateDevexWeights(const sparseBAVector &pivotRow, double pivotValue,
		BAIndex enterIdx, BAIndex leaveIdx) {
	
	// exact (squared) devex weight of entering column
	double devex = data.ctx.reduce(pivotDevex);
	int doreset = 0;
	if (data.ctx.assignedScenario(enterIdx.scen)) {
		double ratio = devex/edgeWeights[enterIdx];
		if ((nIter - devexIter) > 25 && (ratio > 3.0 || ratio < 0.333333)) {
			doreset = 1;
		}
	}
	if (enterIdx.scen != -1) doreset = data.ctx.reduce(doreset);
	if (doreset) {
		resetDevex(); return;
	}

	double invpivot = 1.0/pivotValue;
	const vector<int> localScen = x.localScenarios();
	for (unsigned k = 0; k < localScen.size(); k++) {
		int scen = localScen[k];
		const CoinIndexedVector& alpha1 = pivotRow.getVec(scen).v;
		const denseFlagVector<variableState> &states1 = states.getVec(scen);
		double * COIN_RESTRICT weights1 = &edgeWeights.getVec(scen)[0];
		int alphaNelts = alpha1.getNumElements();
		const int *alphaIdx = alpha1.getIndices();
		const double *alphaElts = alpha1.denseVector();

		for (int j = 0; j < alphaNelts; j++) {
			int idx = alphaIdx[j];
			if (states1[idx] == Basic) continue;
			double value = alphaElts[idx]*invpivot;
			weights1[idx] = max(0.99*weights1[idx],value*value*devex);
		}
	}
	if (data.ctx.assignedScenario(leaveIdx.scen)) {
		edgeWeights[leaveIdx] = max(1.0,invpivot*invpivot*devex);
	}
}

void BALPSolverPrimal::go() {
	//doreport = true;
	if (status == Uninitialized) {
		assert(0 && "Primal Phase I not implemented. Must specify starting basis.");
	}
	if (status == LoadedFromFile) {
	  if (data.ctx.mype() == 0) PIPS_ALG_LOG_SEV(info) << "first invert";
	  if (startTime == 0.0) startTime = MPI_Wtime();
	  initialize();
	}
	bool basischange = false ;
	if (startTime == 0.0) startTime = MPI_Wtime();
	// infeasible
	if (status == Initialized) {
		makeFeasible(basischange); 
	}
	if (status == DualFeasible) return;
	
	if (status == PrimalFeasible) {
		resetDevex(); 		
		for (; nIter < 100000000; nIter++) {
			//if (nIter % reportFrequency == 0) doreport = 1;
			if (phase1) {
			  if (doreport && data.ctx.mype() == 0) {
				  //printf("Iteration %d (Phase I)\n---------\nObjval:\t%f\n",nIter,objval);
			    PIPS_ALG_LOG_SEV(info) << "Iteration " << nIter <<  " (Phase I)\n---------";
			    PIPS_ALG_LOG_SEV(info) << "Objval:\t" << nIter;
			  }
			} else {
				if (doreport) {
					//largestPrimalError = calculateLargestPrimalError();
					if (data.ctx.mype() == 0) {
					  //printf("Iteration %d\n---------\n",nIter);
					  PIPS_ALG_LOG_SEV(info) << "Iteration " << nIter << "\n---------";
					}
				}
			}
			if ( (nIter - lastReinvert >= reinvertFrequency) && nIter > 0) {
				initialize(); // reinvert
				if (status == PrimalFeasible) { 
					lastGoodReinvert = nIter;
					if ((nIter - lastBadReinvert)/reinvertFrequency >= 10) {
						reinvertFrequency = min(reinvertFrequency + 5, reinvertFrequency_backup);
					}
				}
			}
			if (status == Initialized || status == DualFeasible) {
				lastBadReinvert = nIter;
				int nSteps = nIter - lastReinvert;
				if (reinvertFrequency_backup > 1) {
					if (reinvertFrequency == reinvertFrequency_backup) reinvertFrequency /= 2;
					else if (nSteps > 2) reinvertFrequency = nSteps - 2;
					else reinvertFrequency = 1;
				}
				forcePrimalFeasible();
			} 
			
			if (status != Optimal) iterate();
			
			if (status == ProvenUnbounded) {
				// reinvert to make sure
				if (data.ctx.mype() == 0) PIPS_ALG_LOG_SEV(info) << "Got unbounded, but going to retry";
				initialize();
				iterate();
				if (status == ProvenUnbounded) {
				  PIPS_ALG_LOG_SEV(info) <<"UNBOUNDED";
				  return;
				}
			}
			if (status == Optimal) {
				if (wasOptimal != nIter) {
					if (data.ctx.mype() == 0) 
					  PIPS_ALG_LOG_SEV(info) << "might be optimal, let's check";
					wasOptimal = nIter;
					nIter--; 
					initialize();
					continue;
				} else {
					break; // really optimal
				}
			}
		
			

			//if (doreport && data.ctx.mype() == 0) cout << endl;
			doreport = 0;
			
		}
		if (status == Optimal) {
		  double execTime = MPI_Wtime() - startTime;
		  if (!phase1) removePerturbation();
		  if (data.ctx.mype() == 0 && status == Optimal) {
		    PIPS_ALG_LOG_SEV(summary) << "Optimal!!! Objective value: " << objval << " Iterations: " << nIter;
		    //printf("Solution time: %f sec\n", execTime);
		    PIPS_ALG_LOG_SEV(info)<<"Pivot record:";
		    PIPS_ALG_LOG_SEV(info)<<"First replaces first: " << replaceFirst;
		    PIPS_ALG_LOG_SEV(info)<<"First replaces second: " << firstReplaceSecond;
		    PIPS_ALG_LOG_SEV(info)<<"Second replaces first: " << secondReplaceFirst;
		    PIPS_ALG_LOG_SEV(info)<<"Second replaces second (same scenario): " << replaceSecondSelf;
		    PIPS_ALG_LOG_SEV(info)<<"Second replaces second (other scenario): " << replaceSecondOther;
		    PIPS_ALG_LOG_SEV(info)<<"Time in FTRAN: " << boost::format("%f sec (%.1f%%)") 
		      % ftranTime % (100*ftranTime/execTime);
		    PIPS_ALG_LOG_SEV(info)<<boost::format("Time in INVERT: %f sec (%.1f%%)") 
		      % invertTime % (100*invertTime/execTime);
		    PIPS_ALG_LOG_SEV(info)<<boost::format("Time in FTRAN-DSE: %f sec (%.1f%%)") 
		      % ftranDSETime % (100*ftranDSETime/execTime);
		    PIPS_ALG_LOG_SEV(info)<<boost::format("Time in PRICE: %f sec (%.1f%%)") 
		      % priceTime % (100*priceTime/execTime);
		    PIPS_ALG_LOG_SEV(info)<<boost::format("Time selecting entering variable: %f sec (%.1f%%)") 
		      % selectEnteringTime % (100*selectEnteringTime/execTime);
		    PIPS_ALG_LOG_SEV(info)<<boost::format("Time selecting leaving variable: %f sec (%.1f%%)") 
		      % selectLeavingTime % (100*selectLeavingTime/execTime);
		    PIPS_ALG_LOG_SEV(info)<<boost::format("Time updating iterates: %f sec (%.1f%%)") 
		      % updateIteratesTime % (100*updateIteratesTime/execTime);
		    PIPS_ALG_LOG_SEV(info)<<boost::format("Time replacing column in basis: %f sec (%.1f%%)") 
		      % updateColumnTime % (100*updateColumnTime/execTime);
		    double recordedTime = ftranTime + btranTime + invertTime + 
		      ftranDSETime + priceTime + selectEnteringTime +
		      selectLeavingTime + updateIteratesTime + updateColumnTime;
		    
		    PIPS_ALG_LOG_SEV(info)<<boost::format("Time in other: %f sec (%.1f%%)") 
		      % (execTime - recordedTime) % (100*(execTime - recordedTime)/execTime);
		  }
		} else {
		  if (data.ctx.mype() == 0) 
		    PIPS_ALG_LOG_SEV(info) << "Didn't terminate after 100,000,000 iterations. Last objective: " << objval;
		}
	}
	


}

void BALPSolverPrimal::forcePrimalFeasible() {
	
	const vector<int> &localScen = x.localScenarios();
	double maxperturb = 0.0;

	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		const denseVector &x1 = x.getVec(scen);
		denseVector &l1 = lPerturb.getVec(scen);
		denseVector &u1 = uPerturb.getVec(scen);
		int nvar = data.dims.numVars(scen);
		for (int i = 0; i < nvar; i++) {
			if (x1[i] > u1[i] + primalTol) {
				maxperturb = max(maxperturb, x1[i]-u1[i]);
				u1[i] = x1[i];
			}
			if (x1[i] < l1[i] - primalTol) {
				maxperturb = max(maxperturb, l1[i]-x1[i]);
				l1[i] = x1[i];
			}
		}
	}
	double allmax;
	MPI_Allreduce(&maxperturb,&allmax,1,MPI_DOUBLE,MPI_MAX,data.ctx.comm());
	if (data.ctx.mype() == 0) 
	  PIPS_ALG_LOG_SEV(info)<<boost::format("Lost primal feasibility during reinversion. Added perturbations. Largest: %.7g") % allmax;
	if (status == Initialized) {
		status = PrimalFeasible;
	} else {
		assert(status == DualFeasible);
		status = Optimal;
	}

}

// reset reference framework for Devex (RECTIFY)
void BALPSolverPrimal::resetDevex() {
	const vector<int> &localScen = basicIdx.localScenarios();
	for (unsigned k = 0; k < localScen.size(); k++) {
		int scen = localScen[k];
		const denseFlagVector<variableState> &states1 = states.getVec(scen);
		denseVector &weights1 = edgeWeights.getVec(scen);
		denseFlagVector<bool> &reference1 = reference.getVec(scen);

		int nvar = data.dims.numVars(scen);
		for (int i = 0; i < nvar; i++) {
			weights1[i] = 1.0;
			reference1[i] = (states1[i] != Basic);
		}	
	}
	devexIter = nIter;

}

void BALPSolverPrimal::makeFeasible(bool &basischange) {
	if (status == PrimalFeasible) return;
	assert(0 && "Phase 1 not implemented");

}

void BALPSolverPrimal::removePerturbation() {
	didperturb = data.ctx.reduce(didperturb);
	if (!didperturb) {
	  if (data.ctx.mype() == 0) PIPS_ALG_LOG_SEV(info) << "No perturbations used";
		return;
	}
	if (data.ctx.mype() == 0) PIPS_ALG_LOG_SEV(info) << "Did perturbations";

	// may want to print some statistics on the perturbations
	lPerturb.copyFrom(data.l);
	uPerturb.copyFrom(data.u);
	initialize(false);
	if (data.ctx.mype() == 0) {
	  if (status == DualFeasible) {
	    PIPS_ALG_LOG_SEV(info) << "Became primal infeasible after removing perturbations, oops";
	  }
	}
}


#include <sstream>
#include <iomanip>



void BALPSolverPrimal::writeStatus(const string &filebase, bool appendIterateNumber,
	bool writeBasisOnly) {
	stringstream ss;
	ss << filebase;
	if (appendIterateNumber) ss << nIter << "iter";
	string s = ss.str();

	const vector<int> &localScen = basicIdx.localScenarios();
	
	for (unsigned k = 0; k < localScen.size(); k++) {
		int scen = localScen[k];
		if (scen == -1 && data.ctx.mype() != 0) continue;
		stringstream ss2;
		// first-stage is 0, second-stage uses 1-based indices
		ss2 << s << scen+1;
		ofstream f(ss2.str().c_str());
		//const denseVector &dse1 = dse.getVec(scen);
		const vector<int> basicIdx1 = basicIdx.getVec(scen);
		int nbasicThis = 0;
		for (unsigned i = 0; i < basicIdx1.size(); i++) {
			if (basicIdx1[i] >= 0) nbasicThis++;
		}
		f << nbasicThis;
		if (writeBasisOnly) f << " BasisOnly\n";
		else f << " FullDump " << nIter << "\n";
		f.precision(15);
		int nvar = data.dims.numVars(scen);
		const denseFlagVector<variableState> &states1 = states.getVec(scen);
		const denseVector &l1 = myLB().getVec(scen);
		const denseVector &u1 = myUB().getVec(scen);
		const denseVector &c1 = myObjective().getVec(scen);
		for (int i = 0; i < nvar; i++) {
			f << fixed << i;
			if (states1[i] == Basic) {
				f << " Basic";
			} else if (states1[i] == AtLower) {
				f << " AtLower";
			} else { 
				f << " AtUpper";
			}
			if (!writeBasisOnly) {
				f << " " << scientific << 1.0 << " " << l1[i] << " " << u1[i] << " " << c1[i];
			}
			f << "\n";
			
		}
		f.close();
	}
}

void BALPSolverPrimal::loadStatus2(int scen, istream& f) {
	int nvar = data.dims.numVars(scen);
	int nbasic1;
	//denseVector &dse1 = dse.getVec(scen);
	vector<int> &basicIdx1 = basicIdx.getVec(scen);
	basicIdx1.clear();
	denseFlagVector<variableState> &states1 = states.getVec(scen);
	denseVector &l1 = lPerturb.getVec(scen);
	denseVector &u1 = uPerturb.getVec(scen);

	string line;
	getline(f,line);
	enum { BasisOnly, FullDump, OldFormat } format;
	if (line.find("BasisOnly") != string::npos) {
		format = BasisOnly;
	} else if (line.find("FullDump") != string::npos) {
		format = FullDump;
	} else {
		format = OldFormat;
	}
	istringstream iss(line);
	iss >> nbasic1;
	if (format != OldFormat) {
		string s;
		iss >> s;
		if (format == FullDump) {
			int iter;
			iss >> iter;
			if (scen == -1) nIter = iter+1;
		}
	}

	int idx; double d1, d2; string status;
	for (int i = 0; i < nvar; i++) {
		assert(!f.eof());
		f >> idx;
		f >> status;
		if (format == OldFormat) f >> d1;
		else if (format == FullDump) {
			f >> d1;
			f >> l1[idx];
			f >> u1[idx];
			f >> d2;
		}
		if (status == "Basic") {
			states1[idx] = Basic;
			//dse1[basicIdx1.size()] = d;
			basicIdx1.push_back(idx);
		} else if (status == "AtLower") {
			states1[idx] = AtLower;
		} else if (status == "AtUpper") {
			states1[idx] = AtUpper;
		}
	}
	assert(static_cast<int>(basicIdx1.size()) == nbasic1);
}
