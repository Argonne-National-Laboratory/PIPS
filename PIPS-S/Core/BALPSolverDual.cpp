#include "BALPSolverDual.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include "CoinTime.hpp"
#include "PIPSSInterface.hpp"

using namespace std;

// structs for MPI
struct doubledouble { double d1,d2; };
struct doubleint { double d; int i;};

BALPSolverDual::BALPSolverDual(const BAData &data) : BALPSolverBase(data), DSEPricing(true),
	didperturb(false)
{


	// DSE weights are usually identified with rows, but
	// we need to associate them with (basic) primals
	dse.allocate(data.dims,data.ctx, BasicVector);
	cPerturb.allocate(data.dims,data.ctx, PrimalVector);
	cPerturb.copyFrom(data.c);

	primalInfeas.allocate(data.dims,data.ctx, BasicVector);
	

}


void BALPSolverDual::rebuildIndices() {

	if (!DSEPricing || status == Uninitialized) {
		setupIndices();
		return;
	}

	int nbasic = 0;
	const vector<int> &localScen = basicIdx.localScenarios();
	// map from global index to index in DSE vector, need this for copying
	int *map = 0, mapCapacity = 0;
	denseVector tmpDse(dse.getVec(-1).length());
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		vector<int> &basicIdx2 = basicIdx.getVec(scen);
		int nvar = data.dims.numVars(scen);
		
		// pass through old indices to form map
		// alternatively could maintain and update this map at each iteration
		if (nvar > mapCapacity) {
			if (map) delete [] map;
			mapCapacity = nvar;
			map = new int[mapCapacity];
			// for safety
			for (int k = 0; k < mapCapacity; k++) map[k] = -1;
		}
		for (unsigned j = 0; j < basicIdx2.size(); j++) {
			int k = basicIdx2[j];
			assert(k < mapCapacity);
			if (k >= 0) map[k] = j;
		}
		const denseFlagVector<variableState> &states2 = states.getVec(scen);
		denseVector &dse2 = dse.getVec(scen);
		basicIdx2.clear();
		if (dse2.length() > tmpDse.length()) {
			tmpDse.deallocate();
			tmpDse.allocate(dse2.length());
		}
		for (int j = 0; j < nvar; j++) {
			if (states2[j] == Basic) {
				assert(map[j] >= 0); // for safety
				tmpDse[basicIdx2.size()] = dse2[map[j]];
				map[j] = -1;
				//printf("%d %d %d %d\n", scen, j, nbasic,basicIdx2.size());
				basicIdx2.push_back(j);
				nbasic++;
			}
		}
		dse2.swap(tmpDse);
	}
	if (map) delete [] map;
	
	int nbasic1 = basicIdx.getFirstStageVec().size();
	nbasic = data.ctx.reduce(nbasic-nbasic1) + nbasic1;
	assert(nbasic == data.dims.totalCons());
	//if (data.ctx.mype() == 0) cout << "nbasic first stage: " << nbasic1 << endl;

}


void BALPSolverDual::flipBounds() {
	checkDualFeasible(infeasList);
	int didflip = 0;

	const vector<int> &localScen = infeasList.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		denseFlagVector<variableState> &states1 = states.getVec(scen);
		const denseFlagVector<constraintType> &vartype1 = data.vartype.getVec(scen);
		const vector<int> &infeas1 = infeasList.getVec(scen);	
		
		for (unsigned j = 0; j < infeas1.size(); j++) {
			int idx = infeas1[j];
			if (vartype1[idx] == Range || vartype1[idx] == Fixed) {
				if (states1[idx] == AtLower) states1[idx] = AtUpper;
				else if (states1[idx] == AtUpper) states1[idx] = AtLower;
				didflip++;
			}
		}
	}
	didflip = data.ctx.reduce(didflip);
	if (didflip) initialize(false);
}

void BALPSolverDual::initializeIndexedInfeasibilities() {
	assert(status != Uninitialized && status != LoadedFromFile);
	primalInfeas.clear();

	int ninfeas = 0;
	double suminf = 0.0;
	double inf0 = 0.0;

//	FILE *f = fopen("../balplu","w");
	const vector<int> &localScen = primalInfeas.localScenarios();
	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		CoinIndexedVector &primalInfeas1 = primalInfeas.getVec(scen);
		const denseVector &x1 = x.getVec(scen);
		const denseVector &l1 = data.l.getVec(scen);
		const denseVector &u1 = data.u.getVec(scen);
		const vector<int> &basicIdx1 = basicIdx.getVec(scen);
		BAIndex baidx;
		baidx.scen = scen;
		for (unsigned i = 0; i < basicIdx1.size(); i++) {
			int idx = basicIdx1[i];
			baidx.idx = idx;
			infeasType t = checkInfeas(baidx); // there is a better way to organize checkInfeas
			if (t != NotInfeasible) {
				double infeas;
				if (t == Below) {
					infeas = l1[idx]-x1[idx];
				} else {
					infeas = u1[idx]-x1[idx];
				}
				if (scen != -1) suminf += fabs(infeas);
				else inf0 += fabs(infeas);
				if (fabs(l1[idx]-u1[idx]) < 1e-5) infeas *= 100000.0;
				infeas *= infeas;
				assert(infeas);
				primalInfeas1.insert(i,infeas);
				ninfeas++;
			//	fprintf(f,"%s %e %e %e\n",data.names.getVec(scen)[idx].c_str(),x1[idx],l1[idx],u1[idx]);
			}
		}
	}

	int ninfeas1 = primalInfeas.getFirstStageVec().getNumElements();
	double allsuminf = data.ctx.reduce(suminf) + inf0;
	ninfeas = data.ctx.reduce(ninfeas - ninfeas1) + ninfeas1;
	if (data.ctx.mype() == 0) printf("Primal inf %.8g (%d) Elapsed %.1f\n",allsuminf,ninfeas,MPI_Wtime()-startTime);

//	fclose(f);

}


BAIndex BALPSolverDual::price() const {

	const vector<int> &localScen = primalInfeas.localScenarios();

	long double rmax = 0;
	BAIndex maxidx = {-1,-1};
	// don't need to worry about variables that
	// are actually feasible, because they can't be the max
	for (unsigned k = 0; k < localScen.size(); k++) {
		int scen = localScen[k];
		const CoinIndexedVector &primalInfeas1 = primalInfeas.getVec(scen);
		const int *idx = primalInfeas1.getIndices();
		const double *vec = primalInfeas1.denseVector();
		int ninfeas = primalInfeas1.getNumElements();
		const denseVector &dse1 = dse.getVec(scen);
		if (DSEPricing) {
			for (int i = 0; i < ninfeas; i++) {
				int j = idx[i];
				long double r = vec[j];
				if (r > rmax*dse1[j]) {
					rmax = r/dse1[j];
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

	// could also bcast delta here to save a communication

	return maxidx;
}

BAIndex BALPSolverDual::ratioHarris(const sparseBAVector& alpha2, double delta0) {
	
	const vector<int> &localScen = primalInfeas.localScenarios();
	BAContainer<vector<int> > &Q = infeasList;
	Q.clear();
	double thetaMax = 1e20;
	int qsize = 0;
	// pass 1
	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		const CoinIndexedVector& alpha = alpha2.getVec(scen).v;
		const double *alphaElts = alpha.denseVector();
		const int *alphaIdx = alpha.getIndices();
		int nnzAlpha = alpha.getNumElements();
		const denseFlagVector<variableState> &states1 = states.getVec(scen);
		const denseFlagVector<constraintType> &vartype1 = data.vartype.getVec(scen);
		vector<int> &Q1 = Q.getVec(scen);
		const denseVector& d1 = d.getVec(scen);
		
		for (int r = 0; r < nnzAlpha; r++) {
			int idx = alphaIdx[r];
			bool add = false;

			if (vartype1[idx] == Fixed || states1[idx] == Basic) continue;
			else if (states1[idx] == AtLower && alphaElts[idx] > pivotTol) add = true;
			else if (states1[idx] == AtUpper && alphaElts[idx] < -pivotTol) add = true;
			else if (vartype1[idx] == Free && (alphaElts[idx] > pivotTol || alphaElts[idx] < -pivotTol)) add = true;
			
			//printf("%s alpha: %e d: %e\n",data.names.getVec(scen)[idx].c_str(),alphaElts[idx],d1[idx]);
			if (add) {
				Q1.push_back(idx); qsize++;
				double ratio;
				
				if (alphaElts[idx] < 0) {
					ratio = (d1[idx] - dualTol)/alphaElts[idx];
				} else {
					ratio = (d1[idx] + dualTol)/alphaElts[idx];
				}
				//printf("00 ratio: %f (%d,%d)\n",ratio,scen,idx); 
				//printf("ratio: %e thetaMax: %e alpha2[i]: %e d[idx]: %e\n",ratio,thetaMax,alpha2[idx],d[idx]);
				if (ratio < thetaMax) {
					thetaMax = ratio;
				}
			}	
		}
	}

	// reduce (minimum) thetaMax here
	double tmax = thetaMax;
	MPI_Allreduce(&tmax,&thetaMax,1,MPI_DOUBLE,MPI_MIN,data.ctx.comm());

	BAIndex enter = { -1, -1 };
	if (thetaMax == 1e20) {
		//alpha2.print();
		if (data.ctx.mype() == 0) printf("unbounded?\n");
		return enter; // unbounded
	}

	
	// pass 2
	double maxAlpha = 0;
	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		const CoinIndexedVector& alpha = alpha2.getVec(scen).v;
		const double *alphaElts = alpha.denseVector();
		vector<int> &Q1 = Q.getVec(scen);
		const denseVector& d1 = d.getVec(scen);

		for (unsigned r = 0; r < Q1.size(); r++) {
			int idx = Q1[r];
			double ratio = d1[idx]/alphaElts[idx];
			if (ratio <= thetaMax) {
				double absAlpha = abs(alphaElts[idx]);
				//printf("ratio: %f absAlpha: %f\n", ratio, absAlpha);
				if (absAlpha > maxAlpha) {
					maxAlpha = absAlpha;
					enter.scen = scen;
					enter.idx = idx;
				}
			}/* else {
				printf("too big: %e\n",ratio);
			}*/
		}
	}
	doubleint my = {maxAlpha, data.ctx.mype()}, best;
	MPI_Allreduce(&my,&best,1,MPI_DOUBLE_INT,MPI_MAXLOC,data.ctx.comm());
	MPI_Bcast(&enter,1,MPI_2INT,best.i,data.ctx.comm());
	

	//printf("return normally: %d %d (qsize %d, thetaMax %e)\n", enter.scen, enter.idx,qsize,thetaMax);
	return enter;

}


/* Bound-flipping Ratio Test
   We need to be rework things here to allow for distributed vectors.
   Koberstein's BFRT requires too much communication,
   while ours has the same amount as the Harris ratio test.
   See notes in code for description of the procedure.

*/
struct ratioDeltaPair {
		ratioDeltaPair(const double &d, const double &d2) : ratio(d), delta(d2) {}
		bool operator<(const ratioDeltaPair &p) const { return (ratio < p.ratio); }
		double ratio, delta;
};
BAIndex BALPSolverDual::ratioBFRT(const sparseBAVector& alpha2, double delta0) {
	if (delta0 <= 1e-3) return ratioHarris(alpha2,delta0);
	
	const vector<int> &localScen = primalInfeas.localScenarios();
	BAContainer<vector<int> > Q(data.dims,data.ctx, PrimalVector);
	BAContainer<vector<int> > Qtilde(data.dims,data.ctx, PrimalVector); // *will* flip
	BAContainer<vector<int> > Qnew(data.dims,data.ctx, PrimalVector); // *won't* flip
	double thetaMax = 1e20;
	int qsize = 0; // number of elements in Q

	double minimumSlope = (delta0+1e-3)-(localScen.size()-1)*(delta0+1e-3)/data.dims.numScenarios();

	// Do Koberstein BFRT for local variables (including first-stage), allowing
	// slope to decrease to minimumSlope

	// Koberstein "Phase 1"
	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		const CoinIndexedVector& alpha = alpha2.getVec(scen).v;
		const double *alphaElts = alpha.denseVector();
		const int *alphaIdx = alpha.getIndices();
		int nnzAlpha = alpha.getNumElements();
		const denseFlagVector<variableState> &states1 = states.getVec(scen);
		const denseFlagVector<constraintType> &vartype1 = data.vartype.getVec(scen);
		vector<int> &Q1 = Q.getVec(scen);
		vector<int> &Qnew1 = Qnew.getVec(scen);
		const denseVector& d1 = d.getVec(scen);
		
		for (int r = 0; r < nnzAlpha; r++) {
			int idx = alphaIdx[r];
			bool add = false;

			if (vartype1[idx] == Fixed || states1[idx] == Basic) continue;
			else if (states1[idx] == AtLower && alphaElts[idx] > pivotTol) add = true;
			else if (states1[idx] == AtUpper && alphaElts[idx] < -pivotTol) add = true;
			else if (vartype1[idx] == Free && (alphaElts[idx] > pivotTol || alphaElts[idx] < -pivotTol)) add = true;
			/*if (states1[idx] == AtLower) printf("lower ");
			if (states1[idx] == AtUpper) printf("upper ");
			printf("%s (%d,%d) alpha: %e d: %e\n",data.names.getVec(scen)[idx].c_str(),scen,idx,alphaElts[idx],d1[idx]);*/
			if (add) {
				Q1.push_back(idx); qsize++;
				Qnew1.push_back(idx);
				double ratio;
				
				if (alphaElts[idx] < 0) {
					ratio = (d1[idx] - dualTol)/alphaElts[idx];
				} else {
					ratio = (d1[idx] + dualTol)/alphaElts[idx];
				}
				//printf("00 ratio: %f (%d,%d)\n",ratio,scen,idx); 
				//printf("ratio: %e thetaMax: %e alpha2[i]: %e d[idx]: %e\n",ratio,thetaMax,alpha2[idx],d[idx]);
				if (ratio < thetaMax) {
					thetaMax = ratio;
				}
			}	
		}
	}

	//printf("after phase 1, qsize = %d thetaMax %f\n",qsize, thetaMax);

	// Koberstein "Phase 2"
	// linear search through small groups to find those which will go infeasible/flip
	double thetaD = 3*(thetaMax+zeroTol);
	double deltaThis = 0.0;
	int qtildeSize = 0;
	
	vector<int> QnewTemp;
	while (delta0 - deltaThis >= minimumSlope && qsize > 0) {
		delta0 -= deltaThis;
		deltaThis = 0;
		thetaMax = 1e20;
		qsize = 0;
		qtildeSize = 0;
		// those that would go infeasible/flip if step length is thetaD
		for (unsigned j = 0; j < localScen.size(); j++) {
			int scen = localScen[j];
			const double *alphaElts = alpha2.getVec(scen).v.denseVector();
			const denseFlagVector<constraintType> &vartype1 = data.vartype.getVec(scen);
			vector<int> &Qtilde1 = Qtilde.getVec(scen);
			vector<int> &Qnew1 = Qnew.getVec(scen);
			const denseVector& d1 = d.getVec(scen);
			const denseVector& u1 = data.u.getVec(scen);
			const denseVector& l1 = data.l.getVec(scen);
			Qtilde1.clear();
			for (unsigned i = 0; i < Qnew1.size(); i++) {
				bool add = false;
				int idx = Qnew1[i];
				if (alphaElts[idx] > 0.0) {
					if (d1[idx] - thetaD*alphaElts[idx] < -dualTol) {
						add = true;
					}
				} else {
					if (d1[idx] - thetaD*alphaElts[idx] > dualTol) {
						add = true;
					}
				}
				if (add) {
					Qtilde1.push_back(idx); qtildeSize++;
					if (vartype1[idx] == Range) {
						deltaThis += (u1[idx]-l1[idx])*fabs(alphaElts[idx]);
					} else {
						// can't pass this variable
						deltaThis = 1e20;
					}
				} else {
					QnewTemp.push_back(idx); qsize++;
					double ratio;
					if (alphaElts[idx] < 0.0) {
						ratio = (d1[idx] - dualTol)/alphaElts[idx];
					} else {
						ratio = (d1[idx] + dualTol)/alphaElts[idx];
					}
					if (ratio < thetaMax) {
						thetaMax = ratio;
					}
				}
			}
			Qnew1.swap(QnewTemp);
			QnewTemp.clear();
		}
		// thetaMax is next breakpoint step amoung those which won't go infeasible
		// so thetaMax >= thetaD
		//printf("thetaD: %f thetaMax: %f will flip: %d won't flip: %d\n",thetaD,thetaMax,qtildeSize, qsize);
		//printf("delta0: %e, deltaThis: %e, minimumSlope: %e\n",delta0, deltaThis,minimumSlope); 
		assert(thetaMax >= thetaD || deltaThis >= 1e20);
		thetaD = 2*thetaMax;
	}
				

	// Koberstein "Phase 3" (modified)
	// careful pass through Qtilde from last step of phase 2
	// to find the actual breakpoint
	// don't perform harris here, unlike Koberstein
	// we build a list of (ratio,deltadiff) pairs and then sort them
	// there's probably a more efficient way
	// hopefully qtildeSize is small
	
	vector<ratioDeltaPair> pairs;
	pairs.reserve(qtildeSize);
	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		const double *alphaElts = alpha2.getVec(scen).v.denseVector();
		const denseFlagVector<constraintType> &vartype1 = data.vartype.getVec(scen);
		vector<int> &Qtilde1 = Qtilde.getVec(scen);
		const denseVector& u1 = data.u.getVec(scen);
		const denseVector& l1 = data.l.getVec(scen);
		const denseVector& d1 = d.getVec(scen);
		for (unsigned i = 0; i < Qtilde1.size(); i++) {
			int idx = Qtilde1[i];
			double ratio;
			if (alphaElts[idx] < 0.0) {
				ratio = (d1[idx] - dualTol)/alphaElts[idx];
			} else {
				ratio = (d1[idx] + dualTol)/alphaElts[idx];
			}
			double deltaDiff = (vartype1[idx] == Range) ? (u1[idx]-l1[idx])*fabs(alphaElts[idx]) : 1e20;
			pairs.push_back(ratioDeltaPair(ratio,deltaDiff));
		}
	}
	assert(static_cast<int>(pairs.size()) == qtildeSize);
	sort(pairs.begin(),pairs.end());

	thetaD = 1e20;
	deltaThis = 0.0;
	for (unsigned i = 0; i < pairs.size(); i++) {
		const ratioDeltaPair& p = pairs[i];
		deltaThis += p.delta;
		thetaD = p.ratio;
		if (delta0 - deltaThis <= minimumSlope) {
			break;
		}
	}


	// reduce (min) thetaD here
	double td = thetaD;
	MPI_Allreduce(&td,&thetaD,1,MPI_DOUBLE,MPI_MIN,data.ctx.comm());

	BAIndex enter = { -1, -1 };
	if (thetaD == 1e20) {
		printf("unbounded?\n");
		return enter;
	}

	// Now, find best pivot value among those with thetaD-2*dualTol/abs(alpha)<=ratio<=thetaD
	// the choice of 2 is a bit ad-hoc. it's a tradeoff between step size and numerical stability
	// unlike usual BFRT we don't have the option here of going to a previous breakpoint for better stability

	double maxAlpha = 0;
	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		const CoinIndexedVector& alpha = alpha2.getVec(scen).v;
		const double *alphaElts = alpha.denseVector();
		vector<int> &Q1 = Q.getVec(scen);
		const denseVector& d1 = d.getVec(scen);

		for (unsigned r = 0; r < Q1.size(); r++) {
			int idx = Q1[r];
			double ratio = d1[idx]/alphaElts[idx];
			double absAlpha = abs(alphaElts[idx]);
			if (thetaD-2.0*dualTol/absAlpha <= ratio && ratio <= thetaD) {
				//printf("ratio: %f absAlpha: %f\n", ratio, absAlpha);
				if (absAlpha > maxAlpha) {
					maxAlpha = absAlpha;
					enter.scen = scen;
					enter.idx = idx;
				}
			}
		}
	}
	doubleint my = { maxAlpha, data.ctx.mype() }, best;
	MPI_Allreduce(&my,&best,1,MPI_DOUBLE_INT,MPI_MAXLOC,data.ctx.comm());
	MPI_Bcast(&enter,1,MPI_2INT,best.i,data.ctx.comm());

	return enter;
}


bool BALPSolverDual::quickCheckPrimalFeasible() const {
	
	const vector<int> &localScen = primalInfeas.localScenarios();
	int found = 0;
	for (unsigned j = 0; j < localScen.size(); j++) {
		int scen = localScen[j];
		const CoinIndexedVector &primalInfeas1 = primalInfeas.getVec(scen);
		const int *idx = primalInfeas1.getIndices();
		const double *vec = primalInfeas1.denseVector();
		int ninfeas = primalInfeas1.getNumElements();
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

// iteration of phase II dual simplex
// calls bound-flipping ratio test
// then updates according to algorithm 5
void BALPSolverDual::iterate() {
	assert(status == DualFeasible || status == Optimal);
	if (quickCheckPrimalFeasible()) {
		status = Optimal;
		return;
	}

	double t = MPI_Wtime();
	BAIndex leave = price();
	selectLeavingTime += MPI_Wtime() - t;
	// global index of leaving variable
	BAIndex leaveGlobal;
	leaveGlobal.scen = leave.scen;
	leaveGlobal.idx = -1; // dummy value
	double delta;
	infeasType leaveType = Below; // dummy value
	if (data.ctx.assignedScenario(leave.scen)) {
		int leaveIdx = basicIdx.getVec(leave.scen)[leave.idx];
		assert(leaveIdx < data.dims.numVars(leave.scen));
		// only assigned proc needs global index and leaveType
		leaveGlobal.idx = leaveIdx;
		leaveType = checkInfeas(leaveGlobal);
		if (leaveType == NotInfeasible) {
			if (data.vartype[leaveGlobal] == Fixed) {
				leaveType = Below;
			} else {
				assert(0);
			}
		}
		if (leaveType == Below) {
			delta = x[leaveGlobal] - data.l[leaveGlobal];
			assert(delta < 0.0);
		} else {
			delta = x[leaveGlobal] - data.u[leaveGlobal];
			assert(delta > 0.0);
		}
		if (data.ctx.owner(leave.scen) == data.ctx.mype() && doreport) {
			printf("Selected %d,%d (%s) to leave (current: %e, delta: %e)\n", leaveGlobal.scen,leaveGlobal.idx,data.names[leaveGlobal].c_str(),x[leaveGlobal],delta);
			if (leaveType == Below) {
				printf("LB infeasible: l = %e\n",data.l[leaveGlobal]);
			} else {
				printf("UB infeasible: u = %e\n",data.u[leaveGlobal]);
			}
		}

	}
	if (leave.scen != -1) MPI_Bcast(&delta,1,MPI_DOUBLE,data.ctx.owner(leave.scen),data.ctx.comm());


	// BTRAN
	sparseBAVector &rho = btranVec;
	rho.clear();
	if (data.ctx.assignedScenario(leaveGlobal.scen)) {
		rho.getVec(leaveGlobal.scen).v.insert(leave.idx,1.0);
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

	// leaving from below
	if (delta < 0.0) {
		alpha.negate();
		// alpha is now alphaTilde, convenient for ratio test
	}

	// ratio test
	

	double absdelta = abs(delta);
	t = MPI_Wtime();
	BAIndex enterIdx = ratioHarris(alpha,absdelta);
	selectEnteringTime += MPI_Wtime() - t;
	//BAIndex enterIdx = ratioBFRT(alpha,absdelta);
	if (enterIdx.idx == -1) {
		status = ProvenInfeasible; // dual unbounded
		return;
	}
	


	sparseBAVector &aq = ftranVec;
	aq.clear();
	data.getCol(aq,enterIdx);

	// FTRAN
	t = MPI_Wtime();
	la->ftran(aq);
	ftranTime += MPI_Wtime() -t;
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

	// broadcast thetad?

	// EXPAND (section 6.2.2.3)
	double thetad, btranPivot;
	if (data.ctx.assignedScenario(enterIdx.scen)) {
		if (d[enterIdx]/alpha[enterIdx] < 0) {
			didperturb = true;
			if (delta < 0.0) {
				thetad = -1.0e-12;
			} else {
				thetad = 1.0e-12;
			}
			double diff = thetad*alpha[enterIdx] - d[enterIdx];
			d[enterIdx] = thetad*alpha[enterIdx];
			cPerturb[enterIdx] += diff;
			if (doreport && data.ctx.mype() == data.ctx.owner(enterIdx.scen)) 
				printf("did EXPAND\n");
		} else {
			if (delta < 0.0)
				thetad = -d[enterIdx]/alpha[enterIdx];
			else
				thetad = d[enterIdx]/alpha[enterIdx];
		}
		if (doreport && data.ctx.mype() == data.ctx.owner(enterIdx.scen)) {
			printf("Selected %d,%d (%s) to enter (d: %e, alpha: %e)\n",enterIdx.scen,enterIdx.idx,data.names[enterIdx].c_str(),d[enterIdx],alpha[enterIdx]);
		}
	}
	if (delta < 0.0) alpha.negate(); // get alpha back

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
	
	// basis change and update
	t = MPI_Wtime();
	updateDuals(alpha,leaveGlobal,enterIdx,thetad);
	updateIteratesTime += MPI_Wtime() - t;

	double thetap, ftranPivot;
	if (data.ctx.assignedScenario(leave.scen)) {
		double deltaNew = (leaveType == Below) ? (x[leaveGlobal] - data.l[leaveGlobal]) : (x[leaveGlobal] - data.u[leaveGlobal]);
		thetap = deltaNew/aq[leave];
	}
		
	if (leave.scen != -1) {
		doubledouble dd;
		if (data.ctx.mype() == data.ctx.owner(leave.scen)) {
			dd.d1 = aq[leave];
			dd.d2 = thetap;
		}
		MPI_Bcast(&dd,sizeof(doubledouble),MPI_BYTE,data.ctx.owner(leave.scen),data.ctx.comm());
		ftranPivot = dd.d1;
		thetap = dd.d2;
	} else {
		ftranPivot = aq[leave];
	}

	

	if (doreport && data.ctx.mype() == data.ctx.owner(leave.scen)) printf("pivot from FTRAN: %e (%d,%d)\n",aq[leave],leave.scen,leave.idx);
		if (fabs(ftranPivot) < 1e-5 && la->nUpdates() > 0) {
		// drop iteration and reinvert
		//if (data.ctx.mype() == 0) printf("Bad pivot value (%e), dropping and reinverting\n", ftranPivot);
		initialize();
		return;
	}

	BAIndex enter = { enterIdx.scen, -1 };
	if (data.ctx.assignedScenario(enter.scen)) {
		enter.idx = static_cast<int>(basicIdx.getVec(enterIdx.scen).size());
	}

	t  = MPI_Wtime();
	updatePrimals(aq,enterIdx,enter,leave,thetap);
	updateIteratesTime += MPI_Wtime() - t;

	if (data.ctx.assignedScenario(leave.scen)) {
		assert(basicIdx.getVec(leave.scen)[leave.idx] != -1); 
		basicIdx.getVec(leave.scen)[leave.idx] = -1; // indicate invalid index
		if (delta < 0.0) {
			x[leaveGlobal] = data.l[leaveGlobal];
			states[leaveGlobal] = AtLower;
		} else {
			x[leaveGlobal] = data.u[leaveGlobal];
			states[leaveGlobal] = AtUpper;
		}
	}
	if (data.ctx.assignedScenario(enter.scen)) {
		states[enterIdx] = Basic;
		basicIdx.getVec(enterIdx.scen).push_back(enterIdx.idx);
	}

	// can experiment with btranPivot vs. ftranPivot here
	if (DSEPricing) updateDSE(rho, aq, enter, btranPivot);



	recordPivot(enterIdx.scen,leaveGlobal.scen);

	t = MPI_Wtime();
	int status = la->replaceColumnNoF(enter,leave,aq);
	updateColumnTime += MPI_Wtime() - t;
	/*if (status != 0 && data.ctx.mype() == 0) {
		printf("bad pivot, going to reinvert\n");
	}*/
	// perform additional stability test (6.29)
	/*if (abs(alpha[enterIdx] - aq[leave]) > 1e-9*(1+abs(aq[leave]))) {
		printf("stability test on pivots failed, going to reinvert\n");
		status = 1;
	}*/
	if (fabs(ftranPivot-btranPivot) > 1e-8*(1+fabs(ftranPivot)) || fabs(ftranPivot) < 1e-4) {
		//if (data.ctx.mype() == 0) printf("Warning, pivot values from FTRAN and BTRAN differ significantly: %.7g %.7g\n",ftranPivot,btranPivot);
		// keep iteration but reinvert
		nIter += 1;
		initialize();
		nIter -= 1;
		return;
	}
	if (status != 0) {
		initialize();
	}





}


// updates duals and performs bounds flipping (therefore updating primals also)
// returns change in objective from bounds flips in objdelta
void BALPSolverDual::updateDuals(const sparseBAVector &alpha, BAIndex leaveIdx, BAIndex enterIdx, const double thetad) {

	if (d.hasScenario(leaveIdx.scen)) d[leaveIdx] = -thetad;
	if (d.hasScenario(enterIdx.scen)) d[enterIdx] = 0;
	sparseBAVector &atilde = ftranVec2; // rhs for FTRAN
	atilde.clear();
	int nFlipped = 0;

	const vector<int> localScen = x.localScenarios();
	for (unsigned k = 0; k < localScen.size(); k++) {
		int scen = localScen[k];
		const CoinIndexedVector& alpha1 = alpha.getVec(scen).v;
		denseFlagVector<variableState> &states1 = states.getVec(scen);
		const denseFlagVector<constraintType> &vartype1 = data.vartype.getVec(scen);
		//const denseVector &l1 = data.l.getVec(scen);
		//const denseVector &u1 = data.u.getVec(scen);
		//denseVector &x1 = x.getVec(scen);
		denseVector &d1 = d.getVec(scen);
		denseVector &c1 = cPerturb.getVec(scen);

		int alphaNelts = alpha1.getNumElements();
		const int *alphaIdx = alpha1.getIndices();
		const double *alphaElts = alpha1.denseVector();
		for (int j = 0; j < alphaNelts; j++) {
			int idx = alphaIdx[j];
			//BAIndex baidx = { scen, idx };
			if (states1[idx] == Basic || (scen == enterIdx.scen && idx == enterIdx.idx)) continue;
			double dnew = d1[idx] - thetad*alphaElts[idx];
			if (vartype1[idx] != Fixed) {
				/*
				// we check now if bounds need to be flipped
				if (vartype1[idx] == Range) {
					// use 0 or dualTol here?
					// these could cause negative objective changes if flipped by mistake
					d1[idx] = dnew;
					if (states1[idx] == AtLower && dnew < -dualTol) {
						states1[idx] = AtUpper;
						x1[idx] = u1[idx];
						data.addColToVec(atilde, baidx, u1[idx]-l1[idx]);
						nFlipped++;
					} else if (states1[idx] == AtUpper && dnew > dualTol) {
						states1[idx] = AtLower;
						x1[idx] = l1[idx];
						data.addColToVec(atilde, baidx, l1[idx]-u1[idx]);
						nFlipped++;
					}
				} else {*/
					// now deal with infeasibilities with cost shifting

					if (states1[idx] == AtLower || vartype1[idx] == Free) {
						if (dnew >= -dualTol) {
							d1[idx] = dnew;
						} else {
							double delta = -dnew - dualTol;
							c1[idx] += delta;
							d1[idx] = -dualTol;
							if (doreport) printf("did EXPAND lower\n");
						}
					}
					if (states1[idx] == AtUpper || vartype1[idx] == Free) {
						if (dnew <= dualTol) {
							d1[idx] = dnew;
						} else {
							double delta = -dnew + dualTol;
							c1[idx] += delta;
							d1[idx] = dualTol;
							if (doreport) printf("did EXPAND upper\n");
						}
					}
				//}
			} else {
				d1[idx] = dnew;
			}
					
		}
	}
	//nFlipped = data.ctx.reduce(nFlipped);
	if (nFlipped) {
		// need to do an ftran
		if (doreport && data.ctx.mype() == 0) printf("did %d flips\n", nFlipped);
		la->ftran(atilde);
		for (unsigned k = 0; k < localScen.size(); k++) {
			int scen = localScen[k];
			const CoinIndexedVector& atilde1 = atilde.getVec(scen).v;
			denseVector &x1 = x.getVec(scen);
			const vector<int> &basicIdx1 = basicIdx.getVec(scen);

			int atildeNelts = atilde1.getNumElements();
			const int *atildeIdx = atilde1.getIndices();
			const double *atildeElts = atilde1.denseVector();
			for (int j = 0; j < atildeNelts; j++) {
				int i = atildeIdx[j];
				int idx = basicIdx1[i];
				if (idx == -1) { assert(atildeElts[i] == COIN_INDEXED_REALLY_TINY_ELEMENT); continue; }
				x1[idx] -= atildeElts[i];
				updateInfeas(scen,i,idx);
			}


		}

	} else {
		if (doreport && data.ctx.mype() == 0) printf("no flips\n");
	}


}

// update DSE weights ("Dual algorithm I")
// see sections 3.3 and 8.2.2.1 in koberstein thesis
// enter is new basic index of entering column
void BALPSolverDual::updateDSE(sparseBAVector &rho, const sparseBAVector &aq, BAIndex enter, double pivot) {
	double dseleave = rho.dotWithSelf()/pivot/pivot;
	sparseBAVector &tau = rho; // overwrite
	double t = MPI_Wtime();
	la->ftran(tau);
	ftranDSETime += MPI_Wtime() - t;
	t = MPI_Wtime();
	//DSEDensity += tau.density(); 
	double kappa = -2.0/pivot; // premultiply tau by kappa?
	const vector<int> &localScen = x.localScenarios();
	for (unsigned k = 0; k < localScen.size(); k++) {
		int scen = localScen[k];
		const CoinIndexedVector &aq1 = aq.getVec(scen).v;
		denseVector &dse1 = dse.getVec(scen);
		const int *aqIdx = aq1.getIndices();
		const double *aqElts = aq1.denseVector();
		int aqNelts = aq1.getNumElements();
		const double *tauElts = tau.getVec(scen).v.denseVector();
		for (int j = 0; j < aqNelts; j++) {
			int i = aqIdx[j];
			dse1[i] += aqElts[i]*(aqElts[i]*dseleave + kappa*tauElts[i]);
			dse1[i] = max(dse1[i],1.0e-4);
		}
	}
	updateIteratesTime += MPI_Wtime() - t;
	if (data.ctx.assignedScenario(enter.scen)) dse[enter] = dseleave;
	// old weight is still in dse[leave], but it shouldn't be accessed
}

void BALPSolverDual::go() {
	//doreport = true;
	if (status == Uninitialized) {
		startTime = MPI_Wtime();
		setSlackBasis();
	}
	if (status == LoadedFromFile) {
		if (startTime == 0.0) startTime = MPI_Wtime();
		initialize();
	}
	if (startTime == 0.0) startTime = MPI_Wtime();
	// infeasible
	if (status == Initialized) {
		makeFeasible();
	}
	if (status == PrimalFeasible) return;
	
	if (status == DualFeasible) {
		//if (DSEPricing) initializeDSE(!basischange && slackbasis);
		if (DSEPricing) initializeDSE(true);
		for (; nIter < 100000000; nIter++) {
			//if (nIter % reportFrequency == 0) doreport = 1;
			if (phase1) {
				if (doreport && data.ctx.mype() == 0) 
					printf("Iteration %d (Phase I)\n---------\nObjval:\t%f\n",nIter,objval);
			} else {
				if (doreport) {
					//largestPrimalError = calculateLargestPrimalError();
					if (data.ctx.mype() == 0) {
						printf("Iteration %d\n---------\n",nIter);
					}
				}
			}
			if ( (nIter - lastReinvert >= reinvertFrequency) && nIter > 0) {
					initialize(); // reinvert
					if (status == DualFeasible) {
						lastGoodReinvert = nIter;
						if ((nIter - lastBadReinvert)/reinvertFrequency >= 10) {
							reinvertFrequency = min(reinvertFrequency + 5, reinvertFrequency_backup);
						}
					}
			}
			if (status == Initialized || status == PrimalFeasible) {
				lastBadReinvert = nIter;
				int nSteps = nIter - lastGoodReinvert;
				if (reinvertFrequency_backup > 1) {
					if (reinvertFrequency == reinvertFrequency_backup) reinvertFrequency /= 2;
					else if (nSteps > 2) reinvertFrequency = nSteps - 2;
					else reinvertFrequency = 1;
				}
				forceDualFeasible();
			}
			if (status != Optimal) iterate();
			//if ((nIter % 1000) == 0) writeBasis("../basis/basis");
			if (nIter && dumpEvery && (nIter%dumpEvery) == 0) { 
				double t = MPI_Wtime();
				writeStatus(outputName,true,false);
				startTime += MPI_Wtime() - t;
			}
			

			if (status == ProvenInfeasible) {
				// reinvert to make sure
				if (data.ctx.mype() == 0) printf("Got unbounded, but going to retry\n");
				initialize();
				iterate();
				if (status == ProvenInfeasible) {
					if (data.ctx.mype() == 0) printf("UNBOUNDED\n");
					return;
				}
			}
			if (status == Optimal) {
				if (wasOptimal != nIter) {
					if (data.ctx.mype() == 0) printf("might be optimal, let's check\n");
					wasOptimal = nIter;
					nIter--;
					initialize();
					continue;
				} else {
					break; // really optimal
				}
			}

			if (doreport && data.ctx.mype() == 0) cout << endl;
			doreport = 0;
			
		}
		if (status == Optimal) {
			double execTime = MPI_Wtime() - startTime;
			if (!phase1) removePerturbation();
			//printf("Average FTRAN-DSE Result Density: %f\n",DSEDensity/nIter);
			//printf("Average FTRAN Result Density: %f\n",ftranDensity/nIter);
			//printf("Average BTRAN Result Density: %f\n",btranDensity/nIter);
			if (data.ctx.mype() == 0 && status == Optimal) {
				printf("Optimal!!! Objective value: %f Iterations: %d\n", objval,nIter);
				//printf("Solution time: %f sec\n", execTime);
				printf("Pivot record:\n");
				printf("First replaces first: %d\n",replaceFirst);
				printf("First replaces second: %d\n",firstReplaceSecond);
				printf("Second replaces first: %d\n",secondReplaceFirst);
				printf("Second replaces second (same scenario): %d\n",replaceSecondSelf);
				printf("Second replaces second (other scenario): %d\n",replaceSecondOther);

				printf("Time in FTRAN: %f sec (%.1f%%)\n",ftranTime,100*ftranTime/execTime);
				printf("Time in BTRAN: %f sec (%.1f%%)\n",btranTime,100*btranTime/execTime);
				printf("Time in INVERT: %f sec (%.1f%%)\n",invertTime,100*invertTime/execTime);
				printf("Time in FTRAN-DSE: %f sec (%.1f%%)\n",ftranDSETime,100*ftranDSETime/execTime);
				printf("Time in PRICE: %f sec (%.1f%%)\n",priceTime,100*priceTime/execTime);
				printf("Time selecting entering variable: %f sec (%.1f%%)\n",selectEnteringTime,
					100*selectEnteringTime/execTime);
				printf("Time selecting leaving variable: %f sec (%.1f%%)\n",selectLeavingTime,
					100*selectLeavingTime/execTime);
				printf("Time updating iterates: %f sec (%.1f%%)\n",updateIteratesTime,
					100*updateIteratesTime/execTime);
				printf("Time replacing column in basis: %f sec (%.1f%%)\n",updateColumnTime,
					100*updateColumnTime/execTime);
				double recordedTime = ftranTime + btranTime + invertTime +
					ftranDSETime + priceTime + selectEnteringTime +
					selectLeavingTime + updateIteratesTime + updateColumnTime;
				printf("Time in other: %f sec (%.1f%%)\n",execTime - recordedTime,
					100*(execTime - recordedTime)/execTime);
				
			}

		} else {
			if (data.ctx.mype() == 0) 
				printf("Didn't terminate after 100,000,000 iterations. Last objective: %f\n", objval);
		}
	}
	


}

void BALPSolverDual::forceDualFeasible() {
	assert(status != Uninitialized && status != LoadedFromFile);

	const vector<int> &localScen = x.localScenarios();
	double maxperturb = 0.0;
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		const denseFlagVector<variableState> &states1 = states.getVec(scen);
		const denseFlagVector<constraintType> &vartype1 = data.vartype.getVec(scen);
		denseVector &d1 = d.getVec(scen);
		denseVector &c1 = cPerturb.getVec(scen);
		int nvar = data.dims.numVars(scen);
		for (int idx = 0; idx < nvar; idx++) {
			if (states1[idx] == Basic) continue;

			if (states1[idx] == AtLower || vartype1[idx] == Free) {
				if (d1[idx] < -dualTol) {
					double delta = -d1[idx] - dualTol;
					c1[idx] += delta;
					d1[idx] = -dualTol;
					maxperturb = max(maxperturb,delta);
				}
			} if (states1[idx] == AtUpper || vartype1[idx] == Free) {
				if (d1[idx] > dualTol) {
					double delta = -d1[idx] + dualTol;
					c1[idx] += delta;
					d1[idx] = dualTol;
					maxperturb = max(maxperturb,-delta);
				}
			}
		}
	}

	double allmax;
	MPI_Allreduce(&maxperturb,&allmax,1,MPI_DOUBLE,MPI_MAX,data.ctx.comm());
	if (data.ctx.mype() == 0) printf("Lost dual feasibility during reinversion. Added perterbations. Largest: %.7g\n",allmax);
	if (status == Initialized) {
		status = DualFeasible;
	} else {
		assert(status == PrimalFeasible);
		status = Optimal;
	}
}

void BALPSolverDual::makeFeasible() {
	if (status == DualFeasible) return;
	assert(status == Initialized);

	initialize(false); // make sure dual values are accurate
	flipBounds();
	if (checkDualFeasible(infeasList)) {
		if (data.ctx.mype() == 0) printf("was infeasible but got feasible after bounds flipping\n");
		status = DualFeasible; 
		return;
	}
	if (phase1) {
		// phase 1 infeasible, shouldn't get here
		assert(0);
	}
	/* From one of Julian Hall's students:

	Assuming one is using the computation form of LP model:
	min   c^T x
		 A x  = 0
	      L <= x <= U,
	the subproblem approach mentioned by Koberstein could be implemented  
	by changing all the bounds of the model:

	lower variable [lower, inf  ) -> [ 0, 1];
	upper variable ( -inf, upper] -> [-1, 0];
	boxed variable [lower, upper] -> [ 0, 0];       // Tend to be non-basic
	free  variable ( -inf, inf  ) -> [-1000, 1000]; // Tend to be basic

	This subproblem model is for so called "scaled dual infeasibility",  
	and the optimal value should be 0.
	*/

	// if memory becomes an issue, delete local "la" and reallocate after finished
	PIPSSInterface lp2(data,PIPSSInterface::useDual);
	BAData &dCopy = lp2.d;

	const vector<int> &localScen = x.localScenarios();
	for (unsigned k = 0; k < localScen.size(); k++) {
		int scen = localScen[k];
		denseVector &l = dCopy.l.getVec(scen);
		denseVector &u = dCopy.u.getVec(scen);
		denseFlagVector<constraintType> vartype = dCopy.vartype.getVec(scen);
		int nvar = data.dims.numVars(scen);
		for (int i = 0; i < nvar; i++) {
			constraintType t = vartype[i];
			vartype[i] = Range;
			if (t == Range || t == Fixed) {
				l[i] = 0.;
				u[i] = 0.;
				vartype[i] = Fixed;
			} else if (t == LB) {
				l[i] = 0.;
				u[i] = 1.;
			} else if (t == UB) {
				l[i] = -1.;
				u[i] = 0.;
			} else if (t == Free) {
				l[i] = -1000.;
				u[i] = 1000.;
			}
		}
	}
	
	lp2.setStates(states);
	if (data.ctx.mype() == 0) printf("Starting phase 1\n");
	lp2.go();
	assert(lp2.getStatus() == Optimal);
	if (fabs(lp2.getObjective()) > 1e-5) {
		if (data.ctx.mype() == 0) printf("Problem appears dual infeasible\n");
		status = ProvenUnbounded; return;
	}
	// degeneracy in the solution could lead to invalid states, should check for this.
	setStates(lp2.getStatesRef());
	initialize();

}

// enterIdx: global index of entering variable
// enter: new basic index of entering variable
// leave: old basic index of leaving variable
void BALPSolverDual::updatePrimals(const sparseBAVector &aq, BAIndex enterIdx, BAIndex enter, BAIndex leave, double thetap) {
	const vector<int> &localScen = x.localScenarios();
	for (unsigned k = 0; k < localScen.size(); k++) {
		int scen = localScen[k];
		const CoinIndexedVector &aq1 = aq.getVec(scen).v;
		denseVector &x1 = x.getVec(scen);
		const vector<int> &basicIdx1 = basicIdx.getVec(scen);
		const int *aqIdx = aq1.getIndices();
		const double *aqElts = aq1.denseVector();
		int aqNelts = aq1.getNumElements();

		for (int j = 0; j < aqNelts; j++) {
			int i = aqIdx[j];
			int idx = basicIdx1[i];
			if (idx == -1) { assert(aqElts[i] <= COIN_INDEXED_TINY_ELEMENT); continue; }
			x1[idx] -= thetap*aqElts[i];
			// this can be done more efficiently
			updateInfeas(scen,i,idx);
		}
	}
	if (x.hasScenario(enterIdx.scen)) {
		x[enterIdx] += thetap;
		updateInfeas(enter.scen,enter.idx,enterIdx.idx);
	}
	if (x.hasScenario(leave.scen)) {
		primalInfeas.getVec(leave.scen).denseVector()[leave.idx] = FEASIBLE_TINY;
	}
}

void BALPSolverDual::initializeDSE(bool slackbasis) {
	assert(status != Uninitialized);


	// assumes that all indices in basicIdx are valid,
	// so, can only call this after a fresh factorization

	if (slackbasis) {
		const vector<int> &localScen = x.localScenarios();
		for (unsigned k = 0; k < localScen.size(); k++) {
			int scen = localScen[k];
			denseVector &dse1 = dse.getVec(scen);
			int nbasic = basicIdx.getVec(scen).size();
			for (int i = 0; i < nbasic; i++) dse1[i] = 1;
		}
	} else {
		sparseBAVector &rho = btranVec;
		for (int scen = -1; scen < data.dims.numScenarios(); scen++) {
			CoinIndexedVector &rho1 = rho.getVec(scen).v;
			int nbasic = basicIdx.getVec(scen).size();
			for (int i = 0; i < nbasic; i++) { 
				rho.clear();
				rho1.insert(i,1.0);
				la->btran(rho);
				if (data.ctx.assignedScenario(scen)) dse.getVec(scen)[i] = rho.dotWithSelf();
			}
		}
	}
}

void BALPSolverDual::removePerturbation() {
	didperturb = data.ctx.reduce(didperturb);
	if (!didperturb) {
		if (data.ctx.mype() == 0) printf("No perturbations used\n");
		return;
	}
	if (data.ctx.mype() == 0) printf("Did perturbations\n");

	// may want to print some statistics on the perturbations
	//perturb.print();
	cPerturb.copyFrom(data.c);
	initialize(false);
	if (data.ctx.mype() == 0) {
		printf("Objective after removing perturbations: %f\n", objval);
		if (status == PrimalFeasible) {
			printf("Became dual infeasible after removing perturbations, oops\n");
		}
	}
}


#include <sstream>
#include <iomanip>



void BALPSolverDual::writeStatus(const string &filebase, bool appendIterateNumber,
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
		//ofstream f(ss2.str().c_str(), writeBasisOnly ? ios_base::out : ios_base::binary);
		ofstream f(ss2.str().c_str());
		const denseVector &dse1 = dse.getVec(scen);
		const vector<int> basicIdx1 = basicIdx.getVec(scen);
		int nbasicThis = 0;
		int nvar = data.dims.numVars(scen);
		vector<double> tmpDse(nvar,1.0);
		for (unsigned i = 0; i < basicIdx1.size(); i++) {
			if (basicIdx1[i] >= 0) {
				nbasicThis++;
				tmpDse[basicIdx1[i]] = dse1[i];
			}
		}
		// don't use endl here. need to buffer for parallel IO
		f << nbasicThis;
		if (writeBasisOnly) f << " BasisOnly\n";
		else f << " FullDump " << nIter << "\n";
		f.precision(std::numeric_limits<double>::digits10+3);
		const denseFlagVector<variableState> &states1 = states.getVec(scen);
		const denseVector &l1 = myLB().getVec(scen);
		const denseVector &u1 = myUB().getVec(scen);
		const denseVector &c1 = myObjective().getVec(scen);
		for (int i = 0; i < nvar; i++) {
			f << i;
			if (states1[i] == Basic) {
				f << " Basic";
			} else if (states1[i] == AtLower) {
				f << " AtLower";
			} else { 
				f << " AtUpper";
			}
			if (!writeBasisOnly) f << " " << tmpDse[i] << " " << l1[i] << " " << u1[i] << " " << c1[i];
			f << "\n";
		}
		f.close();
	}
}



void BALPSolverDual::loadStatus2(int scen, istream& f) {
	int nvar = data.dims.numVars(scen);
	int nbasic1;
	denseVector &dse1 = dse.getVec(scen);
	vector<int> &basicIdx1 = basicIdx.getVec(scen);
	basicIdx1.clear();
	denseFlagVector<variableState> &states1 = states.getVec(scen);
	denseVector &c1 = cPerturb.getVec(scen);

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
	int idx; double d1, d2, d3; string status;
	for (int i = 0; i < nvar; i++) {
		assert(!f.eof());
		f >> idx;
		f >> status;
		if (format == OldFormat) {
			f >> d1;
		} else if (format == BasisOnly) {
			d1 = 1.0;
		} else {
			f >> d1;
			f >> d2;
			f >> d3;
			f >> c1[idx];
		}
		if (status == "Basic") {
			states1[idx] = Basic;
			dse1[basicIdx1.size()] = max(d1,1e-4);
			basicIdx1.push_back(idx);
		} else if (status == "AtLower") {
			states1[idx] = AtLower;
		} else if (status == "AtUpper") {
			states1[idx] = AtUpper;
		}
	}
	assert(static_cast<int>(basicIdx1.size()) == nbasic1);
}

void BALPSolverDual::setSlackBasis() {

	const vector<int> &localScen = data.ctx.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		denseFlagVector<variableState> &s = states.getVec(scen);
		const denseFlagVector<constraintType> &t = data.vartype.getVec(scen);
		int nvarreal = data.dims.inner.numVars(scen);
		int ncons = data.dims.numCons(scen);
		for (int j = 0; j < nvarreal; j++) {
			if (t[j] != UB) {
				s[j] = AtLower;
			} else {
				s[j] = AtUpper;
			}
		}
		for (int j = 0; j < ncons; j++) {
			s[j+nvarreal] = Basic;
		}
	}
	status = LoadedFromFile;
	setupIndices();
}


