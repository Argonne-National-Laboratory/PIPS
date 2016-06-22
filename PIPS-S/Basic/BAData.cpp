#include "BAData.hpp"
#include <cmath>
#include <boost/bind.hpp>
#include "PIPSLogging.hpp"
using boost::bind; // change to std::bind with C++11
using namespace std;

namespace {

template<typename T1, typename T2> void mergeColAndRow(T1& v, const T2 &col, const T2 &row) {
	for (unsigned i = 0; i < col.size(); i++) {
		v[i] = col[i];
	}
	for (unsigned i = 0; i < row.size(); i++) {
		v[i+col.size()] = row[i];
	}

}

template<typename BAVec, typename T1, typename Col2, typename Row2>
	void formBAVector(BAVec &v, const T1 &c1, Col2 c2, const T1 &r1, Row2 r2, const BAContext& ctx) {

	const vector<int> localScen = ctx.localScenarios();
	mergeColAndRow(v.getFirstStageVec(), c1, r1);

	for (unsigned i = 1; i < localScen.size(); i++) {
		int scen = localScen[i];
		mergeColAndRow(v.getSecondStageVec(scen), c2(scen), r2(scen));
	}

}


// use a lambda for this when we can use C++11
class emptyRowVector {
public:
	emptyRowVector(stochasticInput &i) : i(i) {}
	vector<double> operator()(int s) { return vector<double>(i.nSecondStageCons(s),0.0); }

private:
	stochasticInput &i;
};

void checkConstraintType(const denseVector &L, const denseVector &U, denseFlagVector<constraintType> &T) {
	int n = L.length();
	assert(U.length() == n && T.length() == n);

	for (int i = 0; i < n; i++) {
		if (L[i] < -1e20) {
			if (U[i] > 1e20) {
				T[i] = Free;
			} else {
				T[i] = UB;
			}
		} else {
			if (U[i] > 1e20) {
				T[i] = LB;
			} else {
				T[i] = Range;
			}
		}
		if (T[i] == Range && abs(L[i] - U[i]) < 1e-8) {
			T[i] = Fixed;
		}
	}

}

}

BAData::BAData(stochasticInput &input, BAContext &ctx) : ctx(ctx) {
	int nscen = input.nScenarios();
	ctx.initializeAssignment(nscen); // must do this first
	dims = BADimensions(input,ctx);

	const vector<int> &localScen = ctx.localScenarios();

	l.allocate(dims, ctx, PrimalVector);
	u.allocate(dims, ctx, PrimalVector);
	c.allocate(dims, ctx, PrimalVector);
	vartype.allocate(dims, ctx, PrimalVector);
	names.allocate(dims, ctx, PrimalVector);


	formBAVector(l, input.getFirstStageColLB(),
			bind(&stochasticInput::getSecondStageColLB,&input,_1),
			input.getFirstStageRowLB(),
			bind(&stochasticInput::getSecondStageRowLB,&input,_1),ctx);

	formBAVector(u, input.getFirstStageColUB(),
			bind(&stochasticInput::getSecondStageColUB,&input,_1),
			input.getFirstStageRowUB(),
			bind(&stochasticInput::getSecondStageRowUB,&input,_1),ctx);

	formBAVector(c, input.getFirstStageObj(),
			bind(&stochasticInput::getSecondStageObj,&input,_1),
			vector<double>(input.nFirstStageCons(),0.0),
			emptyRowVector(input), ctx);

	formBAVector(names, input.getFirstStageColNames(),
			bind(&stochasticInput::getSecondStageColNames,&input,_1),
			input.getFirstStageRowNames(),
			bind(&stochasticInput::getSecondStageRowNames,&input,_1),ctx);

	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		checkConstraintType(l.getVec(scen),u.getVec(scen),vartype.getVec(scen));
	}

	Acol.reset(new CoinPackedMatrix(input.getFirstStageConstraints()));
	Arow.reset(new CoinPackedMatrix());
	Arow->reverseOrderedCopyOf(*Acol);

	Tcol.resize(nscen); Trow.resize(nscen);
	Wcol.resize(nscen); Wrow.resize(nscen);

	/*
	We can save memory by not duplicating the constraint matrices, but
	as soon as we add individual scenario cuts these need to be duplicated anyway.
	One could think about adding identical cuts (with different RHSs)
	to each scenario to restore this optimization.

	onlyBoundsVary = input.onlyBoundsVary();
	if (onlyBoundsVary) {
		Tcol[0].reset(new CoinPackedMatrix(input.getLinkingConstraints(0)));
		Trow[0].reset(new CoinPackedMatrix());
		Trow[0]->reverseOrderedCopyOf(*Tcol[0]);

		Wcol[0].reset(new CoinPackedMatrix(input.getSecondStageConstraints(0)));
		Wrow[0].reset(new CoinPackedMatrix());
		Wrow[0]->reverseOrderedCopyOf(*Wcol[0]);

		for (int i = 1; i < nscen; i++) {
			Tcol[i] = Tcol[0];
			Trow[i] = Trow[0];
			Wcol[i] = Wcol[0];
			Wrow[i] = Wrow[0];
		}
	} else {*/
		for (int i = 0; i < nscen; i++) {
			if (!ctx.assignedScenario(i)) continue;
			Tcol[i].reset(new CoinPackedMatrix(input.getLinkingConstraints(i)));
			Trow[i].reset(new CoinPackedMatrix());
			Trow[i]->reverseOrderedCopyOf(*Tcol[i]);

			Wcol[i].reset(new CoinPackedMatrix(input.getSecondStageConstraints(i)));
			Wrow[i].reset(new CoinPackedMatrix());
			Wrow[i]->reverseOrderedCopyOf(*Wcol[i]);
		}
	/*
	}*/

	out1Send.reserve(dims.numFirstStageVars());

}

BAData::BAData(const BAData &d) : dims(d.dims.inner), ctx(d.ctx) {
	l.allocate(dims, ctx, PrimalVector);
	u.allocate(dims, ctx, PrimalVector);
	c.allocate(dims, ctx, PrimalVector);
	vartype.allocate(dims, ctx, PrimalVector);
	names.allocate(dims, ctx, PrimalVector);

	onlyBoundsVary = d.onlyBoundsVary;

	Acol = d.Acol;
	Arow = d.Arow;
	Tcol = d.Tcol;
	Trow = d.Trow;
	Wcol = d.Wcol;
	Wrow = d.Wrow;

	l.copyFrom(d.l);
	u.copyFrom(d.u);
	c.copyFrom(d.c);
	vartype.copyFrom(d.vartype);
	names.copyFrom(d.names);

	out1Send.reserve(dims.numFirstStageVars());


}

BAData::~BAData() {
	// shared_ptr handles the deletes
}


// assume v is cleared
void BAData::getCol(sparseBAVector &v, BAIndex idx) const {
	int scen = idx.scen;
	int col = idx.idx;
	const vector<int> &localScen = v.localScenarios();

	if (scen == -1) { // first stage
		CoinIndexedVector &v1 = v.getFirstStageVec().v;
		double *v1Elts = v1.denseVector();
		int *v1Idx = v1.getIndices();
		int nvarreal = dims.inner.numFirstStageVars();
		if (col < nvarreal) {
			int nnz = 0;
			const double *AcolElts = Acol->getElements();
			const int *AcolIdx = Acol->getIndices();
			CoinBigIndex start = Acol->getVectorFirst(col);
			CoinBigIndex end = Acol->getVectorLast(col);
			for (CoinBigIndex q = start; q < end; q++) {
				int row = AcolIdx[q];
				v1Elts[row] = AcolElts[q];
				v1Idx[nnz++] = row;
			}
			v1.setNumElements(nnz);

			for (unsigned j = 1; j < localScen.size(); j++) {
				int s = localScen[j];
				const double *TcolElts = Tcol[s]->getElements();
				const int *TcolIdx = Tcol[s]->getIndices();
				CoinIndexedVector &v2 = v.getSecondStageVec(s).v;
				double *v2Elts = v2.denseVector();
				int *v2Idx = v2.getIndices();
				nnz = 0;
				start = Tcol[s]->getVectorFirst(col);
				end = Tcol[s]->getVectorLast(col);
				for (CoinBigIndex q = start; q < end; q++) {
					int row = TcolIdx[q];
					v2Elts[row] = TcolElts[q];
					v2Idx[nnz++] = row;
				}
				v2.setNumElements(nnz);
			}
		} else {
			int k = col - nvarreal;
			v1Elts[k] = -1.0;
			v1Idx[0] = k;
			v1.setNumElements(1);
		}

	} else {
		if (!ctx.assignedScenario(scen)) return;
		CoinIndexedVector &v2 = v.getSecondStageVec(scen).v;
		double *v2Elts = v2.denseVector();
		int *v2Idx = v2.getIndices();
		int nvarreal = dims.inner.numSecondStageVars(scen);
		if (col < nvarreal) {
			const double *WcolElts = Wcol[scen]->getElements();
			const int *WcolIdx = Wcol[scen]->getIndices();
			int nnz = 0;
			CoinBigIndex start = Wcol[scen]->getVectorFirst(col);
			CoinBigIndex end = Wcol[scen]->getVectorLast(col);
			for (CoinBigIndex q = start; q < end; q++) {
				int row = WcolIdx[q];
				v2Elts[row] = WcolElts[q];
				v2Idx[nnz++] = row;
			}
			v2.setNumElements(nnz);
		} else {
			int k = col - nvarreal;
			v2Elts[k] = -1.0;
			v2Idx[0] = k;
			v2.setNumElements(1);
		}
	}
}

void BAData::addColToVec(sparseBAVector &v, BAIndex idx, double mult) const {
	assert(0); // need to fix handling of slacks
	int scen = idx.scen;
	int col = idx.idx;
	const vector<int> &localScen = v.localScenarios();

	if (scen == -1) { // first stage
		CoinIndexedVector &v1 = v.getFirstStageVec().v;
		const double *AcolElts = Acol->getElements();
		const int *AcolIdx = Acol->getIndices();
		CoinBigIndex start = Acol->getVectorFirst(col);
		CoinBigIndex end = Acol->getVectorLast(col);
		for (CoinBigIndex q = start; q < end; q++) {
			int row = AcolIdx[q];
			v1.quickAdd(row,mult*AcolElts[q]);
		}

		for (unsigned j = 1; j < localScen.size(); j++) {
			int s = localScen[j];
			const double *TcolElts = Tcol[s]->getElements();
			const int *TcolIdx = Tcol[s]->getIndices();
			CoinIndexedVector &v2 = v.getSecondStageVec(s).v;
			start = Tcol[s]->getVectorFirst(col);
			end = Tcol[s]->getVectorLast(col);
			for (CoinBigIndex q = start; q < end; q++) {
				int row = TcolIdx[q];
				v2.quickAdd(row,mult*TcolElts[q]);
			}
		}

	} else {
		if (!ctx.assignedScenario(scen)) return;
		CoinIndexedVector &v2 = v.getSecondStageVec(scen).v;
		const double *WcolElts = Wcol[scen]->getElements();
		const int *WcolIdx = Wcol[scen]->getIndices();
		CoinBigIndex start = Wcol[scen]->getVectorFirst(col);
		CoinBigIndex end = Wcol[scen]->getVectorLast(col);
		for (CoinBigIndex q = start; q < end; q++) {
			int row = WcolIdx[q];
			v2.quickAdd(row,mult*WcolElts[q]);
		}
	}



}

// multiplies "in" vector by constraint matrix
// set elements corresponding to basic columns to zero if you want to multiply with nonbasic columns
// uses dot product with rows approach, since "in" won't be hyper-sparse
void BAData::multiply(const denseBAVector &in, sparseBAVector &out, const BAFlagVector<variableState>& s) const {

	const denseVector &in1 = in.getFirstStageVec();
	CoinIndexedVector &out1 = out.getFirstStageVec().v;
	double *out1Elts = out1.denseVector();
	int *out1Idx = out1.getIndices();
	int nvar1Real = dims.inner.numFirstStageVars();
	int ncons1 = dims.numFirstStageCons();

	const vector<int> &localScen = in.localScenarios();
	for (unsigned j = 1; j < localScen.size(); j++) {
		int scen = localScen[j];

		const denseVector &in2 = in.getSecondStageVec(scen);
		CoinIndexedVector &out2 = out.getSecondStageVec(scen).v;
		double *out2Elts = out2.denseVector();
		int *out2Idx = out2.getIndices();

		int nvar2Real = dims.inner.numSecondStageVars(scen);
		int ncons2 = dims.numSecondStageCons(scen);
		int nnz2 = 0;

		const double *WrowElts = Wrow[scen]->getElements();
		const int *WrowIdx = Wrow[scen]->getIndices();
		const double *TrowElts = Trow[scen]->getElements();
		const int *TrowIdx = Trow[scen]->getIndices();


		for (int i = 0; i < ncons2; i++) {
			double work = 0.0;
			// W part
			CoinBigIndex start = Wrow[scen]->getVectorFirst(i);
			CoinBigIndex end = Wrow[scen]->getVectorLast(i);
			for (CoinBigIndex q = start; q < end; q++) {
				int col = WrowIdx[q];
				work += WrowElts[q]*in2[col];
			}
			int slackIdx = nvar2Real + i;
			work -= in2[slackIdx];
			// T part
			start = Trow[scen]->getVectorFirst(i);
			end = Trow[scen]->getVectorLast(i);
			for (CoinBigIndex q = start; q < end; q++) {
				int col = TrowIdx[q];
				work += TrowElts[q]*in1[col];
			}
			if (abs(work) > 1e-13) {
				out2Elts[i] = work;
				out2Idx[nnz2++] = i;
			} else {
				out2Elts[i] = 0.0;
			}
		}

		out2.setNumElements(nnz2);

	}

	const double *ArowElts = Arow->getElements();
	const int *ArowIdx = Arow->getIndices();
	int nnz1 = 0;

	// A block
	for (int i = 0; i < ncons1; i++) {
		CoinBigIndex start = Arow->getVectorFirst(i);
		CoinBigIndex end = Arow->getVectorLast(i);
		for (CoinBigIndex q = start; q < end; q++) {
			int col = ArowIdx[q];
			out1Elts[i] += ArowElts[q]*in1[col];
		}
		int slackIdx = nvar1Real + i;
		out1Elts[i] -= in1[slackIdx];
		if (abs(out1Elts[i]) > 1e-13) {
			out1Idx[nnz1++] = i;
		} else {
			out1Elts[i] = 0.0;
		}
	}
	out1.setNumElements(nnz1);

}

//#define PIPSPROF

void BAData::multiplyT(const sparseBAVector &in, sparseBAVector &out) const {
	// assume out is cleared

	const CoinIndexedVector &in1 = in.getFirstStageVec().v;
	CoinIndexedVector &out1 = out.getFirstStageVec().v;
	const double *in1Elts = in1.denseVector();
	const int* in1Idx = in1.getIndices();
	int nvarReal1 = dims.inner.numFirstStageVars();
	int nnzIn1 = in1.getNumElements();


	const vector<int> &localScen = in.localScenarios();

#ifdef PIPSPROF
	MPI_Barrier(ctx.comm());
	double local_t = MPI_Wtime();
#endif

	for (unsigned i = 1; i < localScen.size(); i++) {
		int scen = localScen[i];

		const CoinIndexedVector &in2 = in.getSecondStageVec(scen).v;
		CoinIndexedVector &out2 = out.getSecondStageVec(scen).v;
		const double *in2Elts = in2.denseVector();
		const int* in2Idx = in2.getIndices();

		int nvarReal2 = dims.inner.numSecondStageVars(scen);
		int nnzIn2 = in2.getNumElements();

		const double *WrowElts = Wrow[scen]->getElements();
		const int *WrowIdx = Wrow[scen]->getIndices();
		const double *TrowElts = Trow[scen]->getElements();
		const int *TrowIdx = Trow[scen]->getIndices();

		for (int j = 0; j < nnzIn2; j++) {
			int row = in2Idx[j];
			double mult = in2Elts[row];
			// W block -- add mult*(row) to out2
			CoinBigIndex start = Wrow[scen]->getVectorFirst(row);
			CoinBigIndex end = Wrow[scen]->getVectorLast(row);
			for (CoinBigIndex q = start; q < end; q++) {
				int col = WrowIdx[q];
				out2.quickAdd(col,mult*WrowElts[q]);
			}
			// slack, guaranteed to be zero in "out" until now,
			int slackIdx = nvarReal2 + row;
			out2.quickInsert(slackIdx,-mult);

			// T block -- goes into out1
			start = Trow[scen]->getVectorFirst(row);
			end = Trow[scen]->getVectorLast(row);
			for (CoinBigIndex q = start; q < end; q++) {
				int col = TrowIdx[q];
				out1Send.quickAdd(col,mult*TrowElts[q]);
			}
		}
	}
#ifdef PIPSPROF
	local_t = MPI_Wtime() - local_t;
	MPI_Barrier(ctx.comm());
	double comm_t = MPI_Wtime();
#endif

	MPI_Allreduce(out1Send.denseVector(), out1.denseVector(), dims.numFirstStageVars(), MPI_DOUBLE, MPI_SUM, ctx.comm());
#ifdef PIPSPROF
	comm_t = MPI_Wtime() - comm_t;
	double first_t = MPI_Wtime();
#endif

	out1Send.clear();

	// more efficient way to do it if T blocks are all the same

	// A block
	const double *ArowElts = Arow->getElements();
	const int *ArowIdx = Arow->getIndices();
	double *out1Elts = out1.denseVector();
	for (int j = 0; j < nnzIn1; j++) {
		int row = in1Idx[j];
		double mult = in1Elts[row];
		CoinBigIndex start = Arow->getVectorFirst(row);
		CoinBigIndex end = Arow->getVectorLast(row);
		for (CoinBigIndex q = start; q < end; q++) {
			int col = ArowIdx[q];
			out1Elts[col] += mult*ArowElts[q];
		}

		// slack
		int slackIdx = nvarReal1 + row;
		out1Elts[slackIdx] = -mult;

	}
	out1.scan(1e-13);

#ifdef PIPSPROF
	first_t = MPI_Wtime() - first_t;
	double totallocaltime;
	double maxlocaltime;
	MPI_Reduce(&local_t,&totallocaltime,1,MPI_DOUBLE,MPI_SUM,0,ctx.comm());
	MPI_Reduce(&local_t,&maxlocaltime,1,MPI_DOUBLE,MPI_MAX,0,ctx.comm());
	int nproc = ctx.nprocs();
	double idealtotal = totallocaltime/nproc + first_t;
	if (ctx.mype() == 0) {
		printf("PRICE ideal total: %g actual total: %g load imbalance: %f comm. cost: %f 1st stage cost: %f\n",
			idealtotal,maxlocaltime+comm_t + first_t,
			100.*(maxlocaltime-totallocaltime/nproc)/idealtotal,
			100.*(comm_t)/idealtotal,
			100.*(first_t)/idealtotal);
	}

#endif
}


void BAData::addSecondStageRow(const CoinPackedVectorBase& elts1, const CoinPackedVectorBase &elts2, int scen, double lb, double ub) {

	assert(scen >= 0 && scen < dims.numScenarios());
	if (ctx.assignedScenario(scen)) {
		Trow[scen]->appendRow(elts1);
		Tcol[scen]->reverseOrderedCopyOf(*Trow[scen]);
		Wrow[scen]->appendRow(elts2);
		Wcol[scen]->reverseOrderedCopyOf(*Wrow[scen]);

		int nvar2 = dims.inner.numSecondStageVars(scen);
		int ncons2 = dims.numSecondStageCons(scen);
		denseVector newL(nvar2+ncons2+1);
		denseVector &oldL = l.getSecondStageVec(scen);
		newL.copyBeginning(&oldL[0],nvar2+ncons2);
		newL[nvar2+ncons2] = lb;
		oldL.swap(newL);

		denseVector newU(nvar2+ncons2+1);
		denseVector &oldU = u.getSecondStageVec(scen);
		newU.copyBeginning(&oldU[0],nvar2+ncons2);
		newU[nvar2+ncons2] = ub;
		oldU.swap(newU);

		denseVector newC(nvar2+ncons2+1);
		denseVector &oldC = c.getSecondStageVec(scen);
		newC.copyBeginning(&oldC[0],nvar2+ncons2);
		newC[nvar2+ncons2] = 0.0;
		oldC.swap(newC);

		dims.inner.addSecondStageRow(scen);

		assert(oldU.length() == nvar2+ncons2+1);
		assert(oldL.length() == nvar2+ncons2+1);

	}

	vartype.deallocate();
	names.deallocate();

	vartype.allocate(dims, ctx, PrimalVector);

	// we're dropping the names here by not copying them
	// TODO: fix this
	names.allocate(dims, ctx, PrimalVector);

	const vector<int> &localScen = ctx.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		checkConstraintType(l.getVec(scen),u.getVec(scen),vartype.getVec(scen));
	}

}


	const CoinShallowPackedVector BAData::retrieveARow(int index)const{
		return Arow->getVector(index);

	}

	const CoinShallowPackedVector BAData::retrieveWRow(int index,int scen)const {
		return Wrow[scen]->getVector(index);

	}

	const CoinShallowPackedVector BAData::retrieveTRow(int index,int scen)const{
		return Trow[scen]->getVector(index);

	}

	const CoinShallowPackedVector BAData::retrieveACol(int index)const{
		return Acol->getVector(index);
	}

	const CoinShallowPackedVector BAData::retrieveWCol(int index,int scen)const{
		return Wcol[scen]->getVector(index);

	}

	const CoinShallowPackedVector BAData::retrieveTCol (int index,int scen)const{
		return Tcol[scen]->getVector(index);

	}

	void BAData::addSecondStageConsecutiveRows(const std::vector<CoinPackedVector*> &v1, const std::vector<CoinPackedVector*> &v2, int scenario, std::vector<double> &lb, std::vector<double> &ub, int nRows){


	CoinPackedVectorBase * const * elts1= (CoinPackedVectorBase * const *) &v1[0];
	CoinPackedVectorBase * const * elts2= (CoinPackedVectorBase * const *) &v2[0];
	assert(scenario >= 0 && scenario < dims.numScenarios());
	if (ctx.assignedScenario(scenario)) {
		Trow[scenario]->appendRows(nRows,elts1);
		Tcol[scenario]->reverseOrderedCopyOf(*Trow[scenario]);
		Wrow[scenario]->appendRows(nRows,elts2);
		Wcol[scenario]->reverseOrderedCopyOf(*Wrow[scenario]);

		int nvar2 = dims.inner.numSecondStageVars(scenario);
		int ncons2 = dims.numSecondStageCons(scenario);
		denseVector newL(nvar2+ncons2+nRows);
		denseVector &oldL = l.getSecondStageVec(scenario);
		newL.copyBeginning(&oldL[0],nvar2+ncons2);
		for(int i=0; i< nRows; i++)	newL[nvar2+ncons2+i] = lb[i];
		oldL.swap(newL);

		denseVector newU(nvar2+ncons2+nRows);
		denseVector &oldU = u.getSecondStageVec(scenario);
		newU.copyBeginning(&oldU[0],nvar2+ncons2);
		for(int i=0; i< nRows; i++) newU[nvar2+ncons2+i] = ub[i];
		oldU.swap(newU);

		denseVector newC(nvar2+ncons2+nRows);
		denseVector &oldC = c.getSecondStageVec(scenario);
		newC.copyBeginning(&oldC[0],nvar2+ncons2);
		for(int i=0; i< nRows; i++) newC[nvar2+ncons2+i] = 0.0;
		oldC.swap(newC);

		for(int i=0; i< nRows; i++) dims.inner.addSecondStageRow(scenario);

		assert(oldU.length() == nvar2+ncons2+nRows);
		assert(oldL.length() == nvar2+ncons2+nRows);

	}
	vartype.deallocate();
	names.deallocate();

	vartype.allocate(dims, ctx, PrimalVector);

	// we're dropping the names here by not copying them
	// TODO: fix this
	names.allocate(dims, ctx, PrimalVector);

	const vector<int> &localScen = ctx.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		checkConstraintType(l.getVec(scen),u.getVec(scen),vartype.getVec(scen));
	}
}


void BAData::addFirstStageRow(const CoinPackedVectorBase& elts1, double lb, double ub) {

	assert(lb<=ub);
	Arow->appendRow(elts1);
	Acol->reverseOrderedCopyOf(*Arow);
	int nvar = dims.inner.numFirstStageVars();
	int ncons = dims.numFirstStageCons();
	denseVector newL(nvar+ncons+1);
	denseVector &oldL = l.getFirstStageVec();
	newL.copyBeginning(&oldL[0],nvar+ncons);
	newL[nvar+ncons] = lb;
	oldL.swap(newL);
	denseVector newU(nvar+ncons+1);
	denseVector &oldU = u.getFirstStageVec();
	newU.copyBeginning(&oldU[0],nvar+ncons);
	newU[nvar+ncons] = ub;
	oldU.swap(newU);
	denseVector newC(nvar+ncons+1);
	denseVector &oldC = c.getFirstStageVec();
	newC.copyBeginning(&oldC[0],nvar+ncons);
	newC[nvar+ncons] = 0.0;
	oldC.swap(newC);
	dims.inner.addFirstStageRow();
	assert(oldU.length() == nvar+ncons+1);
	assert(oldL.length() == nvar+ncons+1);
	vartype.deallocate();
	names.deallocate();
	vartype.allocate(dims, ctx, PrimalVector);

	// we're dropping the names here by not copying them
	// TODO: fix this
	names.allocate(dims, ctx, PrimalVector);
	const vector<int> &localScen = ctx.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		checkConstraintType(l.getVec(scen),u.getVec(scen),vartype.getVec(scen));
	}
	out1Send.reserve(dims.numFirstStageVars());

}

void BAData::addFirstStageRows(const std::vector<CoinPackedVector*> &v1, std::vector<double> &lb, std::vector<double> &ub, int nRows) {


	CoinPackedVectorBase * const * elts1= (CoinPackedVectorBase * const *) &v1[0];
	assert(lb.size()==ub.size() && lb.size()==nRows && nRows>0);
	Arow->appendRows(nRows,elts1);
	Acol->reverseOrderedCopyOf(*Arow);

	int nvar2 = dims.inner.numFirstStageVars();
	int ncons2 = dims.numFirstStageCons();
	denseVector newL(nvar2+ncons2+nRows);
	denseVector &oldL = l.getFirstStageVec();
	newL.copyBeginning(&oldL[0],nvar2+ncons2);
	for(int i=0; i< nRows; i++) newL[nvar2+ncons2+i] = lb[i];
	oldL.swap(newL);
	denseVector newU(nvar2+ncons2+nRows);
	denseVector &oldU = u.getFirstStageVec();
	newU.copyBeginning(&oldU[0],nvar2+ncons2);
	for(int i=0; i< nRows; i++)newU[nvar2+ncons2+i] = ub[i];
	oldU.swap(newU);
	denseVector newC(nvar2+ncons2+nRows);
	denseVector &oldC = c.getFirstStageVec();
	newC.copyBeginning(&oldC[0],nvar2+ncons2);
	for(int i=0; i< nRows; i++)newC[nvar2+ncons2+i] = 0.0;
	oldC.swap(newC);

	for(int i=0; i< nRows; i++)dims.inner.addFirstStageRow();

	assert(oldU.length() == nvar2+ncons2+nRows);
	assert(oldL.length() == nvar2+ncons2+nRows);
	vartype.deallocate();
	names.deallocate();
	vartype.allocate(dims, ctx, PrimalVector);

	// we're dropping the names here by not copying them
	// TODO: fix this
	names.allocate(dims, ctx, PrimalVector);
	const vector<int> &localScen = ctx.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		checkConstraintType(l.getVec(scen),u.getVec(scen),vartype.getVec(scen));
	}
	out1Send.reserve(dims.numFirstStageVars());

}

int BAData::addFirstStageColumn( double lb, double ub, double cobj){

	vector<double> elems(Acol->getMinorDim(),0);
	CoinPackedVector elts;
	elts.setFullNonZero(elems.size(),&elems[0]);

	//Assertions
	assert(lb<=ub);

	Acol->appendCol(elts);
	Arow->reverseOrderedCopyOf(*Acol);

	for (int scen=0; scen < Tcol.size(); scen++){
		if (!ctx.assignedScenario(scen)) continue;
		Tcol[scen]->appendCol(elts);
		Trow[scen]->reverseOrderedCopyOf(*Tcol[scen]);
	}


	//increase the size of l, u, c,
	int nvar2 = dims.inner.numFirstStageVars();
	int ncons2 = dims.numFirstStageCons();
	denseVector newL(nvar2+ncons2+1);
	denseVector &oldL = l.getFirstStageVec();
	newL.copyBeginning(&oldL[0],nvar2);
	newL.copyToPosition(&oldL[nvar2],nvar2+1,ncons2);
	newL[nvar2] = lb;
	oldL.swap(newL);

	denseVector newU(nvar2+ncons2+1);
	denseVector &oldU = u.getFirstStageVec();
	newU.copyBeginning(&oldU[0],nvar2);
	newU.copyToPosition(&oldU[nvar2],nvar2+1,ncons2);
	newU[nvar2] = ub;
	oldU.swap(newU);

	denseVector newC(nvar2+ncons2+1);
	denseVector &oldC = c.getFirstStageVec();

	newC.copyBeginning(&oldC[0],nvar2);
	newC.copyToPosition(&oldC[nvar2],nvar2+1,ncons2);
	newC[nvar2] = cobj;
	oldC.swap(newC);

	assert(oldU.length() == nvar2+ncons2+1);
	assert(oldL.length() == nvar2+ncons2+1);

	//Update dims
	dims.inner.addFirstStageVar();

	//check constraint type?
	vartype.deallocate();
	names.deallocate();

	vartype.allocate(dims, ctx, PrimalVector);
	// we're dropping the names here by not copying them
	// TODO: fix this
	names.allocate(dims, ctx, PrimalVector);

	assert(l.getFirstStageVec().length()==u.getFirstStageVec().length());

	const vector<int> &localScen = ctx.localScenarios();

	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		checkConstraintType(l.getVec(scen),u.getVec(scen),vartype.getVec(scen));
	}
	out1Send.reserve(dims.numFirstStageVars());

	return nvar2;
}

int BAData::addSecondStageColumn(int scen,double lb, double ub, double cobj){

	assert(scen >= 0 && scen < dims.numScenarios());
	int returnIndex= -1;
	if (ctx.assignedScenario(scen)) {

		vector<double> elems(Wcol[scen]->getMinorDim(),0);
		CoinPackedVector elts;
		elts.setFullNonZero(elems.size(),&elems[0]);

		//Assertions
		assert(lb<=ub);

		Wcol[scen]->appendCol(elts);
		Wrow[scen]->reverseOrderedCopyOf(*Wcol[scen]);

		int nvar2 = dims.inner.numSecondStageVars(scen);
		int ncons2 = dims.numSecondStageCons(scen);
		denseVector newL(nvar2+ncons2+1);
		denseVector &oldL = l.getSecondStageVec(scen);
		newL.copyBeginning(&oldL[0],nvar2);
		newL.copyToPosition(&oldL[nvar2],nvar2+1,ncons2);
		newL[nvar2] = lb;
		oldL.swap(newL);

		denseVector newU(nvar2+ncons2+1);
		denseVector &oldU = u.getSecondStageVec(scen);
		newU.copyBeginning(&oldU[0],nvar2);
		newU.copyToPosition(&oldU[nvar2],nvar2+1,ncons2);
		newU[nvar2] = ub;
		oldU.swap(newU);

		denseVector newC(nvar2+ncons2+1);
		denseVector &oldC = c.getSecondStageVec(scen);
		newC.copyBeginning(&oldC[0],nvar2);
		newC.copyToPosition(&oldC[nvar2],nvar2+1,ncons2);
		newC[nvar2] = cobj;
		oldC.swap(newC);

		dims.inner.addSecondStageVar(scen);

		assert(oldU.length() == nvar2+ncons2+1);
		assert(oldL.length() == nvar2+ncons2+1);
		returnIndex=nvar2;
	}

	vartype.deallocate();
	names.deallocate();

	vartype.allocate(dims, ctx, PrimalVector);

	// we're dropping the names here by not copying them
	// TODO: fix this
	names.allocate(dims, ctx, PrimalVector);

	const vector<int> &localScen = ctx.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		checkConstraintType(l.getVec(scen),u.getVec(scen),vartype.getVec(scen));
	}
	return returnIndex;
}


void BAData::deleteLastFirstStageRows(int nRows) {

	assert(nRows>=0 && nRows<dims.numFirstStageCons());
	int rowStart=dims.numFirstStageCons()-nRows;
	vector<int> indices(nRows);
	for (int i=0; i<indices.size();i++)indices[i]=rowStart+i;
	Arow->deleteRows(nRows,&indices[0]);
	Acol->reverseOrderedCopyOf(*Arow);

	int nvar2 = dims.inner.numFirstStageVars();
	int ncons2 = dims.numFirstStageCons();
	denseVector newL(nvar2+ncons2-nRows);
	denseVector &oldL = l.getFirstStageVec();
	newL.copyBeginning(&oldL[0],nvar2+rowStart);
	oldL.swap(newL);

	denseVector newU(nvar2+ncons2-nRows);
	denseVector &oldU = u.getFirstStageVec();
	newU.copyBeginning(&oldU[0],nvar2+rowStart);
	oldU.swap(newU);

	denseVector newC(nvar2+ncons2-nRows);
	denseVector &oldC = c.getFirstStageVec();
	newC.copyBeginning(&oldC[0],nvar2+rowStart);
	oldC.swap(newC);

	for(int i=0; i< nRows; i++)dims.inner.removeFirstStageRow();

	assert(oldU.length() == nvar2+rowStart);
	assert(oldL.length() == nvar2+rowStart);

	vartype.deallocate();
	names.deallocate();

	vartype.allocate(dims, ctx, PrimalVector);

	// we're dropping the names here by not copying them
	// TODO: fix this
	names.allocate(dims, ctx, PrimalVector);

	const vector<int> &localScen = ctx.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		checkConstraintType(l.getVec(scen),u.getVec(scen),vartype.getVec(scen));
	}
	out1Send.reserve(dims.numFirstStageVars());

}

void BAData::deleteLastSecondStageConsecutiveRows(int scenario, int nRows){

	assert(scenario >= 0 && scenario < dims.numScenarios());
	if (ctx.assignedScenario(scenario)) {
		assert(nRows>=0 && nRows<dims.numSecondStageCons(scenario));

		int rowStart=dims.numSecondStageCons(scenario)-nRows;

		vector<int> indices(nRows);
		for (int i=0; i<indices.size();i++)indices[i]=rowStart+i;

		Trow[scenario]->deleteRows(nRows,&indices[0]);
		Tcol[scenario]->reverseOrderedCopyOf(*Trow[scenario]);
		Wrow[scenario]->deleteRows(nRows,&indices[0]);
		Wcol[scenario]->reverseOrderedCopyOf(*Wrow[scenario]);

		int nvar2 = dims.inner.numSecondStageVars(scenario);
		int ncons2 = dims.numSecondStageCons(scenario);
		denseVector newL(nvar2+ncons2-nRows);
		denseVector &oldL = l.getSecondStageVec(scenario);
		newL.copyBeginning(&oldL[0],nvar2+rowStart);
		oldL.swap(newL);

		denseVector newU(nvar2+ncons2-nRows);
		denseVector &oldU = u.getSecondStageVec(scenario);
		newU.copyBeginning(&oldU[0],nvar2+rowStart);
		oldU.swap(newU);

		denseVector newC(nvar2+ncons2-nRows);
		denseVector &oldC = c.getSecondStageVec(scenario);
		newC.copyBeginning(&oldC[0],nvar2+rowStart);
		oldC.swap(newC);

		for(int i=0; i< nRows; i++) dims.inner.removeSecondStageRow(scenario);
		assert(oldU.length() == nvar2+rowStart);
		assert(oldL.length() == nvar2+rowStart);
	}

	vartype.deallocate();
	names.deallocate();

	vartype.allocate(dims, ctx, PrimalVector);

	// we're dropping the names here by not copying them
	// TODO: fix this
	names.allocate(dims, ctx, PrimalVector);

	const vector<int> &localScen = ctx.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		checkConstraintType(l.getVec(scen),u.getVec(scen),vartype.getVec(scen));
	}

}


void BAData::deleteLastFirstStageColumns(int nCols){

	int ncons = dims.numFirstStageCons();
	int nvar = dims.inner.numFirstStageVars();
	//TODO add guard to make sure we don't delete a column that still has nonzero coefficients in the constraint matrix
	assert(nCols<=nvar);
	int colStart=dims.inner.numFirstStageVars()-nCols;
	vector<int> indices(nCols);
	for (int i=0; i<indices.size();i++)indices[i]=colStart+i;
	Acol->deleteCols(nCols,&indices[0]);
	Arow->reverseOrderedCopyOf(*Acol);

	for (int scen=0; scen < Tcol.size(); scen++){
		if (!ctx.assignedScenario(scen)) continue;
		Tcol[scen]->deleteCols(nCols,&indices[0]);
		Trow[scen]->reverseOrderedCopyOf(*Tcol[scen]);
	}

	//Shrink l, u, and c
	denseVector newL(nvar+ncons-nCols);
	denseVector &oldL = l.getFirstStageVec();
	newL.copyBeginning(&oldL[0],nvar-nCols);
	newL.copyToPosition(&oldL[nvar],nvar-nCols,ncons);
	oldL.swap(newL);

	denseVector newU(nvar+ncons-nCols);
	denseVector &oldU = u.getFirstStageVec();
	newU.copyBeginning(&oldU[0],nvar-nCols);
	newU.copyToPosition(&oldU[nvar],nvar-nCols,ncons);
	oldU.swap(newU);

	denseVector newC(nvar+ncons-nCols);
	denseVector &oldC = c.getFirstStageVec();
	newC.copyBeginning(&oldC[0],nvar-nCols);
	newC.copyToPosition(&oldC[nvar],nvar-nCols,ncons);
	oldC.swap(newC);

	assert(oldU.length() == nvar+ncons-nCols);
	assert(oldL.length() == nvar+ncons-nCols);

	//Update dims
	for (int i=0; i< nCols; i++)dims.inner.removeFirstStageVar();

	//check constraint type?
	vartype.deallocate();
	names.deallocate();

	vartype.allocate(dims, ctx, PrimalVector);
	// we're dropping the names here by not copying them
	// TODO: fix this
	names.allocate(dims, ctx, PrimalVector);

	assert(l.getFirstStageVec().length()==u.getFirstStageVec().length());

	const vector<int> &localScen = ctx.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		checkConstraintType(l.getVec(scen),u.getVec(scen),vartype.getVec(scen));
	}

	out1Send.reserve(dims.numFirstStageVars());

}



void BAData::deleteLastSecondStageConsecutiveColumns(int scen, int nCols){

	assert(scen >= 0 && scen < dims.numScenarios());

	if (ctx.assignedScenario(scen)) {

		int nvar2 = dims.inner.numSecondStageVars(scen);
		int ncons2 = dims.numSecondStageCons(scen);
		//TODO add guard to make sure we don't delete a column that still has nonzero coefficients in the constraint matrix
		assert(nCols<=nvar2);

		int colStart=dims.inner.numSecondStageVars(scen)-nCols;

		vector<int> indices(nCols);
		for (int i=0; i<indices.size();i++)indices[i]=colStart+i;

		Wcol[scen]->deleteCols(nCols,&indices[0]);
		Wrow[scen]->reverseOrderedCopyOf(*Wcol[scen]);

		//Shrink l, u, and c
		denseVector newL(nvar2+ncons2-nCols);
		denseVector &oldL = l.getSecondStageVec(scen);
		newL.copyBeginning(&oldL[0],nvar2-nCols);
		newL.copyToPosition(&oldL[nvar2],nvar2-nCols,ncons2);
		oldL.swap(newL);

		denseVector newU(nvar2+ncons2-nCols);
		denseVector &oldU = u.getSecondStageVec(scen);
		newU.copyBeginning(&oldU[0],nvar2-nCols);
		newU.copyToPosition(&oldU[nvar2],nvar2-nCols,ncons2);
		oldU.swap(newU);

		denseVector newC(nvar2+ncons2-nCols);
		denseVector &oldC = c.getSecondStageVec(scen);
		newC.copyBeginning(&oldC[0],nvar2-nCols);
		newC.copyToPosition(&oldC[nvar2],nvar2-nCols,ncons2);
		oldC.swap(newC);

		assert(oldU.length() == nvar2+ncons2-nCols);
		assert(oldL.length() == nvar2+ncons2-nCols);

		//Update dims
		for (int i=0; i< nCols; i++)dims.inner.removeSecondStageVar(scen);

	}
	//check constraint type?
	vartype.deallocate();
	names.deallocate();

	vartype.allocate(dims, ctx, PrimalVector);
	// we're dropping the names here by not copying them
	// TODO: fix this
	names.allocate(dims, ctx, PrimalVector);

	assert(l.getSecondStageVec(scen).length()==u.getSecondStageVec(scen).length());

	const vector<int> &localScen = ctx.localScenarios();
	for (unsigned i = 0; i < localScen.size(); i++) {
		int scen = localScen[i];
		checkConstraintType(l.getVec(scen),u.getVec(scen),vartype.getVec(scen));
	}

}

