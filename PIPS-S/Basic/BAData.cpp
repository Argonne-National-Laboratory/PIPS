#include "BAData.hpp"
#include <cmath>

using namespace std;

template<typename T1, typename T2> static void mergeColAndRow(T1& v, const T2 &col, const T2 &row) {
	for (unsigned i = 0; i < col.size(); i++) {
		v[i] = col[i];
	}
	for (unsigned i = 0; i < row.size(); i++) {
		v[i+col.size()] = row[i];
	}

}

template<typename BAVec, typename Col1, typename Col2, typename Row1, typename Row2> 
	static void formBAVector(BAVec &v, Col1 c1, Col2 c2, Row1 r1, Row2 r2, const BAContext& ctx) {
	
	const vector<int> localScen = ctx.localScenarios();
	mergeColAndRow(v.getFirstStageVec(), c1(), r1());

	for (unsigned i = 1; i < localScen.size(); i++) {
		int scen = localScen[i];
		mergeColAndRow(v.getSecondStageVec(scen), c2(scen), r2(scen));
	}

}

// this mess is due to the magic and perveseness of C++
// try cleaning it up when we can use C++11
template<typename T> class wrapper1 {
public:
	wrapper1(stochasticInput& i, vector<T> (stochasticInput::*fun)()) : i(i), fun(fun) {}
	vector<T> operator()() { return (i.*fun)(); }

private:
	stochasticInput &i;
	vector<T> (stochasticInput::*fun)();
};

template<typename T> class wrapper2 {
public:
	wrapper2(stochasticInput& i, vector<T> (stochasticInput::*fun)(int)) : i(i), fun(fun) {}
	vector<T> operator()(int k) { return (i.*fun)(k); }

private:
	stochasticInput &i;
	vector<T> (stochasticInput::*fun)(int);
};

class emptyRowVector {
public:
	emptyRowVector(stochasticInput &i) : i(i) {}
	vector<double> operator()() { return vector<double>(i.nFirstStageCons(),0.0); }
	vector<double> operator()(int s) { return vector<double>(i.nSecondStageCons(s),0.0); }

private:
	stochasticInput &i;
};


BAData::BAData(stochasticInput &input, BAContext &ctx) : dims(BADimensions(input,ctx)), ctx(ctx) {

	int nscen = dims.numScenarios();
	ctx.initializeAssignment(nscen); // must do this first
	
	l.allocate(dims, ctx, PrimalVector);
	u.allocate(dims, ctx, PrimalVector);
	c.allocate(dims, ctx, PrimalVector);
	vartype.allocate(dims, ctx, PrimalVector);
	names.allocate(dims, ctx, PrimalVector);

	
	/*
	formBAVector(l, bind1st(mem_fun(&stochasticInput::getFirstStageColLB),&input),
			bind1st(mem_fun(&stochasticInput::getSecondStageColLB),&input),
			bind1st(mem_fun(&stochasticInput::getFirstStageRowLB),&input),
			bind1st(mem_fun(&stochasticInput::getSecondStageRowLB),&input),ctx);
	formBAVector(l, boost::bind(&stochasticInput::getFirstStageColLB,&input),
			boost::bind(&stochasticInput::getSecondStageColLB,&input),
			boost::bind(&stochasticInput::getFirstStageRowLB,&input),
			boost::bind(&stochasticInput::getSecondStageRowLB,&input),ctx);
	*/

	formBAVector(l, wrapper1<double>(input,&stochasticInput::getFirstStageColLB),
			wrapper2<double>(input,&stochasticInput::getSecondStageColLB),
			wrapper1<double>(input,&stochasticInput::getFirstStageRowLB),
			wrapper2<double>(input,&stochasticInput::getSecondStageRowLB),ctx);
	
	formBAVector(u, wrapper1<double>(input,&stochasticInput::getFirstStageColUB),
			wrapper2<double>(input,&stochasticInput::getSecondStageColUB),
			wrapper1<double>(input,&stochasticInput::getFirstStageRowUB),
			wrapper2<double>(input,&stochasticInput::getSecondStageRowUB),ctx);

	formBAVector(c, wrapper1<double>(input,&stochasticInput::getFirstStageObj),
			wrapper2<double>(input,&stochasticInput::getSecondStageObj),
			emptyRowVector(input), emptyRowVector(input), ctx);
	
	formBAVector(names, wrapper1<string>(input,&stochasticInput::getFirstStageColNames),
			wrapper2<string>(input,&stochasticInput::getSecondStageColNames),
			wrapper1<string>(input,&stochasticInput::getFirstStageRowNames),
			wrapper2<string>(input,&stochasticInput::getSecondStageRowNames),ctx);

	Acol = new CoinPackedMatrix(input.getFirstStageConstraints());
	
	/*
	Acol = data.getAmat();
	Arow = new CoinPackedMatrix();
	Arow->reverseOrderedCopyOf(*Acol);

	// add slacks to column copy
	// only used for INVERT
	int nslacks1 = nrows1;
	vector<double> slackvals1(nslacks1,-1);
	vector<int> start1(nslacks1), len1(nslacks1,1), ind1(nslacks1);
	for (int i = 0; i < nslacks1; i++) { 
		start1[i] = i;
		ind1[i] = i;
	}
	CoinPackedMatrix slackMat1(true,nslacks1,nslacks1,nslacks1,&slackvals1[0],&ind1[0],&start1[0],&len1[0]);
	if (appendSlacks) Acol->rightAppendPackedMatrix(slackMat1);

			
	len1.assign(nslacks1,0);
	start1.assign(nslacks1,0);

	Tcol.resize(nscen,0); Trow.resize(nscen,0);
	Wcol.resize(nscen,0); Wrow.resize(nscen,0);
	onlyBoundsVary = data.onlyBoundsVary();
	if (onlyBoundsVary) {
		Tcol[0] = data.getTmat(0);
		Trow[0] = new CoinPackedMatrix();
		Trow[0]->reverseOrderedCopyOf(*Tcol[0]);

		Wcol[0] = data.getWmat(0);
		Wrow[0] = new CoinPackedMatrix();
		Wrow[0]->reverseOrderedCopyOf(*Wcol[0]);

		int nslacks2 = dims.numSecondStageCons(0);
		vector<double> slackvals2(nslacks2,-1);
		vector<int> start2(nslacks2), len2(nslacks2,1), ind2(nslacks2);
		for (int i = 0; i < nslacks2; i++) { 
			start2[i] = i;
			ind2[i] = i;
		}
		CoinPackedMatrix slackMat2(true,nslacks2,nslacks2,nslacks2,&slackvals2[0],&ind2[0],&start2[0],&len2[0]);
		if (appendSlacks) Wcol[0]->rightAppendPackedMatrix(slackMat2);

		CoinPackedMatrix emptyMat(true,nslacks2,nslacks1,0,0,0,&start1[0],&len1[0]);
		if (appendSlacks) Tcol[0]->rightAppendPackedMatrix(emptyMat);

		for (int i = 1; i < nscen; i++) {
			Tcol[i] = Tcol[0];
			Trow[i] = Trow[0];
			Wcol[i] = Wcol[0];
			Wrow[i] = Wrow[0];
		}
	} else {
		for (int i = 0; i < nscen; i++) {
			if (!ctx.assignedScenario(i)) continue;
			Tcol[i] = data.getTmat(i);
			Trow[i] = new CoinPackedMatrix();
			Trow[i]->reverseOrderedCopyOf(*Tcol[i]);

			Wcol[i] = data.getWmat(i);
			Wrow[i] = new CoinPackedMatrix();
			Wrow[i]->reverseOrderedCopyOf(*Wcol[i]);

			int nslacks2 = dims.numSecondStageCons(i);
			vector<double> slackvals2(nslacks2,-1);
			vector<int> start2(nslacks2), len2(nslacks2,1), ind2(nslacks2);
			for (int j = 0; j < nslacks2; j++) { 
				start2[j] = j;
				ind2[j] = j;
			}
			CoinPackedMatrix slackMat2(true,nslacks2,nslacks2,nslacks2,&slackvals2[0],&ind2[0],&start2[0],&len2[0]);
			if (appendSlacks) Wcol[i]->rightAppendPackedMatrix(slackMat2);

			CoinPackedMatrix emptyMat(true,nslacks2,nslacks1,0,0,0,&start1[0],&len1[0]);
			if (appendSlacks) Tcol[i]->rightAppendPackedMatrix(emptyMat);
		}

	}

	slackbasis = !data.hasStartingBasis();
	if (!slackbasis) { 
		stateCol.allocate(dims.inner, ctx, PrimalVector);
		stateRow.allocate(dims.inner, ctx, DualVector);
		data.getStartingBasis(stateCol,stateRow);
	}
	out1Send.reserve(dims.numFirstStageVars());
	*/

}

BAData::BAData(const BAData &d) : dims(d.dims.inner), ctx(d.ctx) {

	l.allocate(dims, ctx, PrimalVector);
	u.allocate(dims, ctx, PrimalVector);
	c.allocate(dims, ctx, PrimalVector);
	vartype.allocate(dims, ctx, PrimalVector);
	names.allocate(dims, ctx, PrimalVector);

	slackbasis = d.slackbasis;
	onlyBoundsVary = d.onlyBoundsVary;

	// note that we make a SHALLOW copy of the constraint matrices
	// this means that d cannot be deleted while this instance is being used
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

	if (!slackbasis) { 
		stateCol.allocate(dims.inner, ctx, PrimalVector);
		stateRow.allocate(dims.inner, ctx, DualVector);
		stateCol.copyFrom(d.stateCol);
		stateRow.copyFrom(d.stateRow);
	}

	out1Send.reserve(dims.numFirstStageVars());


}

BAData::~BAData() {
	/*
	delete Acol;
	delete Arow;
	int nscen = dims.numScenarios();
	if (onlyBoundsVary) {
		delete Tcol[0];
		delete Trow[0];
		delete Wcol[0];
		delete Wrow[0];
	} else {
		for (int i = 0; i < nscen; i++) {
			if (Tcol[i]) {
				delete Tcol[i];
				delete Trow[i];
				delete Wcol[i];
				delete Wrow[i];
			}
		}
	}*/
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
		if (!ctx.assignedScenario(scen)) return;
		CoinIndexedVector &v2 = v.getSecondStageVec(scen).v;
		double *v2Elts = v2.denseVector();
		int *v2Idx = v2.getIndices();
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
	}
}

void BAData::addColToVec(sparseBAVector &v, BAIndex idx, double mult) const {
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

void BAData::getStartingBasis(BAFlagVector<variableState> &cols, bool &slackbasis) const {
	
	slackbasis = this->slackbasis;
	int nscen = dims.numScenarios();
	if (slackbasis) {
		for (int scen = -1; scen < nscen; scen++) {
			if (!ctx.assignedScenario(scen)) continue;
			denseFlagVector<variableState> &state = cols.getVec(scen);
			const denseFlagVector<constraintType> &vartype1 = vartype.getVec(scen);
			int ncols = dims.inner.numVars(scen);
			int ncolsSlacks = dims.numVars(scen);
			for (int i = 0; i < ncols; i++) {
				if (vartype1[i] != UB) {
					state[i] = AtLower;
				} else {
					state[i] = AtUpper;
				}
			}
			for (int i = ncols; i < ncolsSlacks; i++) {
				state[i] = Basic;
			}
		}
	}
	assert(slackbasis);
	/*else {
		// already filled stateCol and stateRow
	
		int nrows1 = dims.numFirstStageCons();
		int ncols1 = dims.inner.numFirstStageVars();
		mergeColAndRow(cols.getFirstStageVec(),stateCol.getFirstStageVec(),stateRow.getFirstStageVec(),nrows1,ncols1);

		for (int i = 0; i < nscen; i++) {
			if (!ctx.assignedScenario(i)) continue;
			int nrows2 = dims.numSecondStageCons(i);
			int ncols2 = dims.inner.numSecondStageVars(i);
			
			mergeColAndRow(cols.getSecondStageVec(i),stateCol.getSecondStageVec(i),stateRow.getSecondStageVec(i), nrows2, ncols2);
		}
	}*/

}
