#include "BAPFIPar.hpp"
#include <algorithm>

using namespace std;

BAPFIPar::BAPFIPar(const BAData &d, int max_updates) : data(d), npfi(0), maxUpdates(max_updates) {
	nScen = d.dims.numScenarios();
	enter.resize(max_updates);
	leave.resize(max_updates);
	etas.resize(max_updates);
	pivotEntries = new double[max_updates*max_updates];
	pivotRow = new double[max_updates];
	std::fill(pivotEntries,pivotEntries+maxUpdates*maxUpdates,0.0);
	
	pivotRhs.reserve(max_updates);
	pivotEntriesIdx.resize(max_updates);
	for (int i = 0; i < max_updates; i++) {
		pivotEntriesIdx[i].reserve(max_updates/5);
		etas[i].allocate(data.dims,data.ctx);
	}
	localLeaving.resize(data.ctx.nprocs());
	maxLocalLeaving = 0;

	const vector<int> &localScen = data.ctx.localScenarios();
	scenToLocalIdx.resize(nScen+1,-1);
	int idx = 0;
	for (unsigned i = 0; i < localScen.size(); i++) {
		scenToLocalIdx.at(localScen[i]+1) = idx++;
	}

}

BAPFIPar::~BAPFIPar() {
	delete [] pivotEntries;
	delete [] pivotRow;
}


void BAPFIPar::ftranPFI(sparseBAVector &rhs) {
	assert(rhs.vectorType() == BasicVector);
	
	int nprocs = data.ctx.nprocs();
	int mype = data.ctx.mype();
	double *pivotElts = pivotRhs.denseVector();
	
	/*for (int i = 0; i < npfi; i++) {
		pivotElts[i] = rhs[leave[i]];
	}*/

	// gather leaving elements from rhs

	// pick out first-stage elements (which we already have)
	for (unsigned i = 0; i < firstStageLeaving.size(); i++) {
		int idx = firstStageLeaving[i];
		pivotElts[idx] = rhs[leave[idx]];
	}
	// pack elements which we own, using pivotRow as the send buffer
	double *send = pivotRow;
	for (unsigned i = 0; i < localLeaving[mype].size(); i++) {
		int idx = localLeaving[mype][i];
		send[i] = rhs[leave[idx]];
	}
	ftranRecv.resizeAndDestroy(nprocs*maxLocalLeaving);
	double *recv = ftranRecv.ptr();
	
	MPI_Allgather(send,maxLocalLeaving,MPI_DOUBLE,recv, maxLocalLeaving, MPI_DOUBLE, data.ctx.comm());
	
	// unpack
	// don't worry about setting up indices for pivotRhs b/c they're not used (yet)
	for (int p = 0; p < nprocs; p++) {
		for (unsigned i = 0; i < localLeaving[p].size(); i++) {
			int idx = localLeaving[p][i];
			pivotElts[idx] = recv[p*maxLocalLeaving+i];
		}
	}

	miniFtranSparsish(pivotRhs);

	// now take linear combination of the etas with the local data	
	const vector<int> &localScen = rhs.localScenarios();
	const int *pivotIdx = pivotRhs.getIndices();
	int pivotNnz = pivotRhs.getNumElements();
	for (int j = 0; j < pivotNnz; j++) {
		int i = pivotIdx[j];
		double pivot = -pivotElts[i];
		assert(pivot);
		for (unsigned k = 0; k < localScen.size(); k++) {
			int scen = localScen[k];
			CoinIndexedVector &rhsVec = rhs.getVec(scen).v;
			const packedVector &eta = etas[i].getPackedVec(k);
			const double * COIN_RESTRICT etaElts = eta.denseVector();
			const int * COIN_RESTRICT etaIdx = eta.getIndices();
			int etaNnz = eta.getNumElements();
			for (int r = 0; r < etaNnz; r++) {
				rhsVec.quickAdd(etaIdx[r],pivot*etaElts[r]);
			}
		}
	}

	for (int i = 0; i < npfi; i++) {
		if (data.ctx.assignedScenario(leave[i].scen)) rhs.getVec(leave[i].scen).v.zero(leave[i].idx);
	}
	/*
	for (unsigned k = 0; k < localScen.size(); k++) {
		int scen = localScen[k];
		CoinIndexedVector &rhsVec = rhs.getVec(scen).v;
		printf("FTRAN from scen %d has %d nonzeros\n",scen,rhsVec.getNumElements());
	}*/


}



// fill-in occurs only in leaving entries
void BAPFIPar::btranPFISimplex(sparseBAVector &rhs) {
	assert(leaving == leave[npfi]);

	double * COIN_RESTRICT r = pivotRhs.denseVector();
	int * COIN_RESTRICT rIdx = pivotRhs.getIndices();
	memset(r,0,(npfi+1)*sizeof(double));
	r[npfi] = 1.0;
	rIdx[0] = npfi;
	int nnz = 1;
	
	// here we perform dot products with columns by looping
	// through nonzeros in rhs
	// this can also be written as a dense triangular solve
	for (int i = npfi-1; i >= 0; i--) {
		assert(r[i] == 0.0);
		double d = 0.0;
		for (int k = 0; k < nnz; k++) {
			int idx = rIdx[k];
			//printf("update: idx: (%d,%d) d: %.10e r: %.10e entry: %.10e\n",leave[idx].scen,leave[idx].idx,d,r[idx],pivotEntries[maxUpdates*i+idx]);
			d += r[idx]*pivotEntries[maxUpdates*i+idx];
		}
		if (fabs(d) > 1e-13) {
			if (data.ctx.assignedScenario(leave[i].scen))
				rhs.getVec(leave[i].scen).v.quickAddNonZero(leave[i].idx,-d);
			r[i] = -d;
			rIdx[nnz++] = i;
			//printf("BTRAN-PFI product %d: %e (%d,%d)\n",i,-d,leave[i].scen,leave[i].idx);
		}
		if (enter[i].scen == -1) {
			// because we don't loop through first-stage nonzeros in BALinearAlgebra::btran,
			// need to set to zero here
			rhs.getVec(enter[i].scen).v.denseVector()[enter[i].idx] = 0.0;
		} else if (data.ctx.assignedScenario(enter[i].scen)) {
			rhs.getVec(enter[i].scen).v.zero(enter[i].idx);
		}
	}
	// leave pivotRhs messy

}

extern "C" void dtrsv_(const char*UPLO,const char*TRANS,const char*DIAG,const int*N,const double*A,const int*LDA,double*RHS,int*INCRHS);

// because new columns come in at new indices, never old indices, an eta vector will never have 
// nonzero entries corresponding to variables that already left
// this gives us a triangular structure in the pivotal entries of the eta's,
// hence, we can actually apply the etas as a triangular solve.
// questionable whether this is faster or not
void BAPFIPar::miniFtranDense(CoinIndexedVector &rhs) const {
	
	double *r = rhs.denseVector();
	int one = 1;

	dtrsv_("L","N","U",&npfi,pivotEntries,&maxUpdates,r,&one);
	
	int *idx = rhs.getIndices();
	int nnz = 0;
	for (int i = 0; i < npfi; i++) {
		if (fabs(r[i]) > 1e-13) {
			idx[nnz++] = i;
		} else {
			r[i] = 0.0;
		}
	}
	rhs.setNumElements(nnz);

}

// exploit sparsity in pivot entries (but not in rhs)
void BAPFIPar::miniFtranSparsish(CoinIndexedVector &rhs) const {

	double * COIN_RESTRICT r = rhs.denseVector();
	int * COIN_RESTRICT rIdx = rhs.getIndices();
	int nnz = 0;

	for (int i = 0; i < npfi; i++) {
		const int* colIdx = &pivotEntriesIdx[i][0];
		int colNnz = pivotEntriesIdx[i].size();
		double pivot = r[i];
		if (fabs(pivot) > 1e-13) {
			rIdx[nnz++] = i;
		} else {
			r[i] = 0.0;
			continue;
		}
		for (int k = 0; k < colNnz; k++) {
			int row = colIdx[k];
			assert(row > i);
			// note negative sign
			r[row] -= pivot*pivotEntries[maxUpdates*i+row];
		}
	}
	rhs.setNumElements(nnz);
}



void BAPFIPar::newEta(sparseBAVector &ftranVec, BAIndex in, BAIndex out) {
	assert(out == leaving);
	enter[npfi] = in;

	double pivotInv;
	if (data.ctx.assignedScenario(out.scen)) {
		double etaPivot = ftranVec[out];
		ftranVec.getVec(out.scen).v.zero(out.idx);
		pivotInv = 1.0/etaPivot;
	}
	if (out.scen != -1) MPI_Bcast(&pivotInv,1,MPI_DOUBLE,data.ctx.owner(out.scen),data.ctx.comm()); 

	// so that "etaize" works correctly
	if (data.ctx.assignedScenario(in.scen)) ftranVec.getVec(in.scen).v.insert(in.idx,-1.0); 

	const vector<int> &localScen = ftranVec.localScenarios();
	// order indices and "etaize" local vectors
	for (unsigned k = 0; k < localScen.size(); k++) {
		int scen = localScen[k];
		CoinIndexedVector &vec = ftranVec.getVec(scen).v;
		const double * COIN_RESTRICT elts = vec.denseVector();
		int ftranNnz = vec.getNumElements();
		int nvar = vec.capacity();
		packedVector &packedVec = etas[npfi].getPackedVec(k);
		packedVec.resizeAndDestroy(ftranNnz);

		int * COIN_RESTRICT idx = packedVec.getIndices();
		double * COIN_RESTRICT packedElts = packedVec.denseVector();
		int nnz = 0;
		if (ftranNnz < 0.03*nvar) { // sort indices
			int *ftranIdx = vec.getIndices();
			std::sort(ftranIdx,ftranIdx+ftranNnz);
			for (int j = 0; j < ftranNnz; j++) {
				int i = ftranIdx[j];
				double val = elts[i];
				if (val != COIN_INDEXED_REALLY_TINY_ELEMENT) {
					packedElts[nnz] = val*pivotInv;
					idx[nnz++] = i;
				}
			}
		} else { // pass through full vector
			for (int i = 0; i < nvar; i++) {
				double val = elts[i];
				if (val && val != COIN_INDEXED_REALLY_TINY_ELEMENT) {
					packedElts[nnz] = val*pivotInv; // note: didn't negate
					idx[nnz++] = i;
				}
			}
		}
		packedVec.setNumElements(nnz);

	}
	

	// there are no leaving elements to pick out yet
	npfi++;
	assert(npfi <= maxUpdates);
}


void BAPFIPar::setLeaving(BAIndex leaving) {
	this->leaving = leaving;
	leave[npfi] = leaving;

	int owner = data.ctx.owner(leaving.scen);
	if (data.ctx.assignedScenario(leaving.scen)) {
		for (int i = 0; i < npfi; i++) {
			const packedVector &eta = etas[i].getPackedVec(scenToLocalIdx[leaving.scen+1]);
			const double *elts = eta.denseVector();
			const int *idx = eta.getIndices();
			int nnz = eta.getNumElements();
			if (!nnz) { pivotRow[i] = 0.0; continue; }
			// do a binary search through the sorted indices to find the match
			const int *idxPtr = lower_bound(idx,idx+nnz,leaving.idx);
			if (idxPtr != idx+nnz && *idxPtr == leaving.idx) { // has the element
				int pos = static_cast<int>(idxPtr - idx);
				assert(idx[pos] == leaving.idx);
				pivotRow[i] = elts[pos];
			} else {
				pivotRow[i] = 0.0;
			}
			
		}
	}
	// broadcast row of pivot matrix
	// could broadcast packed form, but would require extra communication
	if (leaving.scen != -1) MPI_Bcast(pivotRow,npfi,MPI_DOUBLE,owner,data.ctx.comm());

	// update pivot matrix with new row
	for (int i = 0; i < npfi; i++) {
		if (pivotRow[i]) {
			pivotEntries[maxUpdates*i+npfi] = pivotRow[i];
			pivotEntriesIdx[i].push_back(npfi);
		}
	}

	if (leaving.scen != -1) {
		localLeaving[owner].push_back(npfi);
		if (localLeaving[owner].size() > maxLocalLeaving) {
			maxLocalLeaving = localLeaving[owner].size();
		}
	} else {
		firstStageLeaving.push_back(npfi);
	}
	
}

void BAPFIPar::clear() {
	for (int i = 0; i < npfi; i++) {
		pivotEntriesIdx[i].clear();
		etas[i].clear();
	}
	std::fill(pivotEntries,pivotEntries+maxUpdates*maxUpdates,0.0);
	npfi = 0;
	
	firstStageLeaving.clear();
	for (int i = 0; i < data.ctx.nprocs(); i++) {
		localLeaving[i].clear();
	}
	maxLocalLeaving = 0;
}
