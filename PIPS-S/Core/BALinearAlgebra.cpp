#include "BALinearAlgebra.hpp"

using namespace std;

//#define PIPSPROF

BALinearAlgebra::BALinearAlgebra(const BAData& d) : data(d),
	regions(d.dims.numScenarios()) {

	// strange stuff happens with copy constructors if we have
	// vector<CoinBALPFactorization>
	int nscen = d.dims.numScenarios();
	f.resize(nscen);
	for (int i = 0; i < nscen; i++) {
		if (!data.ctx.assignedScenario(i)) { f[i] = 0; continue; }
		f[i] = new CoinBALPFactorization();
		f[i]->setCollectStatistics(true);
		f[i]->setNumScenariosPerProc(nscen/d.ctx.nprocs());
		regions[i].reserve(data.dims.numSecondStageCons(i));
	}

}

BALinearAlgebra::~BALinearAlgebra() {

	int nscen = data.dims.numScenarios();
	for (int i = 0; i < nscen; i++) {
		if (f[i]) delete f[i];
	}

}

void BALinearAlgebra::reinvert(const BAFlagVector<variableState> &v) {
#ifdef PIPSPROF
	MPI_Barrier(data.ctx.comm());
	double start_t = MPI_Wtime();
#endif

	const denseFlagVector<variableState> &v1 = v.getFirstStageVec();
	int nvar1 = data.dims.numFirstStageVars();
	int ncons1 = data.dims.numFirstStageCons();
	nbasic1 = 0;

	vector<int> colIsBasic1(nvar1);
	vector<int> basicCols1; // indices of basic first-stage variables
	basicCols1.reserve(nvar1);
	for (int i = 0; i < nvar1; i++) {
		if (v1[i] == Basic) {
			colIsBasic1[i] = 0;
			nbasic1++;
			basicCols1.push_back(i);
		}
		else colIsBasic1[i] = -1;

	}

	//printf("nbasic first stage: %d\n", nbasic1);

	int nscen = data.dims.numScenarios();
	int rowsFromThis = 0;	

	for (int i = 0; i < nscen; i++) {
		if (!data.ctx.assignedScenario(i)) continue;

		const denseFlagVector<variableState> &v2 = v.getSecondStageVec(i);
		int nvar2real = data.dims.inner.numSecondStageVars(i);
		int ncons2 = data.dims.numSecondStageCons(i);
		vector<int> colIsBasic2(nvar2real);
		vector<int> rowIsBasic2(ncons2);
		for (int j = 0; j < nvar2real; j++) {
			if (v2[j] == Basic) { colIsBasic2[j] = 0; }
			else colIsBasic2[j] = -1;

		}
		for (int j = 0; j < ncons2; j++) {
			if (v2[j+nvar2real] == Basic) { rowIsBasic2[j] = 0; }
			else rowIsBasic2[j] = -1;
		}
		int status = f[i]->factorizeRect(*data.Wcol[i],*data.Tcol[i],&colIsBasic2[0],&rowIsBasic2[0],&colIsBasic1[0],ncons1);
		assert(status == 0);
		f[i]->checkSparse();
		f[i]->goSparse();
		//f[i]->show_self();
		rowsFromThis += data.dims.numSecondStageCons(i) - f[i]->numberColumns();
	}
#ifdef PIPSPROF
	double local_t = MPI_Wtime() - start_t;
	MPI_Barrier(data.ctx.comm());
	tmp_t = MPI_Wtime();
#endif
	rowsPerProc.resize(data.ctx.nprocs());
	MPI_Allgather(&rowsFromThis,1,MPI_INT,&rowsPerProc[0],1,MPI_INT,data.ctx.comm());
	maxRowsIn = 0; rowsIn = 0;
	for (int i = 0; i < data.ctx.nprocs(); i++) {
		maxRowsIn = max(maxRowsIn, rowsPerProc[i]);
		if (i == data.ctx.mype()) myOffset = rowsIn;
		rowsIn += rowsPerProc[i];
	}
	assert(rowsIn + data.dims.numFirstStageCons() == nbasic1);
	assert(rowsFromThis == rowsPerProc[data.ctx.mype()]);
	ftranRecv.resize(maxRowsIn*data.ctx.nprocs());
	ftranSend.resize(maxRowsIn);
	btranSend.resize(nbasic1);

	reinvertFirstStage(basicCols1);
#ifdef PIPSPROF
	double totallocaltime;
	double maxlocaltime;
	MPI_Reduce(&local_t,&totallocaltime,1,MPI_DOUBLE,MPI_SUM,0,data.ctx.comm());
	MPI_Reduce(&local_t,&maxlocaltime,1,MPI_DOUBLE,MPI_MAX,0,data.ctx.comm());
	int nproc = data.ctx.nprocs();
	double idealtotal = totallocaltime/nproc + tmp2_t;
	if (data.ctx.mype() == 0) {
		printf("INVERT ideal total: %g actual total: %g load imbalance: %f comm. cost: %f 1st stage cost: %f\n",
			idealtotal,nproc*maxlocaltime+tmp_t+tmp2_t,
			100.*(maxlocaltime-totallocaltime/nproc)/idealtotal,
			100.*(tmp_t)/idealtotal,
			100.*(tmp2_t)/idealtotal);
	}
#endif

}


void BALinearAlgebra::reinvertFirstStage(const vector<int>& basicCols1) {
	if (nbasic1 == 0) return;
	if (nbasic1 > region1.capacity()) {
		region1.reserve(nbasic1);
	}

	double *elts = region1.denseVector();
	int *idx = region1.getIndices();
	
	int nscen = data.dims.numScenarios();
	int rowsSoFar = myOffset;
	for (int i = 0; i < nscen; i++) {
		if (!data.ctx.assignedScenario(i)) continue;
		int nbasic2 = f[i]->numberColumns();
		int nrows2 = data.dims.numSecondStageCons(i);
		int nRowsFromThis = nrows2 - nbasic2;
		assert(nRowsFromThis >= 0);
		for (int r = nbasic2; r < nrows2; r++) {
			f[i]->getRowOfT(region1,r);
			//printf("Z row %d:\n",rowsSoFar);
			//row.print();
			int nnz = region1.getNumElements();
			for (int q = 0; q < nnz; q++) {
				int col = idx[q];
				myElements.push_back(elts[col]);
				myIndicesColumn.push_back(col);
				myIndicesRow.push_back(rowsSoFar);
			}
			region1.clear();
			rowsSoFar++;
		}
	}
	rowsSoFar = rowsIn;
	
	int nMyElements = myElements.size();
	// gather rows from each MPI process
	vector<int> elementsCount(data.ctx.nprocs());
	MPI_Allgather(&nMyElements,1,MPI_INT,&elementsCount[0],1,MPI_INT,data.ctx.comm());
	CoinBigIndex nElements = 0;
	int maxsize = 0;
	for (int i = 0; i < data.ctx.nprocs(); i++) {
		nElements += elementsCount[i];
		maxsize = max(elementsCount[i],maxsize);
		//if (data.ctx.mype() == 0) printf("%d has %d elements\n",i,elementsCount[i]);
	}
	// pad so we don't read out of bounds (Allgather is faster than Allgatherv)
	myElements.reserve(maxsize); myIndicesColumn.reserve(maxsize); myIndicesRow.reserve(maxsize);
	double *recvElts = new double[maxsize*data.ctx.nprocs()];
	int *recvIdxRow = new int[maxsize*data.ctx.nprocs()];
	int *recvIdxCol = new int[maxsize*data.ctx.nprocs()];

	MPI_Allgather(&myElements[0],maxsize,MPI_DOUBLE,recvElts,maxsize,MPI_DOUBLE,data.ctx.comm());	
	MPI_Allgather(&myIndicesColumn[0],maxsize,MPI_INT,recvIdxCol,maxsize,MPI_INT,data.ctx.comm());	
	MPI_Allgather(&myIndicesRow[0],maxsize,MPI_INT,recvIdxRow,maxsize,MPI_INT,data.ctx.comm());

	myElements.clear(); myIndicesColumn.clear(); myIndicesRow.clear();
	allElements.reserve(nElements+1000); allElements.resize(nElements);
	allIndicesRow.reserve(nElements+1000); allIndicesRow.resize(nElements);
	allIndicesColumn.reserve(nElements+1000); allIndicesColumn.resize(nElements);

	CoinBigIndex nEltSoFar = 0;
	for (int i = 0; i < data.ctx.nprocs(); i++) {
		// note assumption of assignment of scenarios, that they're in increasing order wrt procs
		std::copy(recvElts+i*maxsize,recvElts+i*maxsize+elementsCount[i],&allElements[nEltSoFar]);
		std::copy(recvIdxRow+i*maxsize,recvIdxRow+i*maxsize+elementsCount[i],&allIndicesRow[nEltSoFar]);
		std::copy(recvIdxCol+i*maxsize,recvIdxCol+i*maxsize+elementsCount[i],&allIndicesColumn[nEltSoFar]);
		nEltSoFar += elementsCount[i];
	}
	
	delete [] recvElts;
	delete [] recvIdxRow;
	delete [] recvIdxCol;

	// now columns of A matrix
	const double *Aelts = data.Acol->getElements();
	const int *ArowIdx = data.Acol->getIndices();
	const CoinBigIndex *Astarts = data.Acol->getVectorStarts();
	int nvar1Real = data.dims.inner.numFirstStageVars();
	for (int i = 0; i < nbasic1; i++) {
		int col = basicCols1[i];
		if (col < nvar1Real) {
			for (CoinBigIndex j = Astarts[col]; j < Astarts[col+1]; j++) {
				int row = ArowIdx[j];
				allElements.push_back(Aelts[j]);
				allIndicesColumn.push_back(i);
				allIndicesRow.push_back(row+rowsSoFar);
			}
		} else { // slack
			int slackIdx = col - nvar1Real;
			allElements.push_back(-1);
			allIndicesColumn.push_back(i);
			allIndicesRow.push_back(slackIdx+rowsSoFar);
		}
	}
#ifdef PIPSPROF
	tmp2_t = MPI_Wtime();
	tmp_t = tmp2_t - tmp_t;
#endif	
	int status = f0.factorizeSquare(nbasic1,allElements.size(),&allIndicesRow[0],&allIndicesColumn[0],&allElements[0]);
#ifdef PIPSPROF
	tmp2_t = MPI_Wtime() - tmp2_t;
#endif
	assert(status == 0);
	allElements.clear(); allIndicesRow.clear(); allIndicesColumn.clear();
}

// a DUAL-sized vector is the right-hand side for FTRAN
// that is, has elements corresponding to each row
// but the solution is a PRIMAL-sized vector
// in regular simplex these are the same, but
// in our case, e.g. the number of first-stage duals isn't the same
// as the number of first-stage basic primals

// number of basic first-stage primals >= number of first-stage duals
// number of basic second-stage primals <= number of second-stage duals
// v needs space for both
void BALinearAlgebra::ftran(sparseBAVector &v) {

	const vector<int> &localScen = v.localScenarios();
	// see notes
	// step 1, FTRAN-G and pack Z_ir_i values
	int rowsSoFar = 0;
	for (unsigned j = 1; j < localScen.size(); j++) {
		int i = localScen[j];
		CoinIndexedVector &region = regions[i];
		CoinIndexedVector &v2 = v.getSecondStageVec(i).v;
		int nbasic2 = f[i]->numberColumns();
		int nrows2 = data.dims.numSecondStageCons(i);
		int nRowsFromThis = nrows2 - nbasic2;
		
		if (v2.getNumElements() != 0) { 
			f[i]->updateColumnG(&region,&v2);
			std::copy(region.denseVector()+nbasic2,region.denseVector()+nrows2,&ftranSend[rowsSoFar]);
		} else {
			std::fill(&ftranSend[rowsSoFar],&ftranSend[rowsSoFar+nRowsFromThis],0.0);
		}
		rowsSoFar += nRowsFromThis;
	}



	// collect Z_ir_i values
	MPI_Allgather(&ftranSend[0],maxRowsIn,MPI_DOUBLE,&ftranRecv[0],maxRowsIn,MPI_DOUBLE,data.ctx.comm());

	// unpack Z_ir_i values
	double *rhsElts = region1.denseVector();
	int *rhsIdx = region1.getIndices();
	int rhsNnz = 0;
	rowsSoFar = 0;
	for (int p = 0; p < data.ctx.nprocs(); p++) {
		int nRowsFromThis = rowsPerProc[p];
		for (int k = 0; k < nRowsFromThis; k++) {
			double val = ftranRecv[p*maxRowsIn+k];
			if (val) {
				int inRow = rowsSoFar+k;
				rhsElts[inRow] = val;
				rhsIdx[rhsNnz++] = inRow;
			}
		}
		rowsSoFar += nRowsFromThis;
	}

	int nrows1 = data.dims.numFirstStageCons();

	// first-stage rhs at the end
	CoinIndexedVector &v1 = v.getFirstStageVec().v;
	double *v1Elts = v1.denseVector();
	int *v1Idx = v1.getIndices();
	int v1Nnz = v1.getNumElements();
	for (int k = 0; k < v1Nnz; k++) {
		int row = v1Idx[k];
		double value = v1Elts[row];
		v1Elts[row] = 0.0;
		if (value == COIN_INDEXED_REALLY_TINY_ELEMENT) continue;
		assert(row < nrows1);
		int inRow = rowsSoFar+row;
		rhsElts[inRow] = value;
		rhsIdx[rhsNnz++] = inRow;
	}
	v1.setNumElements(0);
	region1.setNumElements(rhsNnz);
	
	// 1st stage solve

	if (nbasic1) f0.updateColumn(&v1,&region1);
	
	// swap would be better here
	v1 = region1;
	region1.clear();
	
	// steps 3 and 4

	for (unsigned j = 1; j < localScen.size(); j++) {
		int i = localScen[j];
		CoinIndexedVector &region = regions[i];
		int *regionIdx = region.getIndices();
		double *regionElts = region.denseVector();
		int nbasic2 = f[i]->numberColumns();
		int nnzRegion = region.getNumElements();
		
		// reuse region with Xr_i, but need to remove irrelevant indices
		int nnz = 0;
		for (int j = 0; j < nnzRegion; j++) {
			int row = regionIdx[j];
			if (row < nbasic2) {
				regionIdx[nnz++] = row;
			} else {
				regionElts[row] = 0.0;
			}
		}
		region.setNumElements(nnz);

		f[i]->multXT(v1,region);
		f[i]->updateColumnUQ(&region,&v.getSecondStageVec(i).v);

	}

}

// a PRIMAL-sized vector is the RHS for BTRAN
void BALinearAlgebra::btran(sparseBAVector &v) {

	const vector<int> &localScen = v.localScenarios();

	CoinIndexedVector &v1 = v.getFirstStageVec().v;
	int *v1Idx = v1.getIndices();
	double *v1Elts = v1.denseVector();
	double *sendbuf = &btranSend[0];
	double *r1Elts = region1.denseVector();
	std::fill(sendbuf, sendbuf+nbasic1, 0.0);

	// BTRAN-U

	for (unsigned j = 1; j < localScen.size(); j++) {
		int i = localScen[j];
		CoinIndexedVector &region = regions[i];
		CoinIndexedVector &v2 = v.getSecondStageVec(i).v;
		if (v2.getNumElements() == 0) continue;
		f[i]->updateColumnTransposeUQ(&region,&v2);
		if (nbasic1) f[i]->multXTTranspose(region, sendbuf);
	}
	
	// form first-stage rhs
	if (nbasic1) MPI_Allreduce(sendbuf, r1Elts, nbasic1, MPI_DOUBLE, MPI_SUM, data.ctx.comm());
	int v1Nnz = 0;
	for (int i = 0; i < nbasic1; i++) {
		double val = v1Elts[i] + r1Elts[i];
		if (fabs(val) > 1e-13) {
			v1Elts[i] = val;
			v1Idx[v1Nnz++] = i;
		} else {
			v1Elts[i] = 0.0;
		}
	}
	v1.setNumElements(v1Nnz);
	std::fill(r1Elts, r1Elts+nbasic1, 0.0);

	

	// step 4, first-stage solve

	
	double *bufvec = v1.denseVector(); // TODO: remove this
	if (nbasic1) f0.updateColumnTranspose(&region1,&v1);
	//region1.checkClear();
	
	// TODO: loop through nonzero elements of v1 instead?
	// step 5, copy out \beta vectors and do BTRAN-G
	int rowsSoFar = myOffset, nnz;
	for (unsigned j = 1; j < localScen.size(); j++) {
		int i = localScen[j];
		int nbasic2 = f[i]->numberColumns();
		int nrows2 = data.dims.numSecondStageCons(i);
		int nRowsFromThis = nrows2 - nbasic2;
		CoinIndexedVector &region = regions[i];
		CoinIndexedVector &v2 = v.getSecondStageVec(i).v;
		double *regionElts = region.denseVector();
		int *regionIdx = region.getIndices();
		nnz = region.getNumElements();
		for (int r = nbasic2; r < nrows2; r++) {
			double value = bufvec[r-nbasic2+rowsSoFar];
			if (value) {
				regionIdx[nnz++] = r;
				regionElts[r] = value;
			} else {
				regionElts[r] = 0.0;
			}
		}
		region.setNumElements(nnz);
		rowsSoFar += nRowsFromThis;
		f[i]->updateColumnTransposeG(&region,&v2);
		//printf("BTRAN from scen %d has %d nonzeros\n",i,v2.getNumElements());
	}

	int nrows1 = data.dims.numFirstStageCons();
	// clear out \beta vectors
	rowsSoFar = nbasic1 - nrows1;
	int n = v1.getNumElements();
	nnz = 0;
	//v1.print();
	memmove(v1Elts,v1Elts+rowsSoFar,nrows1*sizeof(double));
	memset(v1Elts+nrows1,0,rowsSoFar*sizeof(double));
	for (int j = 0; j < n; j++) {
		int i = v1Idx[j];
		int diff = i-rowsSoFar;
		if (diff >= 0) {
			v1Idx[nnz++] = diff;
		//	printf("%d %d %f %f\n",nnz-1,diff,value, v1Elts[diff]);
		}
	}
	v1.setNumElements(nnz);
	//printf("BTRAN from scen -1 has %d nonzeros\n",v1.getNumElements());
	//printf("NEW: %d\n", rowsSoFar);
	//v1.print();
	//v1.checkClean();



}
