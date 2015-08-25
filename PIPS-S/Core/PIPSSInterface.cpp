#include "PIPSSInterface.hpp"
#include "BALPSolverDual.hpp"
#include "BALPSolverPrimal.hpp"
#include "CoinPackedVector.hpp"

using namespace std;

PIPSSInterface::PIPSSInterface(stochasticInput &in, BAContext &ctx, solveType t) : d(in,ctx), boundsChanged(false), st(t) {

	if (t == usePrimal) {
		solver = new BALPSolverPrimal(d);
	} else {
		solver = new BALPSolverDual(d);
	}
	if (d.ctx.mype() == 0) {
	  PIPS_APP_LOG_SEV(summary)<<boost::format("First stage: %d cons %d vars")
	    % d.dims.numFirstStageCons() %  d.dims.inner.numFirstStageVars();
	  PIPS_APP_LOG_SEV(summary)<<boost::format("Second stage: %d cons %d vars %d scenarios")
	    % d.dims.numSecondStageCons(0) % d.dims.inner.numSecondStageVars(0) % d.dims.numScenarios();
	  PIPS_APP_LOG_SEV(info)<<boost::format("reinvert every %d") % solver->reinvertFrequency;
	}



}

PIPSSInterface::PIPSSInterface(const BAData& _d, solveType t) : d(_d), boundsChanged(false), st(t) {

	if (t == usePrimal) {
		solver = new BALPSolverPrimal(d);
	} else {
		solver = new BALPSolverDual(d);
	}
}


PIPSSInterface::~PIPSSInterface() {
	delete solver;
}


void PIPSSInterface::go() {
	double t = MPI_Wtime();
	int mype = d.ctx.mype();

	// Reallocate if bounds changed for primal solve, change to dual solve
	if (boundsChanged && st == usePrimal) {
		BALPSolverBase *solver2 = new BALPSolverDual(d);
		solver2->setStates(solver->getStates());
		solver2->setPrimalTolerance(solver->getPrimalTolerance());
		solver2->setDualTolerance(solver->getDualTolerance());
		solver2->phase1 = solver->phase1;
		delete solver;
		solver = solver2;
		boundsChanged = false;
		st = useDual;
	}

	solver->go();

	// clean up infeasibilities
	int count = 0;
	while (solver->getStatus() != Optimal && count < 10) {
		if (solver->getStatus() == ProvenInfeasible ||
		    solver->getStatus() == ProvenUnbounded) break;
		BALPSolverBase *solver2;
		if (solver->getStatus() == DualFeasible) {
			solver2 = new BALPSolverDual(d);
			if (mype == 0) PIPS_APP_LOG_SEV(info)<<"Switching to dual";
			st = useDual;
		} else {
			assert(solver->getStatus() == PrimalFeasible);
			solver2 = new BALPSolverPrimal(d);
			if (mype == 0) PIPS_APP_LOG_SEV(info)<<"Switching to primal";
			st = usePrimal;
		}
		solver2->setStates(solver->getStates());
		solver2->setPrimalTolerance(solver->getPrimalTolerance());
		solver2->setDualTolerance(solver->getDualTolerance());
		solver2->setReinversionFrequency(5); // be very careful

		solver2->ftranTime = solver->ftranTime;
		solver2->btranTime = solver->btranTime;
		solver2->updateIteratesTime = solver->updateIteratesTime;
		solver2->ftranDSETime = solver->ftranDSETime;
		solver2->invertTime = solver->invertTime;
		solver2->selectEnteringTime = solver->selectEnteringTime;
		solver2->selectLeavingTime = solver->selectLeavingTime;
		solver2->priceTime = solver->priceTime;
		solver2->updateColumnTime = solver->updateColumnTime;
		solver2->startTime = solver->startTime;

		solver2->replaceFirst = solver->replaceFirst;
		solver2->firstReplaceSecond = solver->firstReplaceSecond;
		solver2->secondReplaceFirst = solver->secondReplaceFirst;
		solver2->replaceSecondSelf = solver->replaceSecondSelf;
		solver2->replaceSecondOther = solver->replaceSecondOther;

		solver2->nIter = solver->nIter;
		solver2->phase1 = solver->phase1;

		delete solver; // free memory
		solver2->go();
		solver = solver2;
	}

 	if (solver->getStatus() != Optimal) {
	  //if (mype == 0)
	   // PIPS_APP_LOG_SEV(fatal)<<boost::format("Switched between primal and dual %d times, and still not optimal!")
	     // % count;
	 // MPI_Abort(MPI_COMM_WORLD,1);
	}
	t = MPI_Wtime() - t;

	if (mype == 0)
	  PIPS_APP_LOG_SEV(summary)<<boost::format("Solve took %f seconds") % t;


}

std::vector<double> PIPSSInterface::getFirstStagePrimalColSolution() const {
	const denseVector &x = solver->getPrimalSolution().getFirstStageVec();
	int nvar1real = d.dims.inner.numFirstStageVars();
	return std::vector<double>(&x[0],&x[nvar1real]);
}

std::vector<double> PIPSSInterface::getSecondStagePrimalColSolution(int scen) const {
	assert(d.ctx.assignedScenario(scen));
	const denseVector &x = solver->getPrimalSolution().getSecondStageVec(scen);
	int nvar2real = d.dims.inner.numSecondStageVars(scen);
	return std::vector<double>(&x[0],&x[nvar2real]);
}

std::vector<double> PIPSSInterface::getFirstStageDualColSolution() const {
	const denseVector &x = solver->getDualColSolution().getFirstStageVec();
	int nvar1real = d.dims.inner.numFirstStageVars();
	return std::vector<double>(&x[0],&x[nvar1real]);
}

std::vector<double> PIPSSInterface::getSecondStageDualColSolution(int scen) const {
	assert(d.ctx.assignedScenario(scen));
	const denseVector &x = solver->getDualColSolution().getSecondStageVec(scen);
	int nvar2real = d.dims.inner.numSecondStageVars(scen);
	return std::vector<double>(&x[0],&x[nvar2real]);
}

std::vector<double> PIPSSInterface::getSecondStageDualRowSolution(int scen) const {
	assert(d.ctx.assignedScenario(scen));
	const sparseVector &x = solver->btranVec.getSecondStageVec(scen);
	int ncons2 = d.dims.inner.numSecondStageCons(scen);
	return std::vector<double>(&x[0],&x[ncons2]);
}

void PIPSSInterface::setFirstStageColState(int idx,variableState s) {
	solver->states.getFirstStageVec()[idx] = s;
}

void PIPSSInterface::setFirstStageRowState(int idx,variableState s) {
	int nvarreal = d.dims.inner.numFirstStageVars();
	solver->states.getFirstStageVec()[idx+nvarreal] = s;
}



void PIPSSInterface::setSecondStageColState(int scen, int idx,variableState s) {
	solver->states.getSecondStageVec(scen)[idx] = s;
}

void PIPSSInterface::setSecondStageRowState(int scen, int idx,variableState s) {
	int nvarreal = d.dims.inner.numSecondStageVars(scen);
	solver->states.getSecondStageVec(scen)[idx+nvarreal] = s;
}
void PIPSSInterface::commitStates() {
	solver->status = LoadedFromFile;
	solver->setupIndices();
}

variableState PIPSSInterface::getFirstStageColState(int idx) const {
	return solver->states.getFirstStageVec()[idx];
}

variableState PIPSSInterface::getFirstStageRowState(int idx) const {
	int nvarreal = d.dims.inner.numFirstStageVars();
	return solver->states.getFirstStageVec()[nvarreal+idx];
}

variableState PIPSSInterface::getSecondStageColState(int scen, int idx) const {
	return solver->states.getSecondStageVec(scen)[idx];
}

variableState PIPSSInterface::getSecondStageRowState(int scen, int idx) const {
	int nvarreal = d.dims.inner.numSecondStageVars(scen);
	return solver->states.getSecondStageVec(scen)[nvarreal+idx];
}


void PIPSSInterface::setFirstStageColLB(int idx, double newLb) {
  d.l.getFirstStageVec()[idx] = newLb;
  boundsChanged = true;
}

void PIPSSInterface::setFirstStageColUB(int idx, double newUb) {
  d.u.getFirstStageVec()[idx] = newUb;
  boundsChanged = true;
}

void PIPSSInterface::setSecondStageColLB(int scen, int idx, double newLb) {
        d.l.getSecondStageVec(scen)[idx] = newLb;
        boundsChanged = true;
}

void PIPSSInterface::setSecondStageColUB(int scen, int idx, double newUb) {
        d.u.getSecondStageVec(scen)[idx] = newUb;
        boundsChanged = true;
}

void PIPSSInterface::setLB(const denseBAVector& lb) {
        // use denseBAVector::copyFrom?
        d.l.copyFrom(lb);
	boundsChanged = true;
}

void PIPSSInterface::setUB(const denseBAVector& ub) {
        d.u.copyFrom(ub);
        boundsChanged = true;
}

void PIPSSInterface::addRow(const std::vector<double>& elts1, const std::vector<double> &elts2, int scen, double lb, double ub) {

	CoinPackedVector e1;
	CoinPackedVector e2;
	if (d.ctx.assignedScenario(scen)) {
		e1.setFullNonZero(elts1.size(),&elts1[0]);
		e2.setFullNonZero(elts2.size(),&elts2[0]);
	}
	d.addRow(e1,e2,scen,lb,ub);


}

void PIPSSInterface::commitNewRows() {

	const vector<int> &localScen = d.ctx.localScenarios();

	assert(solver->status == Optimal);

	BALPSolverDual* solver2 = new BALPSolverDual(d);
	solver2->setPrimalTolerance(solver->getPrimalTolerance());
	solver2->setDualTolerance(solver->getDualTolerance());

	solver2->states.getFirstStageVec().copyFrom(solver->states.getFirstStageVec());
	for (unsigned i = 1; i < localScen.size(); i++) {
		int scen = localScen[i];

		denseFlagVector<variableState> &oldStates = solver->states.getSecondStageVec(scen), &newStates = solver2->states.getSecondStageVec(scen);

		copy(&oldStates[0],&oldStates[oldStates.length()],&newStates[0]);
		printf("%d new rows in scenario %d\n",newStates.length()-oldStates.length(),scen);
		for (int k = oldStates.length(); k < newStates.length(); k++) {
			newStates[k] = Basic;
		}
	}

	delete solver;
	solver = solver2;
	st = useDual;

	commitStates();


}


	const CoinShallowPackedVector PIPSSInterface::retrieveARow(int index ) const {
		return d.retrieveARow(index);
	}

	const CoinShallowPackedVector PIPSSInterface::retrieveWRow(int index,int scen) const{
		return d.retrieveWRow(index,scen);
	}

	const CoinShallowPackedVector PIPSSInterface::retrieveTRow(int index,int scen) const{
		return d.retrieveTRow(index,scen);
	}

	const CoinShallowPackedVector PIPSSInterface::retrieveACol(int index) const{ 
		return d.retrieveACol(index);
	}

	const CoinShallowPackedVector PIPSSInterface::retrieveWCol(int index,int scen) const{
		return d.retrieveWCol(index,scen);
	}

	const CoinShallowPackedVector PIPSSInterface::retrieveTCol (int index,int scen) const{
		return d.retrieveTCol(index,scen);
	}

	//returns 1 if feasible, returns -1 if upper bound broken, returns -2 if lower bound broken
	int PIPSSInterface::isRowFeasible(int index, int scen, denseBAVector& solution){
		//Calculate expression,
		//TODO: Change the hardcoded tolerances for a centralized control variable
		double intTol=10e-9;
		double expression=0;
		if (scen<0){
			CoinShallowPackedVector row=d.retrieveARow(index);
			int nElems=row.getNumElements();
			const int *indices=row.getIndices();
			const double *elems=row.getElements(); 

	    	for (int el=0; el<nElems; el++){
	    		expression+=(elems[el]*solution.getFirstStageVec()[indices[el]]);

	    	}
	    	double lb=d.l.getFirstStageVec()[d.dims.inner.numFirstStageVars()+index];
	    	double ub=d.u.getFirstStageVec()[d.dims.inner.numFirstStageVars()+index];
	    	//cout<< std::setprecision(15)<<" In the end our row "<<index<<" had "<<lb<<" "<<expression<<" "<<ub<<endl;
	    	if (expression<lb-intTol) return -2;
	    	if (expression>ub+intTol) return -1;
	    	return 1;

		}
		else{
			assert(scen >= 0 && scen < d.dims.numScenarios());
			if (d.ctx.assignedScenario(scen)){
				CoinShallowPackedVector row=d.retrieveTRow(index,scen);
				int nElems=row.getNumElements();
				const int *indices=row.getIndices();
				const double *elems=row.getElements(); 

		    	for (int el=0; el<nElems; el++){
		    		//cout<<expression<<"+="<<elems[el]*solution.getFirstStageVec()[indices[el]]<<" "<<elems[el]<<" "<<solution.getFirstStageVec()[indices[el]]<<" "<<indices[el]<<endl;
		    		expression+=(elems[el]*solution.getFirstStageVec()[indices[el]]);

		    	}

		    	CoinShallowPackedVector row2=d.retrieveWRow(index,scen);
				int nElems2=row2.getNumElements();
				const int *indices2=row2.getIndices();
				const double *elems2=row2.getElements(); 

		    	for (int el=0; el<nElems2; el++){
		    		//cout<<expression<<"+="<<elems2[el]*solution.getSecondStageVec(scen)[indices2[el]]<<" "<<elems2[el]<<" "<<solution.getSecondStageVec(scen)[indices2[el]]<<" s"<<scen<<" "<<indices2[el]<<endl;
		    		expression+=(elems2[el]*solution.getSecondStageVec(scen)[indices2[el]]);

		    	}
		    	double lb=d.l.getSecondStageVec(scen)[d.dims.inner.numSecondStageVars(scen)+index];
		    	double ub=d.u.getSecondStageVec(scen)[d.dims.inner.numSecondStageVars(scen)+index];
		    	//cout<<std::setprecision(15)<<" In the end our row "<<index<<" "<<scen<<" had "<<lb<<" "<<expression<<" "<<ub<<endl;
	    	
		    	if (expression<lb-intTol) return -2;
		    	if (expression>ub+intTol) return -1;
		    	return 1;

		    }

		    else return 1;

		}


		return 0;
	}


	void PIPSSInterface::commitNewColsAndRows()  {
	
		assert(solver->status == Optimal);
	
		BALPSolverDual* solver2 = new BALPSolverDual(d);
		solver2->setPrimalTolerance(solver->getPrimalTolerance());
		solver2->setDualTolerance(solver->getDualTolerance());
		delete solver;
		solver = solver2;
		st = useDual;
		solver->status = Uninitialized;
		
	}
	
	void PIPSSInterface::generateBetas(sparseBAVector &beta){
		solver->generateBetas(beta);
	}

	void PIPSSInterface::generateNonBasicRow(BAIndex in, sparseBAVector &row){
		solver->generateNonBasicRow(in, row);
	}

	void PIPSSInterface::addFirstStageRow(const std::vector<double>& elts1, double lb , double ub ){
		CoinPackedVector e1;
		e1.setFullNonZero(elts1.size(),&elts1[0]);
		d.addFirstStageRow(e1,lb,ub);
		cout<<"Added first stage row "<<d.dims.inner.numFirstStageVars()<<" and cons "<<d.dims.numFirstStageCons()<<endl;
	}
		
	void PIPSSInterface::addSecondStageRows(const std::vector< std::vector <double> >& elts1, const std::vector< std::vector <double> >&elts2, int scen, std::vector<double> &lb, std::vector<double> &ub, int nRows) {
	 	std::vector<CoinPackedVector*> vectors1(nRows);
		std::vector<CoinPackedVector*> vectors2(nRows);
	 	for (int i=0; i< nRows; i++){
	 		vectors1[i]= new CoinPackedVector();
	 		vectors2[i]= new CoinPackedVector();
	 		vectors1[i]->setFullNonZero(elts1[i].size(),&elts1[i][0]);
			vectors2[i]->setFullNonZero(elts2[i].size(),&elts2[i][0]);
	 	}
		d.addSecondStageConsecutiveRows(vectors1, vectors2, scen, lb, ub, nRows);
		for (int i=0; i< nRows; i++){
			delete vectors1[i];
			delete vectors2[i];
		}
	}

	void PIPSSInterface::addFirstStageRows(const std::vector< std::vector <double> >& elts, std::vector<double> &lb, std::vector<double> &ub, int nRows) {
	 	std::vector<CoinPackedVector*> vectors(nRows);
	 	for (int i=0; i< nRows; i++){
	 		vectors[i]= new CoinPackedVector();
	 		vectors[i]->setFullNonZero(elts[i].size(),&elts[i][0]);	
	 	}
		d.addFirstStageRows(vectors, lb, ub, nRows);
		for (int i=0; i< nRows; i++){
			delete vectors[i];
		}
	}

	void PIPSSInterface::deleteLastFirstStageRows(int nRows){
		d.deleteLastFirstStageRows(nRows);
	}
		
	void PIPSSInterface::deleteLastFirstStageColumns(int nCols){
		d.deleteLastFirstStageColumns(nCols);
	}

	void PIPSSInterface::deleteLastSecondStageConsecutiveRows(int scenario, int nRows){
		d.deleteLastSecondStageConsecutiveRows(scenario,nRows);
	}

	int PIPSSInterface::addFirstStageColumn(double lb, double ub, double c){
		int out= d.addFirstStageColumn(lb,ub,c);
		return out;
	}

	int PIPSSInterface::addSecondStageColumn(int scen,double lb, double ub, double cobj){
		int out= d.addSecondStageColumn(scen,lb,ub,cobj);
		return out;
	}
