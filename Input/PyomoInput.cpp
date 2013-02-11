#include "PyomoInput.hpp"

using namespace std;

PyomoInput::PyomoInput(const string &root_filename, int num_scens) 
{	
  nScenarios_ = num_scens;
  nl_root_filename = root_filename;

  nSecondStageVars_.resize(nScenarios_);
  nSecondStageCons_.resize(nScenarios_);
	
	nFirstStageVars_ = em->getNLocalVars();
	nFirstStageCons_ = em->getNLocalCons();

	for (int i = 0; i < nScenarios_; i++) {
		ExpandedModelInterface *em2 = em->children[i];
		assert(em2->children.size() == 0); // only support one level
		nSecondStageVars_[i] = em2->getNLocalVars();
		nSecondStageCons_[i] = em2->getNLocalCons();
	}

	dimsEqual = true;
	for (int i = 0; i < nScenarios_; i++) {
		if (nSecondStageVars_[i] != nSecondStageVars_[0] ||
		    nSecondStageCons_[i] != nSecondStageCons_[0]) {
		    dimsEqual = false; break;
		}
	}

	// make sure we have a dual block-angular program
	for (int i = 0; i < nScenarios_; i++) {
		ExpandedModelInterface *em2 = em->children[i];
		assert(em->getNzJacobianOfIntersection(em2) == 0);
	}
	
}

PyomoInput::~PyomoInput()
{

}

vector<double> PyomoInput::getFirstStageColLB() {
	vector<double> v(nFirstStageVars_);
	em->getColLowBounds(&v[0]);
	return v;
}

vector<double> PyomoInput::getFirstStageColUB() {
	vector<double> v(nFirstStageVars_);
	em->getColUpBounds(&v[0]);
	return v;
}

vector<double> PyomoInput::getFirstStageObj() {
	vector<double> v(nFirstStageVars_);
	em->getObjGradient(&v[0]);
	return v;
}

vector<string> PyomoInput::getFirstStageColNames() {
	const list<string> &col = em->getLocalVarNames();
	return vector<string>(col.begin(),col.end());
}

vector<double> PyomoInput::getFirstStageRowLB() {
	vector<double> v1(nFirstStageCons_), v2(nFirstStageCons_);
	em->getRowBounds(&v1[0],&v2[0]);
	return v1;
}

vector<double> PyomoInput::getFirstStageRowUB() {
	vector<double> v1(nFirstStageCons_), v2(nFirstStageCons_);
	em->getRowBounds(&v1[0],&v2[0]);
	return v2;
}

vector<string> PyomoInput::getFirstStageRowNames() {
	const list<string> &row = em->getLocalConNames();
	return vector<string>(row.begin(),row.end());
}

vector<double> PyomoInput::getSecondStageColLB(int s) {
	vector<double> v(nSecondStageVars_[s]);
	em->children[s]->getColLowBounds(&v[0]);
	return v;
}

vector<double> PyomoInput::getSecondStageColUB(int s) {
	vector<double> v(nSecondStageVars_[s]);
	em->children[s]->getColUpBounds(&v[0]);
	return v;
}

vector<double> PyomoInput::getSecondStageObj(int s) {
	vector<double> v(nSecondStageVars_[s]);
	em->children[s]->getObjGradient(&v[0]);
	return v;
}

vector<string> PyomoInput::getSecondStageColNames(int s) {
	const list<string> &col = em->children[s]->getLocalVarNames();
	return vector<string>(col.begin(),col.end());
}

vector<double> PyomoInput::getSecondStageRowLB(int s) {
	vector<double> v1(nSecondStageCons_[s]), v2(nSecondStageCons_[s]);
	em->children[s]->getRowBounds(&v1[0],&v2[0]);
	return v1;
}

vector<double> PyomoInput::getSecondStageRowUB(int s) {
	vector<double> v1(nSecondStageCons_[s]), v2(nSecondStageCons_[s]);
	em->children[s]->getRowBounds(&v1[0],&v2[0]);
	return v2;
}

vector<string> PyomoInput::getSecondStageRowNames(int s) {
	const list<string> &row = em->children[s]->getLocalConNames();
	return vector<string>(row.begin(),row.end());
}

CoinPackedMatrix PyomoInput::getFirstStageConstraints() {

	int nnz = em->getNzJacobianOfIntersection(em);
	vector<int> colbeg(nFirstStageVars_+1);
	vector<int> collen(nFirstStageVars_);
	vector<int> rowidx(nnz);
	vector<double> elts(nnz);

	em->getJacobianOfIntersection(em, &colbeg[0], &collen[0], &rowidx[0], &elts[0]);

	vector<CoinBigIndex> colbeg2(colbeg.begin(),colbeg.end());

	return CoinPackedMatrix(true,nFirstStageCons_,nFirstStageVars_,
		nnz, &elts[0], &rowidx[0], &colbeg2[0], &collen[0]);

}

CoinPackedMatrix PyomoInput::getLinkingConstraints(int idx) {

	int nnz = em->children[idx]->getNzJacobianOfIntersection(em);
	vector<int> colbeg(nFirstStageVars_+1);
	vector<int> collen(nFirstStageVars_);
	vector<int> rowidx(nnz);
	vector<double> elts(nnz);

	em->children[idx]->getJacobianOfIntersection(em, &colbeg[0], &collen[0], &rowidx[0], &elts[0]);

	vector<CoinBigIndex> colbeg2(colbeg.begin(),colbeg.end());

	return CoinPackedMatrix(true,nSecondStageCons_[idx],nFirstStageVars_,
		nnz, &elts[0], &rowidx[0], &colbeg2[0], &collen[0]);

}


CoinPackedMatrix PyomoInput::getSecondStageConstraints(int idx) {

	int nnz = em->children[idx]->getNzJacobianOfIntersection(em->children[idx]);
	vector<int> colbeg(nSecondStageVars_[idx]+1);
	vector<int> collen(nSecondStageVars_[idx]);
	vector<int> rowidx(nnz);
	vector<double> elts(nnz);

	em->children[idx]->getJacobianOfIntersection(em->children[idx], &colbeg[0], &collen[0], &rowidx[0], &elts[0]);

	vector<CoinBigIndex> colbeg2(colbeg.begin(),colbeg.end());

	return CoinPackedMatrix(true,nSecondStageCons_[idx],nSecondStageVars_[idx],
		nnz, &elts[0], &rowidx[0], &colbeg2[0], &collen[0]);

}



