#include "proximalBAQP.hpp"


using namespace std;


proximalQPModel::proximalQPModel(int nvar1, bundle_t const &cuts, vector<vector<double> > const &proxCenter, double tau) :
	nvar1(nvar1), cuts(cuts), proxCenter(proxCenter), tau(tau) {}

vector<double> proximalQPModel::getFirstStageColLB() {
	// y free
	return vector<double>(nvar1,-COIN_DBL_MAX);
}

vector<double> proximalQPModel::getFirstStageColUB() {
	return vector<double>(nvar1,COIN_DBL_MAX);
}

vector<double> proximalQPModel::getFirstStageObj() {
	vector<double> obj(nvar1+1,0.);
	return vector<double>(nvar1,0.);
}
vector<string> proximalQPModel::getFirstStageColNames() {
	vector<string> names(nvar1, "y");
	return names;
}

vector<double> proximalQPModel::getFirstStageRowLB() {
	return vector<double>(0);
}
vector<double> proximalQPModel::getFirstStageRowUB() {
	return vector<double>(0);
}

vector<string> proximalQPModel::getFirstStageRowNames() {
	return vector<string>(0);
}


vector<double> proximalQPModel::getSecondStageColLB(int scen) {
	return vector<double>(cuts[scen].size(),0.);
}

vector<double> proximalQPModel::getSecondStageColUB(int scen) {
	return vector<double>(cuts[scen].size(),COIN_DBL_MAX);
}

// f_i(\gamma_i^r) + g^T(\gamma_i^+-\gamma_i^r). also, \hat f_i^r
double cutInfo::computeC(vector<double> const &proxCenter) const {
	assert(evaluatedAt.size() == subgradient.size() && subgradient.size() == proxCenter.size());
	double out = objval;
	for (unsigned i = 0; i < subgradient.size(); i++) {
		out += subgradient[i]*(proxCenter[i]-evaluatedAt[i]);
	}
	return out;
}

vector<double> proximalQPModel::getSecondStageObj(int scen) {
	vector<double> obj;
	obj.reserve(cuts[scen].size());
	for (unsigned i = 0; i < cuts[scen].size(); i++) {
		obj.push_back(-cuts[scen][i].computeC(proxCenter[scen]));
	}
	return obj;
}

vector<string> proximalQPModel::getSecondStageColNames(int scen) {
	vector<string> names(cuts[scen].size(),"u");
	return names;
}

vector<double> proximalQPModel::getSecondStageRowUB(int scen) {
	vector<double> ub(1,1.);
	return ub;
}

vector<double> proximalQPModel::getSecondStageRowLB(int scen) {
	return getSecondStageRowUB(scen);
}

vector<string> proximalQPModel::getSecondStageRowNames(int scen) {
	vector<string> names(1,"theta");
	return names;

}

CoinPackedMatrix proximalQPModel::getFirstStageConstraints() {
	vector<CoinBigIndex> starts(nvar1+1,0);
	return CoinPackedMatrix(true, 0, nvar1,0,0,0,&starts[0],0);
}

CoinPackedMatrix proximalQPModel::getSecondStageConstraints(int scen) {
	/* diagonal blocks are of the form:
	   [  e^T ]
	*/	
	unsigned ncol = cuts[scen].size();
	vector<double> elts(ncol,1.);
	vector<int> idx(ncol,0);
	vector<CoinBigIndex> starts(ncol+1);
	for (unsigned col = 0; col <= ncol; col++) starts[col] = col;


	return CoinPackedMatrix(true,1,ncol,ncol,&elts[0],&idx[0],&starts[0],0);

}

CoinPackedMatrix proximalQPModel::getLinkingConstraints(int scen) {
	// this is easy, empty!

	vector<CoinBigIndex> starts(nvar1+1,0);
	return CoinPackedMatrix(true, 1, nvar1,0,0,0,&starts[0],0);

}

CoinPackedMatrix proximalQPModel::getFirstStageHessian() {
	// [ (1/\tau) I ]
	vector<double> elts(nvar1,1./tau);
	vector<int> idx(nvar1);
	vector<CoinBigIndex> starts(nvar1+1);

	for (int i = 0; i < nvar1; i++) {
		idx[i] = starts[i] = i;
	}
	starts[nvar1] = nvar1;

	return CoinPackedMatrix(true,nvar1,nvar1,nvar1,&elts[0],&idx[0],&starts[0],0);
}

CoinPackedMatrix proximalQPModel::getSecondStageHessian(int scen) {
	// [ (1/\tau) G_iG_i^T ]
	CoinBigIndex nnz = 0;
	int ncol = cuts[scen].size();
	vector<double> elts; elts.reserve(ncol*(ncol+1)/2);
	vector<int> idx; idx.reserve(ncol*(ncol+1)/2);
	vector<CoinBigIndex> starts(ncol+1);

	for (int j = 0; j < ncol; j++) {
		starts[j] = nnz;
		vector<double> const& vj = cuts[scen][j].subgradient;
		for (int i = j; i < ncol; i++) {
			// dot product between ith and jth subgradient
			vector<double> const& vi = cuts[scen][i].subgradient;
			double sum = 0;
			for (int k = 0; k < nvar1; k++) {
				sum += vi[k]*vj[k];
			}
			if (fabs(sum) > 1e-12) {
				sum /= tau;
				elts.push_back(sum);
				idx.push_back(i);
				nnz++;
			}
		}
	}
	starts[ncol] = nnz;

	return CoinPackedMatrix(true,ncol,ncol,nnz,&elts[0],&idx[0],&starts[0],0);
}

CoinPackedMatrix proximalQPModel::getSecondStageCrossHessian(int scen) {
	// [ - 1/(tau \sqrt{N}) G_i ]

	CoinBigIndex nnz = 0;
	int ncol = nvar1;
	int nrow = cuts[scen].size();
	vector<double> elts; elts.reserve(nrow*ncol);
	vector<int> idx; idx.reserve(nrow*ncol);
	vector<CoinBigIndex> starts(ncol+1);

	for (int j = 0; j < ncol; j++) {
		starts[j] = nnz;
		for (int i = 0; i < nrow; i++) {
			double elt = cuts[scen][i].subgradient[j];
			if (fabs(elt) > 1e-12) {
				elt /= -tau*sqrt((double)cuts.size());
				elts.push_back(elt);
				idx.push_back(i);
				nnz++;
			}
		}
	}
	starts[ncol] = nnz;

	return CoinPackedMatrix(true,nrow,ncol,nnz,&elts[0],&idx[0],&starts[0],0);
}
