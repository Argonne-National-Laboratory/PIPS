#include "l1TrustBALP.hpp"
#include <cmath>


using namespace std;


l1TrustModel::l1TrustModel(int nvar1, std::vector<std::vector<cutInfo> > const &cuts, double trustRadius, std::vector<std::vector<double> > const& center) :
	nvar1(nvar1), cuts(cuts), center(center), trustRadius(trustRadius) {}

vector<double> l1TrustModel::getFirstStageColLB() {
	vector<double> l(nvar1+1,-COIN_DBL_MAX);
	l[nvar1] = 0.;

	return l;
}

vector<double> l1TrustModel::getFirstStageColUB() {

	return vector<double>(nvar1+1,COIN_DBL_MAX);
}

vector<double> l1TrustModel::getFirstStageObj() {

	vector<double> obj(nvar1+1,0.);
	obj[nvar1] = trustRadius;
	return obj;
}
vector<string> l1TrustModel::getFirstStageColNames() {

	vector<string> names(nvar1+1);
	return names;
}

vector<double> l1TrustModel::getFirstStageRowLB() {
	return vector<double>(0);
}
vector<double> l1TrustModel::getFirstStageRowUB() {
	return vector<double>(0);
}

vector<string> l1TrustModel::getFirstStageRowNames() {
	return vector<string>(0);
}


vector<double> l1TrustModel::getSecondStageColLB(int scen) {
	return vector<double>(nSecondStageVars(scen),0.);
}

vector<double> l1TrustModel::getSecondStageColUB(int scen) {
	return vector<double>(nSecondStageVars(scen),COIN_DBL_MAX);
}


vector<double> l1TrustModel::getSecondStageObj(int scen) {
	vector<double> obj;
	obj.reserve(nSecondStageVars(scen));
	for (int i = 0; i < nvar1; i++) {
		obj.push_back(-center[scen][i]);
	}
	for (int i = 0; i < nvar1; i++) {
		obj.push_back(center[scen][i]);
	}
	for (unsigned i = 0; i < cuts[scen].size(); i++) {
		obj.push_back(-cuts[scen][i].computeC());
	}
	return obj;
}

vector<string> l1TrustModel::getSecondStageColNames(int scen) {
	vector<string> names(nSecondStageVars(scen));
	return names;
}

vector<double> l1TrustModel::getSecondStageRowUB(int scen) {
	vector<double> ub(nSecondStageCons(scen),0.);
	ub[0] = 1.;
	return ub;
}

vector<double> l1TrustModel::getSecondStageRowLB(int scen) {
	return getSecondStageRowUB(scen);
}

vector<string> l1TrustModel::getSecondStageRowNames(int scen) {
	vector<string> names(nSecondStageCons(scen));
	return names;

}

CoinPackedMatrix l1TrustModel::getFirstStageConstraints() {
	vector<CoinBigIndex> starts(nvar1+2,0);
	return CoinPackedMatrix(true, 0, nvar1+1,0,0,0,&starts[0],0);
}

CoinPackedMatrix l1TrustModel::getSecondStageConstraints(int scen) {
	/* diagonal blocks are of the form:
	   [       e^T   ]
	   [ I -I -G_i^T ]
	   [ I  I        ]
	*/
	
	CoinBigIndex nnz = 0;
	int ncuts = cuts[scen].size();
	int ncol = ncuts + 2*nvar1;
	vector<double> elts; elts.reserve(2*nvar1+ncuts*(nvar1+1));
	vector<int> idx; idx.reserve(2*nvar1+ncuts*(nvar1+1));
	vector<CoinBigIndex> starts(ncol+1);
	for (int i = 0; i < nvar1; i++) {
		starts[i] = nnz;
		elts.push_back(1.);
		elts.push_back(1.);
		idx.push_back(1+i);
		idx.push_back(1+nvar1+i);
		nnz += 2;
	}
	for (int i = 0; i < nvar1; i++) {
		starts[i+nvar1] = nnz;
		elts.push_back(-1.);
		elts.push_back(1.);
		idx.push_back(1+i);
		idx.push_back(1+nvar1+i);
		nnz += 2;
	}

	for (int i = 0; i < ncuts; i++) {
		starts[i+2*nvar1] = nnz++;
		elts.push_back(1.);
		idx.push_back(0);
		vector<double> const& subgrad = cuts[scen][i].subgradient;
		for (unsigned k = 0; k < subgrad.size(); k++) {
			if (fabs(subgrad[k]) > 1e-7) {
				elts.push_back(-subgrad[k]);
				idx.push_back(k+1);
				nnz++;
			}
		}
	}
	starts[ncol] = nnz;

	return CoinPackedMatrix(true,2*nvar1+1,ncol,nnz,&elts[0],&idx[0],&starts[0],0);

}

CoinPackedMatrix l1TrustModel::getLinkingConstraints(int scen) {
	/* [      ]
           [ I    ]
	   [   -e ]
	 */

	vector<double> elts(2*nvar1,1.);
	vector<int> idx; idx.reserve(2*nvar1);
	vector<CoinBigIndex> starts; starts.reserve(nvar1+2);
	for (int i = 0; i < nvar1; i++) {
		starts.push_back(i);
		idx.push_back(i+1);
	}
	starts.push_back(nvar1);
	for (int i = 0; i < nvar1; i++) {
		elts[i+nvar1] = -1.;
		idx.push_back(1+nvar1+i);
	}
	starts.push_back(2*nvar1);

	return CoinPackedMatrix(true,2*nvar1+1,nvar1+1,2*nvar1,&elts[0],&idx[0],&starts[0],0);

}


