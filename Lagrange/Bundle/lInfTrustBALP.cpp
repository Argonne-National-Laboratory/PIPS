#include "lInfTrustBALP.hpp"
#include <cmath>

using namespace std;


lInfTrustModel::lInfTrustModel(int nvar1, bundle_t const &cuts, double trustRadius, std::vector<std::vector<double> > const& center) :
	nvar1(nvar1), cuts(cuts), trustRadius(trustRadius), center(center) {}

vector<double> lInfTrustModel::getFirstStageColLB() {
	vector<double> l(nvar1,-COIN_DBL_MAX);

	return l;
}

vector<double> lInfTrustModel::getFirstStageColUB() {

	return vector<double>(nvar1,COIN_DBL_MAX);
}

vector<double> lInfTrustModel::getFirstStageObj() {

	vector<double> obj(nvar1,0.);
	return obj;
}
vector<string> lInfTrustModel::getFirstStageColNames() {

	vector<string> names(nvar1,"l");
	return names;
}

vector<double> lInfTrustModel::getFirstStageRowLB() {
	return vector<double>(0);
}
vector<double> lInfTrustModel::getFirstStageRowUB() {
	return vector<double>(0);
}

vector<string> lInfTrustModel::getFirstStageRowNames() {
	return vector<string>(0);
}


vector<double> lInfTrustModel::getSecondStageColLB(int scen) {
	// >= 0, <=0 >= 0 (first two are from multipliers on box constraints for gamma)
	vector<double> lb(cuts[scen].size()+2*nvar1,0.);
	fill(lb.begin()+nvar1,lb.begin()+2*nvar1,-COIN_DBL_MIN);
	return lb;
}

vector<double> lInfTrustModel::getSecondStageColUB(int scen) {
	vector<double> ub(cuts[scen].size(),COIN_DBL_MAX);
	fill(ub.begin()+nvar1,ub.begin()+2*nvar1,0.);
	return ub;
}


vector<double> lInfTrustModel::getSecondStageObj(int scen) {
	vector<double> obj;
	obj.reserve(cuts[scen].size()+2*nvar1);
	for (int i = 0; i < nvar1; i++) {
		obj.push_back(-center[scen][i]+trustRadius);
	}
	for (int i = 0; i < nvar1; i++) {
		obj.push_back(-center[scen][i]-trustRadius);
	}
	for (unsigned i = 0; i < cuts[scen].size(); i++) {
		obj.push_back(-cuts[scen][i].computeC());
	}
	
	return obj;
}

vector<string> lInfTrustModel::getSecondStageColNames(int scen) {
	vector<string> names(cuts[scen].size()+2*nvar1,"");
	return names;
}

vector<double> lInfTrustModel::getSecondStageRowUB(int scen) {
	vector<double> ub(nvar1+1,0.);
	ub[0] = 1.;
	return ub;
}

vector<double> lInfTrustModel::getSecondStageRowLB(int scen) {
	return getSecondStageRowUB(scen);
}

vector<string> lInfTrustModel::getSecondStageRowNames(int scen) {
	vector<string> names(nvar1+1);
	names[0] = "t";
	std::fill(names.begin()+1,names.end(),"g");
	return names;

}

CoinPackedMatrix lInfTrustModel::getFirstStageConstraints() {
	vector<CoinBigIndex> starts(nvar1+2,0);
	return CoinPackedMatrix(true, 0, nvar1+1,0,0,0,&starts[0],0);
}

CoinPackedMatrix lInfTrustModel::getSecondStageConstraints(int scen) {
	/* diagonal blocks are of the form:
	   [       e^T   ] 
	   [ I  I -G_i^T ]
	   This is convenient for updating the basis when subgradients are added
	 */
	
	CoinBigIndex nnz = 0;
	unsigned ncol = cuts[scen].size() + 2*nvar1;
	vector<double> elts; elts.reserve(cuts[scen].size()*(nvar1+1)+2*nvar1);
	vector<int> idx; idx.reserve(cuts[scen].size()*(nvar1+1)+2*nvar1);
	vector<CoinBigIndex> starts(ncol+1);
	for (int i = 0; i < nvar1; i++) {
		elts.push_back(1.);
		idx.push_back(i+1);
		starts[i] = nnz++;
	}
	for (int i = 0; i < nvar1; i++) {
		elts.push_back(1.);
		idx.push_back(i+1);
		starts[i+nvar1] = nnz++;
	}
	for (unsigned col = 0; col < cuts[scen].size(); col++) {
		elts.push_back(1.);
		idx.push_back(0);
		starts[col+2*nvar1] = nnz++;
		vector<double> const& subgrad = cuts[scen][col].subgradient;
		for (unsigned k = 0; k < subgrad.size(); k++) {
			if (fabs(subgrad[k]) > 1e-7) {
				elts.push_back(-subgrad[k]);
				idx.push_back(k+1);
				nnz++;
			}
		}
	}

	starts[ncol] = nnz;

	return CoinPackedMatrix(true,nvar1+1,ncol,nnz,&elts[0],&idx[0],&starts[0],0);

}

CoinPackedMatrix lInfTrustModel::getLinkingConstraints(int scen) {
	// [   ]
	// [ I ]

	vector<double> elts(nvar1,1.);
	vector<int> idx; idx.reserve(nvar1);
	vector<CoinBigIndex> starts; starts.reserve(nvar1+1);
	for (int i = 0; i < nvar1; i++) {
		starts.push_back(i);
		idx.push_back(i+1);
	}
	starts.push_back(nvar1);

	return CoinPackedMatrix(true,nvar1+1,nvar1,nvar1,&elts[0],&idx[0],&starts[0],0);

}


