#include "levelBALP.hpp"
#include <cmath>


using namespace std;


levelModel::levelModel(int nvar1, std::vector<std::vector<cutInfo> > const &cuts, double l, std::vector<std::vector<double> > const& center) :
	nvar1(nvar1), cuts(cuts), center(center), l(l) {}

vector<double> levelModel::getFirstStageColLB() {
	// v >= 0, l free
	vector<double> l(nvar1 + 1);
	l[0] = 0.;
	std::fill(l.begin()+1,l.end(),-COIN_DBL_MAX);

	return l;
}

vector<double> levelModel::getFirstStageColUB() {

	return vector<double>(nvar1+1,COIN_DBL_MAX);
}

vector<double> levelModel::getFirstStageObj() {

	vector<double> obj(nvar1+1,0.);
	obj[0] = l;
	return obj;
}
vector<string> levelModel::getFirstStageColNames() {

	vector<string> names(nvar1+1);
	names[0] = "v";
	std::fill(names.begin()+1,names.end(),"l"); // don't worry about numbers for now
	return names;
}

vector<double> levelModel::getFirstStageRowLB() {
	return vector<double>(0);
}
vector<double> levelModel::getFirstStageRowUB() {
	return vector<double>(0);
}

vector<string> levelModel::getFirstStageRowNames() {
	return vector<string>(0);
}


vector<double> levelModel::getSecondStageColLB(int scen) {
	return vector<double>(nSecondStageVars(scen),0.);
}

vector<double> levelModel::getSecondStageColUB(int scen) {
	return vector<double>(nSecondStageVars(scen),COIN_DBL_MAX);
}


vector<double> levelModel::getSecondStageObj(int scen) {
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

vector<string> levelModel::getSecondStageColNames(int scen) {
	vector<string> names(nSecondStageVars(scen));
	return names;
}

vector<double> levelModel::getSecondStageRowUB(int scen) {
	vector<double> ub(nSecondStageCons(scen),0.);
	for (int i = 0; i < nvar1; i++) {
		ub[nvar1+1+i] = 1.;
	}
	return ub;
}

vector<double> levelModel::getSecondStageRowLB(int scen) {
	return getSecondStageRowUB(scen);
}

vector<string> levelModel::getSecondStageRowNames(int scen) {
	vector<string> names(nSecondStageCons(scen));
	return names;

}

CoinPackedMatrix levelModel::getFirstStageConstraints() {
	vector<CoinBigIndex> starts(nvar1+2,0);
	return CoinPackedMatrix(true, 0, nvar1+1,0,0,0,&starts[0],0);
}

CoinPackedMatrix levelModel::getSecondStageConstraints(int scen) {
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

CoinPackedMatrix levelModel::getLinkingConstraints(int scen) {
	/* [ -1   ]
           [    I ]
	   [      ]
	 */

	vector<double> elts(nvar1+1,1.);
	elts[0] = -1.;
	vector<int> idx; idx.reserve(nvar1+1);
	vector<CoinBigIndex> starts; starts.reserve(nvar1+2);
	for (int i = 0; i < nvar1+1; i++) {
		starts.push_back(i);
		idx.push_back(i);
	}
	starts.push_back(nvar1+1);

	return CoinPackedMatrix(true,2*nvar1+1,nvar1+1,nvar1+1,&elts[0],&idx[0],&starts[0],0);

}


