#include "ucRollingModel.hpp"
#include "CbcBALPInterface.hpp"
#include "ClpBALPInterface.hpp"

using namespace std;

int main(int argc, char **argv) {
	MPI_Init(&argc,&argv);

	int nscen = 1;
	int T = 3;
	ucRollingModel model("/sandbox/vzavala/Illinois",nscen,0,T);

	BAContext ctx(MPI_COMM_WORLD);
	ctx.initializeAssignment(nscen);

	ClpBALPInterface cbc(model, ctx);

	
	cbc.go();

	vector<double> sol = cbc.getFirstStagePrimalColSolution();
	vector<double> sol2 = cbc.getSecondStagePrimalColSolution(0);
	
	map<int,int>::const_iterator it;
	double gencost = 0;
	double oncost = 0;
	for (it = model.genMap.begin(); it != model.genMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= T; i++) {
			//printf("Pgen[%d,%d,%d] = %f, yuc[%d,%d] = %f, np_cap[%d] = %f, fuel[%d] = %s\n",i,it->first,0,
			//	sol2[model.Pgen(i,idx)],
			//	i,it->first,sol[model.yuc(i,idx)],it->first,
			//	model.genData[idx].np_cap, it->first,model.genData[idx].fuel.c_str());
			gencost += sol2[model.Pgen(i,idx)]*model.genData[idx].gen_cost;
			oncost += 1e3*sol[model.yuc(i,idx)];
		}
	}

	printf("Gen cost: %f\nOn/off cost: %f\n",gencost,oncost);
	/*
	for (it = model.busMap.begin(); it != model.busMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= T; i++) {
			printf("slackp[%d,%d] = %f, slackm[%d,%d] = %f\n",i,it->first,
			sol2[model.slackp(i,idx)],i,it->first,sol2[model.slackm(i,idx)]);
		}
	}*/
	/*for (it = model.loadMap.begin(); it != model.loadMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= T; i++) {
			printf("SOC[%d,%d] = %f, PSto[%d,%d] = %f\n",i,it->first,sol2[model.SOC(i,idx)],i,it->first,sol2[model.Psto(i,idx)]);
		}
	}*/

	/*int k = 0;
	
	for (int i = 0; i < model.nFirstStageVars()/(T+1); i++) {
		cout << i+1 << "\t";
		for (int j = 0; j <= T; j++) {
			cout << sol[k++] << "\t";
		}
		cout << endl;
	}

	cout << "-------\n";
	//cout << sol2[model.Pgen(0,236)] << "\t" << sol2[model.Pgen(1,236)] << endl;

	const double *rowact = cbc.model.getRowActivity();
	vector<string> names = model.getSecondStageRowNames(0);
	for (int i = 0; i < model.nSecondStageCons(0); i++) {
		cout << names[i] << "\t" << rowact[i] << endl;
	}*/
	MPI_Finalize();

}
