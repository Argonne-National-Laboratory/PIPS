#include "ucRollingModel.hpp"
#include <sstream>
#include <fstream>
#include <cstring>

using namespace std;

ucRollingModel::ucRollingModel(string const& dataRoot, int nscen, int tOffset, int tHorizon,
		MPI_Comm comm) {
	givenInitial = false;
	readData(dataRoot, comm);

}

namespace {
// read file on proc 0 and broadcast to all
string readFile(string const& fname, MPI_Comm comm) { 
	unsigned filelen;
	char *filedata;
	int mype;
	MPI_Comm_rank(comm,&mype);

	if (mype == 0) {
		ifstream f(fname.c_str());
		assert(f.is_open());
		
		string data((istreambuf_iterator<char>(f)),istreambuf_iterator<char>()); // read into string
		f.close();
		filelen = data.length() + 1;
		filedata = new char[filelen];
		memcpy(filedata,data.c_str(),filelen);
		MPI_Bcast(&filelen,1,MPI_UNSIGNED,0,comm);
	} else {
		MPI_Bcast(&filelen,1,MPI_UNSIGNED,0,comm);
		filedata = new char[filelen];
	}
	MPI_Bcast(filedata,filelen,MPI_CHAR,0,comm);
	
	string data(filedata);
	delete [] filedata;

	return data;
}
}

void ucRollingModel::readData(string const& dataRoot, MPI_Comm comm) {

	string line;

	{
		istringstream ss(readFile(dataRoot+"/Lines_data.tab",comm));
		getline(ss,line);
		getline(ss,line);
		while (!getline(ss,line).eof()) {
			istringstream ls(line);
			ls.exceptions(ifstream::failbit | ifstream::badbit);
			linStruct l;
			l.idx = linData.size();
			string name, status;
			ls >> name;
			ls >> status;
			ls >> l.snd_bus;
			ls >> l.rec_bus;
			ls >> l.X;
			ls >> l.V;
			ls >> l.Pmax;
			linData.insert(pair<string,linStruct>(name,l));
		}
	}

	{
		istringstream ss(readFile(dataRoot+"/bus_data.tab",comm));
		getline(ss,line);
		getline(ss,line);
		while (!getline(ss,line).eof()) {
			istringstream ls(line);
			ls.exceptions(ifstream::failbit | ifstream::badbit);
			busStruct b;
			b.idx = busData.size();
			int name;
			ls >> name;
			busData.insert(pair<int,busStruct>(name,b));
		}
	}

	{
		istringstream ss(readFile(dataRoot+"/Gen_data.tab",comm));
		getline(ss,line);
		getline(ss,line);
		while (!getline(ss,line).eof()) {
			istringstream ls(line);
			ls.exceptions(ifstream::failbit | ifstream::badbit);
			genStruct g;
			g.idx = genData.size();
			int name; string alias;
			ls >> name;
			ls >> alias;
			ls >> g.bus_gen;
			ls >> g.np_cap;
			ls >> g.sum_cap;
			ls >> g.win_cap;
			ls >> g.min_cap;
			ls >> g.fuel;
			ls >> g.min_hrate;
			ls >> g.min_power;
			ls >> g.max_ur;
			ls >> g.max_dr;
			genData.insert(pair<int,genStruct>(name,g));
		}
	}

	{
		istringstream ss(readFile(dataRoot+"/fuel_data.tab",comm));
		getline(ss,line);
		getline(ss,line);
		while (!getline(ss,line).eof()) {
			istringstream ls(line);
			ls.exceptions(ifstream::failbit | ifstream::badbit);
			fuelStruct f;
			f.idx = fuelData.size();
			string name, desc;
			ls >> name;
			ls >> desc;
			ls >> f.HV;
			ls >> f.Unitprice;
			fuelData.insert(pair<string,fuelStruct>(name,f));
		}
	}

	{
		istringstream ss(readFile(dataRoot+"/load_load.tab",comm));
		getline(ss,line);
		getline(ss,line);
		while (!getline(ss,line).eof()) {
			istringstream ls(line);
			ls.exceptions(ifstream::failbit | ifstream::badbit);
			loadStruct l;
			l.idx = loadData.size();
			int name;
			ls >> name;
			ls >> l.bus_load;
			loadData.insert(pair<int,loadStruct>(name,l));
		}
	}
	
	// todo: load data

	for (std::map<int,genStruct>::iterator it = genData.begin(); it != genData.end(); ++it) {
		genStruct &g = it->second;
		std::map<string,fuelStruct>::iterator fit = fuelData.find(g.fuel);
		assert(fit != fuelData.end());
		fuelStruct const& f = fit->second;
		g.gen_cost = 1e-3*(g.min_hrate/f.HV)*f.Unitprice;
	}

	{
		istringstream ss(readFile(dataRoot+"/wind_data.tab",comm));
		getline(ss,line);
		getline(ss,line);
		while (!getline(ss,line).eof()) {
			istringstream ls(line);
			ls.exceptions(ifstream::failbit | ifstream::badbit);
			windStruct w;
			w.idx = windData.size();
			int name;
			ls >> name;
			ls >> w.bus_wind;
			ls >> w.wind_share;
			windData.insert(pair<int,windStruct>(name,w));
		}
	}

	
}



vector<double> ucRollingModel::getFirstStageColLB() {
	return vector<double>(0);
}

vector<double> ucRollingModel::getFirstStageColUB() {}
vector<double> ucRollingModel::getFirstStageObj() {}
vector<string> ucRollingModel::getFirstStageColNames() {}
vector<double> ucRollingModel::getFirstStageRowLB() {}
vector<double> ucRollingModel::getFirstStageRowUB() {}
vector<string> ucRollingModel::getFirstStageRowNames() {}

vector<double> ucRollingModel::getSecondStageColLB(int scen) {}
vector<double> ucRollingModel::getSecondStageColUB(int scen) {}
vector<double> ucRollingModel::getSecondStageObj(int scen) {}
vector<string> ucRollingModel::getSecondStageColNames(int scen) {}
vector<double> ucRollingModel::getSecondStageRowUB(int scen) {}
vector<double> ucRollingModel::getSecondStageRowLB(int scen) {}
vector<string> ucRollingModel::getSecondStageRowNames(int scen) {}

CoinPackedMatrix ucRollingModel::getFirstStageConstraints() {}
CoinPackedMatrix ucRollingModel::getSecondStageConstraints(int scen) {}
CoinPackedMatrix ucRollingModel::getLinkingConstraints(int scen) {}



