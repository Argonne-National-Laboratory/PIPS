#include "ucRollingModel.hpp"
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <limits>

#include "CoinPackedVector.hpp"

using namespace std;

const double THETASCALE = 1e5;//1e5;
const double OBJSCALE = 1e2;

ucRollingModel::ucRollingModel(string const& dataRoot, int nscen, int tOffset, int tHorizon,
		MPI_Comm comm) {
	this->nscen = nscen;
	timeOffset = tOffset;
	horizon = tHorizon;
	givenInitial = false;
	sigma = 5.;
	readData(dataRoot, comm);
	initializeVariables();
	

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
		linData.reserve(2522);
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
			linMap.insert(pair<string,int>(name,linData.size()));
			linData.push_back(l);
		}
	}

	{
		istringstream ss(readFile(dataRoot+"/bus_data.tab",comm));
		getline(ss,line);
		getline(ss,line);
		while (!getline(ss,line).eof()) {
			istringstream ls(line);
			ls.exceptions(ifstream::failbit | ifstream::badbit);
			int name;
			ls >> name;
			busMap.insert(pair<int,int>(name,busMap.size()));
		}
	}

	{
		istringstream ss(readFile(dataRoot+"/Gen_data.tab",comm));
		getline(ss,line);
		getline(ss,line);
		genData.reserve(261);
		while (!getline(ss,line).eof()) {
			istringstream ls(line);
			ls.exceptions(ifstream::failbit | ifstream::badbit);
			genStruct g;
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
			genMap.insert(pair<int,int>(name,genData.size()));
			genData.push_back(g);
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
			string name, desc;
			ls >> name;
			ls >> desc;
			ls >> f.HV;
			ls >> f.Unitprice;
			fuelMap.insert(pair<string,int>(name,fuelData.size()));
			fuelData.push_back(f);
		}
	}

	{
		istringstream ss(readFile(dataRoot+"/load_load.tab",comm));
		getline(ss,line);
		getline(ss,line);
		loadData.reserve(870);
		while (!getline(ss,line).eof()) {
			istringstream ls(line);
			ls.exceptions(ifstream::failbit | ifstream::badbit);
			loadStruct l;
			int name;
			ls >> name;
			ls >> l.bus_load;
			loadMap.insert(pair<int,int>(name,loadData.size()));
			loadData.push_back(l);
		}
	}
	
	// for now, read loads for all time periods
	{
		istringstream ss(readFile(dataRoot+"/Loads.dat",comm));
		int t = 0;
		loads.reserve(8760); // should there be 8761 lines?
		while (!getline(ss,line).eof()) {
			istringstream ls(line);
			ls.exceptions(ifstream::failbit | ifstream::badbit);
			loads.push_back(vector<double>());
			loads[t].resize(loadData.size());
			for (unsigned i = 0; i < loadData.size(); i++) {
				ls >> loads[t][i];
			}
			t++;
		}
	}
	for (unsigned i = 0; i < genData.size(); i++) {
		genStruct &g = genData[i];
		std::map<string,int>::iterator fit = fuelMap.find(g.fuel);
		assert(fit != fuelMap.end());
		fuelStruct const& f = fuelData[fit->second];
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
			int name;
			ls >> name;
			ls >> w.bus_wind;
			ls >> w.wind_share;
			windMap.insert(pair<int,int>(name,windData.size()));
			windData.push_back(w);
		}
	}

	{
		istringstream ss(readFile(dataRoot+"/wind_data.dat",comm));
		double w;
		//ss >> w; // seems that first one should be skipped??
		ss >> w;
		while (!ss.fail()) {
			wind_total_determ.push_back(w);
			ss >> w;
		}
		//assert(wind_total_determ.size() == loads.size());
	}

	{
		istringstream ss(readFile(dataRoot+"/SOCprofile.dat",comm));
		double w;
		ss >> w;
		while (!ss.fail()) {
			SOCprof_determ.push_back(w);
			ss >> w;
		}
		assert(loads.size() == SOCprof_determ.size());
	}

	for (unsigned i = 0; i < loadData.size(); i++) {
		loadData[i].SOCcap = (0.05*10e6*2.0/1000.0)/loadData.size();
		loadData[i].SOC_init = loadData[i].SOCcap*SOCprof_determ[0];
	}
				
	wind_total.resize(nscen);

	cout << genData.size() << " generators, " << linData.size() << " lines, " << busMap.size()  << " buses, " << loadData.size() << " loads\n";

}

void ucRollingModel::initializeVariables() {

	yuc.initialize(horizon+1,0,genMap);
	nvar1 = yuc.totalVars();

	nvar2 = 0;
	Pgen.initialize(horizon+1,nvar2,genMap);
	nvar2 += Pgen.totalVars();

	//dPgen.initialize(horizon,nvar2,genMap);
	//nvar2 += dPgen.totalVars();

	Pwind.initialize(horizon+1,nvar2,windMap);
	nvar2 += Pwind.totalVars();

	theta.initialize(horizon+1,nvar2,busMap);
	nvar2 += theta.totalVars();

	slackp.initialize(horizon+1,nvar2,busMap);
	nvar2 += slackp.totalVars();

	slackm.initialize(horizon+1,nvar2,busMap);
	nvar2 += slackm.totalVars();

	P.initialize(horizon+1,nvar2,linMap);
	nvar2 += P.totalVars();

	SOC.initialize(horizon+1,nvar2,loadMap);
	nvar2 += SOC.totalVars();

	Psto.initialize(horizon+1,nvar2,loadMap);
	nvar2 += Psto.totalVars();

	ncons2 = 0;
	ncons2 += (horizon+1)*busMap.size(); // pfeq
	ncons2 += (horizon+1)*linMap.size(); // Peq
	ncons2 += (horizon+1)*genMap.size(); //Pgenequb
	ncons2 += horizon*genMap.size(); // rampdynamics
	ncons2 += horizon*loadMap.size(); // stodynamics
	ncons2 += loadMap.size(); // stoTeq


}

#include <stdint.h>

namespace{
// taken from rand() manpage
// so that generated sequnece is reproducible
/* RAND_MAX assumed to be 32767 */
double  myrand(uint64_t &next) {
	next = next * 1103515245 + 12345;
	return((uint32_t)(next/65536) % 32768)/double(32767);
}

// polar method
double gaussian(double mean, double stdev, uint64_t &state) {
	double u,v,s;
	do {
		u = myrand(state)*2.-1.;
		v = myrand(state)*2.-1.;
		s = u*u + v*v;
	} while (s >= 1 || s == 0.0);
	// throw away extra
	return mean + stdev*u*sqrt(-2.0*log(s)/s);

}
}

// sigma: standard deviation as a percent of the mean
void ucRollingModel::generateWind(int scen, double sigma) {
	vector<double> &v = wind_total[scen];
	if (v.size()) return; // already generated
	if (scen == 0) { v = wind_total_determ; return; }
	uint64_t seed = scen;
	for (unsigned i = 0; i < wind_total_determ.size(); i++) {
		v.push_back(gaussian(wind_total_determ[i],wind_total_determ[i]*sigma/100.,seed));
	}


}

vector<double> ucRollingModel::getFirstStageColLB() {
	vector<double> lb(nvar1,0.0);
	// leave baseload plants ON
	for (map<int,int>::const_iterator it = genMap.begin(); it != genMap.end(); ++it) {
		int idx = it->second;
		const genStruct &g = genData[idx];
		if (g.fuel == "Other" || (g.np_cap >= 500 && (g.fuel == "UR" || g.fuel == "BIT" || g.fuel == "SUB"))) {
			for (int i = 0; i <= horizon; i++) {
				lb[yuc(i,idx)] = 1.0;
			}
		}
	}
	return lb;
}

vector<double> ucRollingModel::getFirstStageColUB() {
	return vector<double>(nvar1,1.0);
}


vector<double> ucRollingModel::getFirstStageObj() { 
	return vector<double>(nvar1,OBJSCALE*1e3); // TODO: vary per generator
}

namespace{
template<typename T> string toStr(string const& name, int t, T key) {
	stringstream str;
	str << name << "[" << t << "," << key << "]";
	return str.str();
}
}

vector<string> ucRollingModel::getFirstStageColNames() {
	vector<string> s(nvar1);
	for (map<int,int>::const_iterator it = genMap.begin(); it != genMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= horizon; i++) {
			s[yuc(i,idx)] = toStr("yuc",i,it->first);
		}
	}
	return s;
}

vector<double> ucRollingModel::getFirstStageRowLB() {
	return vector<double>(0);
}

vector<double> ucRollingModel::getFirstStageRowUB() { 
	return vector<double>(0);
}

vector<string> ucRollingModel::getFirstStageRowNames() {
	return vector<string>(0);
}

vector<double> ucRollingModel::getSecondStageColLB(int scen) {
	generateWind(scen,sigma);
	vector<double> lb(nvar2);

	map<int,int>::const_iterator it;

	for (it = genMap.begin(); it != genMap.end(); ++it) {
		int idx = it->second;
		const genStruct &g = genData[idx];
		for (int i = 0; i <= horizon; i++) {
			lb[Pgen(i,idx)] = 0.;
		}
		if (givenInitial) lb[Pgen(0,idx)] = g.Pgen_init;
	}

	for (it = windMap.begin(); it != windMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= horizon; i++) {
			lb[Pwind(i,idx)] = 0.0;
		}
	}
	
	for (it = busMap.begin(); it != busMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= horizon; i++) {
			lb[theta(i,idx)] = -THETASCALE*3.14/2.;
			lb[slackp(i,idx)] = 0.;
			lb[slackm(i,idx)] = 0.;
		}
	}
	// reference bus
	for (int i = 0; i <= horizon; i++) {
		lb[theta.lookup(i,22671)] = 0.;
	}
	
	for (map<string,int>::const_iterator it = linMap.begin(); it != linMap.end(); ++it) {
		int idx = it->second;
		const linStruct &l = linData[idx];
		for (int i = 0; i <= horizon; i++) {
			lb[P(i,idx)] = -l.Pmax;
		}
	}

	for (it = loadMap.begin(); it != loadMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= horizon; i++) {
			lb[SOC(i,idx)] = 0.;
			lb[Psto(i,idx)] = -numeric_limits<double>::infinity();
		}
		lb[SOC(0,idx)] = loadData[idx].SOC_init;
	}

	return lb;

}
vector<double> ucRollingModel::getSecondStageColUB(int scen) {
	generateWind(scen,sigma);

	vector<double> ub(nvar2);

	map<int,int>::const_iterator it;

	for (it = genMap.begin(); it != genMap.end(); ++it) {
		int idx = it->second;
		const genStruct &g = genData[idx];
		for (int i = 0; i <= horizon; i++) {
			ub[Pgen(i,idx)] = g.np_cap; // *yuc
		}
		if (givenInitial) ub[Pgen(0,idx)] = g.Pgen_init;
	}

	for (it = windMap.begin(); it != windMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= horizon; i++) {
			ub[Pwind(i,idx)] = windData[idx].wind_share*wind_total[scen][i+timeOffset];
		}
	}
	
	for (it = busMap.begin(); it != busMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= horizon; i++) {
			ub[theta(i,idx)] = THETASCALE*3.14/2.;
			ub[slackp(i,idx)] = numeric_limits<double>::infinity();
			ub[slackm(i,idx)] = numeric_limits<double>::infinity();
		}
	}
	// reference bus
	for (int i = 0; i <= horizon; i++) {
		ub[theta.lookup(i,22671)] = 0.;
	}
	
	for (map<string,int>::const_iterator it = linMap.begin(); it != linMap.end(); ++it) {
		int idx = it->second;
		const linStruct &l = linData[idx];
		for (int i = 0; i <= horizon; i++) {
			ub[P(i,idx)] = l.Pmax;
		}
	}

	for (it = loadMap.begin(); it != loadMap.end(); ++it) {
		int idx = it->second;
		loadStruct const &l = loadData[idx];
		for (int i = 0; i <= horizon; i++) {
			ub[SOC(i,idx)] = l.SOCcap*SOCprof_determ[i+timeOffset]; // in AMPL SOCprof data is per scenario but actually all the same
			ub[Psto(i,idx)] = numeric_limits<double>::infinity();
		}
		ub[SOC(0,idx)] = loadData[idx].SOC_init;
	}

	return ub;
}

vector<double> ucRollingModel::getSecondStageObj(int scen) {

	vector<double> obj(nvar2, 0.0);

	map<int,int>::const_iterator it;

	for (it = genMap.begin(); it != genMap.end(); ++it) {
		int idx = it->second;
		const genStruct &g = genData[idx];
		for (int i = 0; i <= horizon; i++) {
			obj[Pgen(i,idx)] = OBJSCALE*g.gen_cost/nscen;
		}
	}

	for (it = busMap.begin(); it != busMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= horizon; i++) {
			obj[slackp(i,idx)] = OBJSCALE*1e3/nscen;
			obj[slackm(i,idx)] = OBJSCALE*1e3/nscen;
		}
	}
	
	return obj;

}

namespace{
template<typename T> string toStr(string const& name, int t, T key, int scen) {
	stringstream str;
	str << name << "[" << t << "," << key << "] (" << scen << ")";
	return str.str();
}
}

vector<string> ucRollingModel::getSecondStageColNames(int scen) {
	vector<string> s(nvar2);
	
	map<int,int>::const_iterator it;
	for (it = genMap.begin(); it != genMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= horizon; i++) {
			s[Pgen(i,idx)] = toStr("Pgen",i,it->first,scen);
			//s[dPgen(i,idx)] = toStr("dPgen",i,it->first,scen);
		}
	}

	for (it = windMap.begin(); it != windMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= horizon; i++) {
			s[Pwind(i,idx)] = toStr("Pwind",i,it->first,scen);
		}
	}
	
	for (it = busMap.begin(); it != busMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= horizon; i++) {
			s[theta(i,idx)] = toStr("theta",i,it->first,scen);
			s[slackp(i,idx)] = toStr("slackp",i,it->first,scen);
			s[slackm(i,idx)] = toStr("slackm",i,it->first,scen);
		}
	}
	
	for (map<string,int>::const_iterator it = linMap.begin(); it != linMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= horizon; i++) {
			s[P(i,idx)] = toStr("P",i,it->first,scen);
		}
	}

	for (it = loadMap.begin(); it != loadMap.end(); ++it) {
		int idx = it->second;
		for (int i = 0; i <= horizon; i++) {
			s[SOC(i,idx)] = toStr("SOC",i,it->first,scen);
			s[Psto(i,idx)] = toStr("Psto",i,it->first,scen);
		}
	}

	return s;

}


vector<double> ucRollingModel::getSecondStageRowUB(int scen) {
	generateWind(scen,sigma);

	vector<double> ub(ncons2);
	
	map<int,int>::const_iterator it;
	
	int r = 0;

	// pfeq
	for (it = busMap.begin(); it != busMap.end(); ++it) {
		//cout << it->first << "\t";
		for (int i = 0; i <= horizon; i++) {
			// = load on this bus - wind to this bus
			double sum = 0.;
			for (unsigned loadidx = 0; loadidx < loadData.size(); loadidx++) {
				if (loadData[loadidx].bus_load == it->first) {
					sum += loads[timeOffset+i][loadidx];
				}
			}
			for (unsigned windidx = 0; windidx < windData.size(); windidx++) {
				if (windData[windidx].bus_wind == it->first) {
					sum -= wind_total[scen][timeOffset+i]*windData[windidx].wind_share;
				}
			}
			ub[r++] = sum;
			//cout << sum << "\t";

		}
		//cout << endl;
	}	
	
	// Peq
	for (unsigned linidx = 0; linidx < linData.size(); linidx++) {
		for (int i = 0; i <= horizon; i++) {
			ub[r++] = 0.;
		}
	}
	
	// Pgenequb
	for (unsigned genidx = 0; genidx < genData.size(); genidx++) {
		for (int i = 0; i <= horizon; i++) {
			ub[r++] = 0.;
		}
	}

	// rampdynamics: -max_dr[i] <=  Pgen[s,t+1,i] - Pgen[s,t,i] <= max_ur[i]
	// we've presolved out dPgen
	for (it = genMap.begin(); it != genMap.end(); ++it) {
		int idx = it->second;
		const genStruct &g = genData[idx];
		for (int i = 0; i < horizon; i++) {
			ub[r++] = g.max_ur;
		}
	}

	// stodynamics
	for (unsigned loadidx = 0; loadidx < loadData.size(); loadidx++) {
		for (int i = 0; i < horizon; i++) {
			ub[r++] = 0.;
		}
	}

	// stoTeq
	for (it = loadMap.begin(); it != loadMap.end(); ++it) {
		int idx = it->second;
		loadStruct const &l = loadData[idx];
		ub[r++] = l.SOCcap*SOCprof_determ[horizon+1+timeOffset];
	}

	assert(r == ncons2);

	return ub;

}

vector<double> ucRollingModel::getSecondStageRowLB(int scen) {
	generateWind(scen,sigma);

	vector<double> lb(ncons2);
	
	map<int,int>::const_iterator it;
	
	int r = 0;
	/*
	for (int i = 0; i <= horizon; i++) {
		cout << "wind_total[" << i << "] = " << wind_total[0][i] << endl;
	}*/

	// pfeq
	for (it = busMap.begin(); it != busMap.end(); ++it) {
		for (int i = 0; i <= horizon; i++) {
			// = load on this bus - wind to this bus
			double sum = 0.;
			for (unsigned loadidx = 0; loadidx < loadData.size(); loadidx++) {
				if (loadData[loadidx].bus_load == it->first) {
					sum += loads[timeOffset+i][loadidx];
				}
			}
			for (unsigned windidx = 0; windidx < windData.size(); windidx++) {
				if (windData[windidx].bus_wind == it->first) {
					sum -= wind_total[scen][timeOffset+i]*windData[windidx].wind_share;
				}
			}
			lb[r++] = sum;
		}
	}	
	
	// Peq
	for (unsigned linidx = 0; linidx < linData.size(); linidx++) {
		for (int i = 0; i <= horizon; i++) {
			lb[r++] = 0.;
		}
	}
	
	// Pgenequb
	for (it = genMap.begin(); it != genMap.end(); ++it) {
		int idx = it->second;
		genStruct const &g = genData[idx];
		for (int i = 0; i <= horizon; i++) {
			lb[r++] = -g.np_cap;
		}
	}

	// rampdynamics
	for (it = genMap.begin(); it != genMap.end(); ++it) {
		int idx = it->second;
		const genStruct &g = genData[idx];
		for (int i = 0; i < horizon; i++) {
			lb[r++] = -g.max_dr;
		}
	}

	// stodynamics
	for (unsigned loadidx = 0; loadidx < loadData.size(); loadidx++) {
		for (int i = 0; i < horizon; i++) {
			lb[r++] = 0.;
		}
	}

	// stoTeq
	for (it = loadMap.begin(); it != loadMap.end(); ++it) {
		lb[r++] = 0.;
	}
	
	assert(r == ncons2);

	return lb;

}

vector<string> ucRollingModel::getSecondStageRowNames(int scen) {

	vector<string> s(ncons2);

	map<int,int>::const_iterator it;
	
	int r = 0;

	// pfeq
	for (it = busMap.begin(); it != busMap.end(); ++it) {
		for (int i = 0; i <= horizon; i++) {
			s[r++] = toStr("pfeq",i,it->first,scen);
		}
	}	
	
	// Peq
	for (map<string,int>::const_iterator linit = linMap.begin(); linit != linMap.end(); ++linit) {
		for (int i = 0; i <= horizon; i++) {
			s[r++] = toStr("Peq",i,linit->first,scen);
		}
	}
	
	// Pgenequb
	for (it = genMap.begin(); it != genMap.end(); ++it) {
		for (int i = 0; i <= horizon; i++) {
			s[r++] = toStr("Pgenequb",i,it->first,scen);
		}
	}

	// rampdynamics
	for (it = genMap.begin(); it != genMap.end(); ++it) {
		for (int i = 0; i < horizon; i++) {
			s[r++] = toStr("rampdyanmics",i,it->first,scen);
		}
	}

	// stodynamics
	for (it = loadMap.begin(); it != loadMap.end(); ++it) {
		for (int i = 0; i < horizon; i++) {
			s[r++] = toStr("stodynamics",i,it->first,scen);
		}
	}

	// stoTeq
	for (it = loadMap.begin(); it != loadMap.end(); ++it) {
		s[r++] = toStr("stoTeq",horizon,it->first,scen);
	}

	assert(r == ncons2);

	return s;

}

CoinPackedMatrix ucRollingModel::getFirstStageConstraints() {
	vector<CoinBigIndex> starts(nvar1+1,0);
	return CoinPackedMatrix(true, 0, nvar1,0,0,0,&starts[0],0);

}

CoinPackedMatrix ucRollingModel::getSecondStageConstraints(int scen) {
	
	vector<CoinPackedVectorBase*> rows(ncons2);
	map<int,int>::const_iterator it;

	int r = 0;

	// pfeq
	for (it = busMap.begin(); it != busMap.end(); ++it) {
		for (int i = 0; i <= horizon; i++) {
			vector<double> elts;
			vector<int> idx;
			for (unsigned linidx = 0; linidx < linData.size(); linidx++) {
				if (linData[linidx].rec_bus == it->first) {
					elts.push_back(1.0);
					idx.push_back(P(i,linidx));
				}
				if (linData[linidx].snd_bus == it->first) {
					elts.push_back(-1.0);
					idx.push_back(P(i,linidx));
				}
			}
			for (unsigned genidx = 0; genidx < genData.size(); genidx++) {
				if (genData[genidx].bus_gen == it->first) {
					elts.push_back(1.0);
					idx.push_back(Pgen(i,genidx));
				}
			}
			for (unsigned windidx = 0; windidx < windData.size(); windidx++) {
				if (windData[windidx].bus_wind == it->first) {
					elts.push_back(-1.0);
					idx.push_back(Pwind(i,windidx));
				}
			}
			for (unsigned loadidx = 0; loadidx < loadData.size(); loadidx++) {
				if (loadData[loadidx].bus_load == it->first) {
					elts.push_back(1.0);
					idx.push_back(Psto(i,loadidx));
				}
			}
			elts.push_back(1.0);
			idx.push_back(slackm(i,it->second));
			elts.push_back(-1.0);
			idx.push_back(slackp(i,it->second));
			rows[r++] = new CoinPackedVector(elts.size(),&idx[0],&elts[0]);
		}
	}	
	
	// Peq
	for (map<string,int>::const_iterator linit = linMap.begin(); linit != linMap.end(); ++linit) {
		linStruct const &l = linData[linit->second];
		for (int i = 0; i <= horizon; i++) {
			double elts[] = {l.X, -l.V*l.V/THETASCALE, l.V*l.V/THETASCALE};
			int idx[] = {P(i,linit->second),theta.lookup(i,l.snd_bus),theta.lookup(i,l.rec_bus)};
			rows[r++] = new CoinPackedVector(3,idx,elts);
		}
	}
	
	// Pgenequb -- linking constraint
	for (it = genMap.begin(); it != genMap.end(); ++it) {
		for (int i = 0; i <= horizon; i++) {
			double one = 1.0;
			int idx = Pgen(i,it->second);
			rows[r++] = new CoinPackedVector(1,&idx,&one);
		}
	}

	// rampdynamics
	for (it = genMap.begin(); it != genMap.end(); ++it) {
		for (int i = 0; i < horizon; i++) {
			double elts[] = {1.0,-1.0};
			int idx[] = {Pgen(i+1,it->second),Pgen(i,it->second)};
			rows[r++] = new CoinPackedVector(2,idx,elts);
		}
	}

	// stodynamics
	for (it = loadMap.begin(); it != loadMap.end(); ++it) {
		for (int i = 0; i < horizon; i++) {
			double elts[] = {1.0,-1.0,-1.0};
			int idx[] = {SOC(i+1,it->second),SOC(i,it->second),Psto(i,it->second)};
			rows[r++] = new CoinPackedVector(3,idx,elts);
		}
	}

	// stoTeq
	for (it = loadMap.begin(); it != loadMap.end(); ++it) {
		double elts[] = {1.0,1.0};
		int idx[] = {SOC(horizon,it->second),Psto(horizon,it->second)};
		rows[r++] = new CoinPackedVector(2,idx,elts);
	}

	assert(r == ncons2);

	CoinPackedMatrix mat;
	assert(mat.isColOrdered());
	mat.setDimensions(0,nvar2);
	mat.appendRows(ncons2,&rows[0]);

	for (int i = 0; i < ncons2; i++) delete rows[i];

	return mat;

}

CoinPackedMatrix ucRollingModel::getLinkingConstraints(int scen) {
	
	vector<CoinPackedVectorBase*> rows(ncons2);
	map<int,int>::const_iterator it;

	int r = 0;

	// pfeq
	for (it = busMap.begin(); it != busMap.end(); ++it) {
		for (int i = 0; i <= horizon; i++) {
			rows[r++] = new CoinPackedVector();
		}
	}	
	
	// Peq
	for (map<string,int>::const_iterator linit = linMap.begin(); linit != linMap.end(); ++linit) {
		for (int i = 0; i <= horizon; i++) {
			rows[r++] = new CoinPackedVector();
		}
	}
	
	// Pgenequb -- linking constraint
	for (it = genMap.begin(); it != genMap.end(); ++it) {
		for (int i = 0; i <= horizon; i++) {
			double cap = -genData[it->second].np_cap;
			int idx = yuc(i,it->second);
			rows[r++] = new CoinPackedVector(1,&idx,&cap);
		}
	}

	// rampdynamics
	for (it = genMap.begin(); it != genMap.end(); ++it) {
		for (int i = 0; i < horizon; i++) {
			rows[r++] = new CoinPackedVector();
		}
	}

	// stodynamics
	for (it = loadMap.begin(); it != loadMap.end(); ++it) {
		for (int i = 0; i < horizon; i++) {
			rows[r++] = new CoinPackedVector();
		}
	}

	// stoTeq
	for (it = loadMap.begin(); it != loadMap.end(); ++it) {
		rows[r++] = new CoinPackedVector();
	}

	assert(r == ncons2);

	CoinPackedMatrix mat;
	assert(mat.isColOrdered());
	mat.setDimensions(0,nvar1);
	mat.appendRows(ncons2,&rows[0]);

	assert(mat.getNumRows() == ncons2);
	
	for (int i = 0; i < ncons2; i++) delete rows[i];

	return mat;

}



