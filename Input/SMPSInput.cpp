#include <sstream>
#include <fstream>
#include <cmath>

#include "SMPSInput.hpp"

using namespace std;

namespace {
// helper function for copying bools and strings from mps file
template<typename T1, typename T2> void copySubset(
	CoinMpsIO const &reader, T2 (CoinMpsIO::*f)(int) const, 
	int from, int to, vector<T1> &out) {
	
	unsigned n = to - from;
	assert(out.size() == n);
	for (unsigned i = 0; i < n; i++) {
		out[i] = (reader.*f)(i+from);
	}
}
}

// parse input files
// minimal error checking, 
// but set exceptions on ifstream so failure isn't silent.
// Some asserts should really be exceptions.
SMPSInput::SMPSInput(string const& cor, string const& tim, string const& sto) :
	corfile(cor), timfile(tim), stofile(sto) {

	reader.readMps(cor.c_str());

	nvar = reader.getNumCols();
	ncons = reader.getNumRows();

	//debugging checks
	for (int i = 0; i < nvar; i++) {
		assert(reader.getColLower()[i] <= reader.getColUpper()[i]);
	}
	for (int i = 0; i < ncons; i++) {
		assert(reader.getRowLower()[i] <= reader.getRowUpper()[i]);
	}

	// parse tim file
	ifstream ft(tim.c_str());
	ft.exceptions(ifstream::failbit | ifstream::badbit);
	string line;
	// skip first two lines
	getline(ft, line);
	assert(line.find("TIME") == 0);
	getline(ft, line);
	assert(line.find("PERIODS") == 0);
	
	getline(ft,line);
	string c, r, period;
	{
		istringstream is(line);
		is >> c;
		is >> r;
		is >> period;
		assert(c == reader.columnName(0));
		assert(r == reader.rowName(0));
	}

	getline(ft,line);
	{	
		istringstream is(line);
		is >> c;
		is >> r;
		is >> period;
	}

	ft.close();

	for (int i = 0; i < nvar; i++) {
		if (reader.isInteger(i)) {
			assert(reader.getColLower()[i] == 0.0);
			assert(reader.getColUpper()[i] == 1.0);
		}
	}	

	// find first column of second stage
	nvar1 = reader.columnIndex(c.c_str());
	assert(nvar1 >= 0 && nvar1 < nvar);

	nvar2 = nvar - nvar1;

	// find first row of second stage
	ncons1 = reader.rowIndex(r.c_str());
	assert(ncons1 >= 0 && ncons1 < ncons);
	ncons2 = ncons - ncons1;

	// load first-stage
	firstStageData.initialize(nvar1,ncons1);
	vector<int> cols(nvar1), rows(ncons1);
	for (int i = 0; i < nvar1; i++) cols[i] = i;
	for (int i = 0; i < ncons1; i++) rows[i] = i;
	
	firstStageData.mat = CoinPackedMatrix(*reader.getMatrixByCol(),
			ncons1, &rows[0], nvar1, &cols[0]);
	
	copy(reader.getColLower(),reader.getColLower()+nvar1,&firstStageData.collb[0]);
	copy(reader.getColUpper(),reader.getColUpper()+nvar1,&firstStageData.colub[0]);
	copy(reader.getObjCoefficients(),reader.getObjCoefficients()+nvar1,&firstStageData.obj[0]);
	
	copy(reader.getRowLower(),reader.getRowLower()+ncons1,&firstStageData.rowlb[0]);
	copy(reader.getRowUpper(),reader.getRowUpper()+ncons1,&firstStageData.rowub[0]);

	copySubset(reader, &CoinMpsIO::columnName, 0, nvar1, firstStageData.colname);
	copySubset(reader, &CoinMpsIO::rowName, 0, ncons1, firstStageData.rowname);
	copySubset(reader, &CoinMpsIO::isInteger, 0, nvar1, firstStageData.isColInteger);


	// load second-stage template
	secondStageTemplate.initialize(nvar2,ncons2);
	cols.resize(nvar2); rows.resize(ncons2);
	for (int i = 0; i < nvar2; i++) cols[i] = i+nvar1;
	for (int i = 0; i < ncons2; i++) rows[i] = i+ncons1;

	secondStageTemplate.mat = CoinPackedMatrix(*reader.getMatrixByCol(),
			ncons2, &rows[0], nvar2, &cols[0]);

	copy(reader.getColLower()+nvar1,reader.getColLower()+nvar,&secondStageTemplate.collb[0]);
	copy(reader.getColUpper()+nvar1,reader.getColUpper()+nvar,&secondStageTemplate.colub[0]);
	copy(reader.getObjCoefficients()+nvar1,reader.getObjCoefficients()+nvar,&secondStageTemplate.obj[0]);

	
	copy(reader.getRowLower()+ncons1,reader.getRowLower()+ncons,&secondStageTemplate.rowlb[0]);
	copy(reader.getRowUpper()+ncons1,reader.getRowUpper()+ncons,&secondStageTemplate.rowub[0]);

	copySubset(reader, &CoinMpsIO::columnName, nvar1, nvar, secondStageTemplate.colname);
	copySubset(reader, &CoinMpsIO::rowName, ncons1, ncons, secondStageTemplate.rowname);
	copySubset(reader, &CoinMpsIO::isInteger, nvar1, nvar, secondStageTemplate.isColInteger);

	continuousrecourse = true;
	for (int i = 0; i < nvar2; i++) {
		if (secondStageTemplate.isColInteger[i]) {
			continuousrecourse = false; break;
		}
	}

	cols.resize(nvar1);
	for (int i = 0; i < nvar1; i++) cols[i] = i;
	TmatTemplate = CoinPackedMatrix(*reader.getMatrixByCol(),
		ncons2, &rows[0], nvar1, &cols[0]);

	// pass through STO file
	ifstream fs(sto.c_str());
	fs.exceptions(ifstream::failbit | ifstream::badbit);

	getline(fs,line);
	assert(line.find("STOCH") == 0);
	getline(fs,line);
	assert(line.find("SCENARIOS") == 0);
	assert(line.find("DISCRETE") != string::npos);
	int nlines = -1;
	onlyboundsvary = true;

	getline(fs,line);
	while (line.find("ENDATA") == string::npos) {
		istringstream ss(line);
		if (line.find("SC") != string::npos) {
			if (nlines >= 0) {
				scenarioLens.push_back(nlines);
			}
			scenarioStarts.push_back(fs.tellg());
			nlines = 0;
			string s;
			double p;
			ss >> s; assert(s == "SC");
			ss >> s;
			ss >> s; assert(s.find("ROOT") != string::npos);
			ss >> p;
			probabilities.push_back(p);
		} else {
			nlines++;
			if (line.find("RHS") == string::npos) {
				ss >> c;
				ss >> r;
				if (reader.rowIndex(r.c_str()) != ncons) { // not objective row
					onlyboundsvary = false;
				}
			}
		}
		getline(fs,line);
	}
	fs.close();
	
	scenarioLens.push_back(nlines);
	assert(scenarioLens.size() == scenarioStarts.size());
	nscen = scenarioLens.size();
	Tmats.resize(nscen);
	scenarioData.resize(nscen);

	probabilitiesequal = true;
	double sum = 0;
	for (int i = 0; i < nscen; i++) {
		// can test for equality b/c they should be identical if
		// read from identical text
		if (probabilities[i] != probabilities[0]) {
			probabilitiesequal = false;
		}
		sum += probabilities[i];
	}
	assert(fabs(sum-1.0) < 1e-4);

}

vector<double> SMPSInput::getSecondStageColLB(int scen) { 
	cacheScenario(scen);
	return scenarioData.at(scen).collb;

}

vector<double> SMPSInput::getSecondStageColUB(int scen) { 
	cacheScenario(scen);
	return scenarioData.at(scen).colub;

}

vector<double> SMPSInput::getSecondStageObj(int scen) { 
	cacheScenario(scen);
	vector<double> obj = scenarioData.at(scen).obj;
	double scale = scenarioProbability(scen);
	for (unsigned i = 0; i < obj.size(); i++) obj[i] *= scale;
	return obj;

}

vector<double> SMPSInput::getSecondStageRowLB(int scen) { 
	cacheScenario(scen);
	
	return scenarioData.at(scen).rowlb;

}

vector<double> SMPSInput::getSecondStageRowUB(int scen) { 
	cacheScenario(scen);

	return scenarioData.at(scen).rowub;

}


vector<string> SMPSInput::getSecondStageColNames(int scen) {
	return scenarioData.at(scen).colname;
}

vector<string> SMPSInput::getSecondStageRowNames(int scen) {
	return scenarioData.at(scen).colname;
}


CoinPackedMatrix SMPSInput::getSecondStageConstraints(int scen) {
	cacheScenario(scen);
	return (onlyboundsvary ? secondStageTemplate.mat : scenarioData[scen].mat);
}


CoinPackedMatrix SMPSInput::getLinkingConstraints(int scen) {
	cacheScenario(scen);
	return (onlyboundsvary ? TmatTemplate : Tmats[scen]);
}

void SMPSInput::cacheScenario(int scen) {
	if (scenarioData[scen].isInitialized()) return;
	scenarioData[scen].initialize(nvar2, ncons2);

	// copy
	// don't need to copy matrices if only bounds vary
	scenarioData[scen] = secondStageTemplate;
	if (!onlyboundsvary) Tmats[scen] = TmatTemplate;

	ifstream fs(stofile.c_str());
	fs.exceptions(ifstream::failbit | ifstream::badbit);
	// seek to where scenario starts
	fs.seekg(scenarioStarts[scen]);
	
	string line, col, row;
	double val;
	for (int i = 0; i < scenarioLens[scen]; i++) {
		getline(fs, line);
		istringstream ss(line);
		ss >> col;
		ss >> row;
		ss >> val;
		
		int r = reader.rowIndex(row.c_str());
		int scenRow = r - ncons1;

		if (col.find("RHS") != string::npos) {
			assert(scenRow >= 0 && scenRow < ncons2);
			if (reader.getRowSense()[r] == 'L') {
				scenarioData[scen].rowub[scenRow] = val;
			} else if (reader.getRowSense()[r] == 'G') {
				scenarioData[scen].rowlb[scenRow] = val;
			} else if (reader.getRowSense()[r] == 'E') {
				assert(scenarioData[scen].rowlb[scenRow] == 
					scenarioData[scen].rowub[scenRow]);
				scenarioData[scen].rowlb[scenRow] = val;
				scenarioData[scen].rowub[scenRow] = val;
			} else {
				// not sure what input looks like when upper/lower bounds are
				// changed separately. don't have an example
				assert(0);
			}

			//printf("%d %d %f\n", scen, scenRow, val);

		} else {
			int c = reader.columnIndex(col.c_str());
			int scenCol = c - nvar1;
			if (r == ncons) { // objective row
				assert(scenCol >= 0 && scenCol < nvar2);
				scenarioData[scen].obj[scenCol] = val;
			} else {  
				assert(!onlyboundsvary);
				assert(scenRow >= 0 && scenRow < ncons2);
				if (c < nvar1) { // T matrix
					// SMPS assumes coefficient must
					// exist, but we don't check here
					Tmats[scen].modifyCoefficient(scenRow, c, val);
				} else { // W mat
					scenarioData[scen].mat.modifyCoefficient(scenRow, scenCol, val);
				}
			}
		}	


		
	}

	// make sure we're at the end
	getline(fs, line);
	assert(line.find("SC") != string::npos || line.find("ENDATA") != string::npos);
	

}
