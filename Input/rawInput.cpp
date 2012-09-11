#include "rawInput.hpp"

#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

rawInput::rawInput(const string &datarootname, int overrideScenarioNumber, MPI_Comm comm) : datarootname(datarootname) {

  unsigned filelen;
  char *filedata;
  MPI_Comm_rank(comm,&mype_);

  if (mype_ == 0) {

    string f0name = datarootname + "0";
    ifstream f0(f0name.c_str());
    assert(f0.is_open());
    
    string data((istreambuf_iterator<char>(f0)),istreambuf_iterator<char>()); // read into string
    f0.close();
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
  
  parseZeroData(data, overrideScenarioNumber);

}

rawInput::rawInput(const std::string &datarootname, 
		   const std::string& zerodata, 
		   int overrideScenarioNumber /*= 0*/, 
		   MPI_Comm comm /*= MPI_COMM_SELF*/)
  : datarootname(datarootname)
{
  parseZeroData(zerodata, overrideScenarioNumber);
}

void rawInput::parseZeroData(const std::string &zerodata, int overrideScenarioNumber)
{
	istringstream f1(zerodata);
	f1.exceptions(ifstream::failbit | ifstream::badbit);
	
	
	f1 >> nFirstStageVars_;
	f1 >> nFirstStageCons_;
	f1 >> nScenariosTrue;

	assert(nScenariosTrue == 1);
	
	string line;
	getline(f1,line);
	getline(f1,line);
	size_t loc = line.find("Only Bounds Vary");
	assert(loc == 0);
	f1 >> nSecondStageVars_;
	f1 >> nSecondStageCons_;

	if (overrideScenarioNumber) {
		nScenarios_ = overrideScenarioNumber;
	} else {
		nScenarios_ = nScenariosTrue;
	}

	localData.resize(nScenarios_);

	firstStageData.initialize(nFirstStageVars_, nFirstStageCons_);
	for (int i = 0; i < nFirstStageVars_; i++) {
		f1 >> firstStageData.colnames[i];
		f1 >> firstStageData.collb[i];
		f1 >> firstStageData.colub[i];
		f1 >> firstStageData.obj[i];
		//printf("%s, %f, %f, %f\n",firstStageData.colnames[i].c_str(),firstStageData.collb[i],firstStageData.colub[i],firstStageData.obj[i]);
	}
	for (int i = 0; i < nFirstStageCons_; i++) {
		f1 >> firstStageData.rownames[i];
		f1 >> firstStageData.rowlb[i];
		f1 >> firstStageData.rowub[i];
	}


	getline(f1,line);
	getline(f1,line);
	loc = line.find("A matrix");
	assert(loc == 0);

	CoinBigIndex nnz;
	vector<CoinBigIndex> starts;
	vector<int> rowIdx;
	vector<double> elts;
	
	f1 >> nnz;
	starts.resize(nFirstStageVars_+1);
	rowIdx.resize(nnz);
	elts.resize(nnz);
	for (int i = 0; i <= nFirstStageVars_; i++) {
		f1 >> starts[i];
	}
	for (int i = 0; i < nnz; i++) {
		f1 >> rowIdx[i];
	}
	for (int i = 0; i < nnz; i++) {
		f1 >> elts[i];
	}
	Amat.copyOf(true, nFirstStageCons_, nFirstStageVars_, nnz, &elts[0], &rowIdx[0], &starts[0], 0);

	getline(f1,line);
	getline(f1,line);
	loc = line.find("W matrix");
	assert(loc == 0);

	f1 >> nnz;
	starts.resize(nSecondStageVars_+1);
	rowIdx.resize(nnz);
	elts.resize(nnz);
	for (int i = 0; i <= nSecondStageVars_; i++) {
		f1 >> starts[i];
	}
	for (int i = 0; i < nnz; i++) {
		f1 >> rowIdx[i];
	}
	for (int i = 0; i < nnz; i++) {
		f1 >> elts[i];
	}
	Wmat.copyOf(true, nSecondStageCons_, nSecondStageVars_, nnz, &elts[0], &rowIdx[0], &starts[0], 0);

	getline(f1,line);
	getline(f1,line);
	loc = line.find("T matrix");
	assert(loc == 0);

	f1 >> nnz;
	starts.resize(nFirstStageVars_+1);
	rowIdx.resize(nnz);
	elts.resize(nnz);
	for (int i = 0; i <= nFirstStageVars_; i++) {
		f1 >> starts[i];
	}
	for (int i = 0; i < nnz; i++) {
		f1 >> rowIdx[i];
	}
	for (int i = 0; i < nnz; i++) {
		f1 >> elts[i];
	}
	Tmat.copyOf(true, nSecondStageCons_, nFirstStageVars_, nnz, &elts[0], &rowIdx[0], &starts[0], 0);

}

void rawInput::scenData::initialize(int nvar, int ncons) {
	collb.resize(nvar);
	colub.resize(nvar);
	obj.resize(nvar);
	rowlb.resize(ncons);
	rowub.resize(ncons);
	rownames.resize(ncons);
	colnames.resize(nvar);
	didLoad = true;
}


void rawInput::loadLocalScenData(int scen) {
	if (localData[scen].didLoad) return;

	cout << "Proc " << mype_ << " rawInput::loadLocalScenData for scenario " << scen << endl;
	

	localData[scen].initialize(nSecondStageVars_,nSecondStageCons_);

	stringstream fname;
	fname << datarootname << scen+1;
	ifstream f(fname.str().c_str());
	if (!f.is_open()) {
		cerr << "Unable to open" << fname.str() << " from process " << mype_ << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	f.exceptions(ifstream::failbit | ifstream::badbit);

	int nvar, ncons;
	f >> nvar;
	f >> ncons;
	assert(nvar == nSecondStageVars_ && ncons == nSecondStageCons_);
	string line;
	getline(f,line);
	getline(f,line);
	size_t loc = line.find("Only Bounds Vary");
	assert(loc == 0);
	double scale = 1.0;
	if (nScenariosTrue != nScenarios_) {
		scale = static_cast<double>(nScenariosTrue)/static_cast<double>(nScenarios_);
		//printf("Rescaling scenario objectives by %f\n",scale);
	}

	for (int i = 0; i < nvar; i++) {
		f >> localData[scen].colnames[i];
		f >> localData[scen].collb[i];
		f >> localData[scen].colub[i];
		f >> localData[scen].obj[i];
		
		// NOTE WE'RE ASSUMING EQUALLY-LIKELY SCENARIOS
		// DON'T OVERRIDE SCENARIO COUNT IF THAT'S NOT THE CASE
		localData[scen].obj[i] *= scale;
	}
	for (int i = 0; i < ncons; i++) {
		f >> localData[scen].rownames[i];
		f >> localData[scen].rowlb[i];
		f >> localData[scen].rowub[i];
	}

}

vector<double> rawInput::getSecondStageColLB(int scen) {
	loadLocalScenData(scen);
	return localData[scen].collb;
}

vector<double> rawInput::getSecondStageColUB(int scen) {
	loadLocalScenData(scen);
	return localData[scen].colub;
}

vector<double> rawInput::getSecondStageObj(int scen) {
	loadLocalScenData(scen);
	return localData[scen].obj;
}

vector<string> rawInput::getSecondStageColNames(int scen) {
	loadLocalScenData(scen);
	return localData[scen].colnames;
}

vector<double> rawInput::getSecondStageRowUB(int scen) {
	loadLocalScenData(scen);
	return localData[scen].rowub;
}

vector<double> rawInput::getSecondStageRowLB(int scen) {
	loadLocalScenData(scen);
	return localData[scen].rowlb;
}

vector<string> rawInput::getSecondStageRowNames(int scen) {
	loadLocalScenData(scen);
	return localData[scen].rownames;
}



